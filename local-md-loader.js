const { readFileSync } = require('fs');
const fs = require('fs');
const { dirname, resolve } = require('path');
const babelParser = require('@babel/parser');
const traverse = require('@babel/traverse').default;
const generate = require('@babel/generator').default;
const t = require('@babel/types');

// 1. 子模块预编译产物入口（确保路径指向子模块 dist/index.js）
const submoduleEntry = resolve(__dirname, './md-renderer-react/dist/index.js');

// 2. 简化子模块加载逻辑：直接用 Node 原生 require（去掉冗余的 createCjsContext）

// 新 ESM 动态加载方式（核心修改）
const loadMdRenderer = async () => {
  try {
    const submoduleEntry = resolve(__dirname, './md-renderer-react/dist/index.js');
    if (!fs.existsSync(submoduleEntry)) {
      throw new Error(
        `子模块编译产物不存在！请先执行 "npm run build:submodule" 编译子模块\n产物路径：${submoduleEntry}`
      );
    }

    // 关键：用 ESM 动态 import() 加载子模块（而非 require）
    // 注意：Node.js 加载本地 ESM 文件需用 file:// 协议或绝对路径（不同 Node 版本可能有差异）
    const submoduleExports = await import(`file://${submoduleEntry}`);
    
    // 验证子模块导出（ESM 导出的 default 可能在 .default 上，需适配）
    const mdRenderer = submoduleExports.mdRenderer || submoduleExports.default?.mdRenderer;
    if (!mdRenderer) {
      throw new Error(
        `子模块未导出 "mdRenderer" 实例（当前导出内容：${Object.keys(submoduleExports).join(', ')}）`
      );
    }

    if (typeof mdRenderer.getRenderedCode !== 'function') {
      throw new Error(
        `mdRenderer 实例缺少 "getRenderedCode" 方法（当前实例方法：${Object.keys(mdRenderer).filter(key => typeof mdRenderer[key] === 'function').join(', ')}）`
      );
    }

    return mdRenderer;
  } catch (err) {
    throw new Error(`加载子模块失败：${err.message}\n子模块入口：${submoduleEntry}`);
  }
};

// 3. 图片路径处理逻辑（完全保留原有功能）
const getImageVarName = (path) => {
  return `img_${Buffer.from(path).toString('base64').replace(/[^a-zA-Z0-9_]/g, '_')}`;
};

const isLocalPath = (path) => {
  return !/^https?:\/\//.test(path) && !path.startsWith('/');
};

// 4. Loader 主逻辑（核心修改：方法名、参数、异步 await）
module.exports = async function (source) {
  const callback = this.async(); // 标记为异步 Loader，必须调用 callback

  try {
    // 加载子模块的 mdRenderer 实例
    const mdRenderer = await loadMdRenderer();
    const mdFilePath = this.resourcePath;
    this.addDependency(mdFilePath); // 告诉 Webpack：md 文件变化时重新构建

    // 读取 Markdown 内容
    const mdContent = readFileSync(mdFilePath, 'utf8');
    const mdDirectory = dirname(mdFilePath);

    // 关键修改1：调用子模块实际方法 getRenderedCode（而非 getReactComponentCode）
    // 关键修改2：子模块方法仅接收 mdContent 一个参数，去掉多余的 options
    // 关键修改3：方法是 async 的，必须加 await 才能拿到结果（否则拿到 Promise）
    const codeResult = await mdRenderer.getRenderedCode(mdContent);

    // 验证子模块返回结果格式
    if (typeof codeResult !== 'object' || codeResult.success === undefined) {
      throw new Error(`子模块方法返回格式异常，预期 { success: boolean, data: string }，实际：${JSON.stringify(codeResult)}`);
    }

    if (!codeResult.success) {
      throw new Error(`React 组件生成失败（子模块内部错误）：${codeResult.data}`);
    }

    // 处理图片路径（原有逻辑不变）
    const ast = babelParser.parse(codeResult.data, {
      sourceType: 'module', // 解析为 ES 模块
      plugins: ['jsx'],     // 支持解析 JSX 语法
    });

    const imageImports = new Map();
    traverse(ast, {
      JSXElement(path) {
        // 只处理 Figure 标签的 path 属性（原有逻辑）
        if (
          path.node.openingElement.name.type === 'JSXIdentifier' &&
          path.node.openingElement.name.name === 'Figure'
        ) {
          const pathAttr = path.node.openingElement.attributes.find(
            (attr) => attr.name && attr.name.name === 'path'
          );

          // 处理本地图片路径（转为 Webpack 可识别的 import）
          if (pathAttr?.value?.value && isLocalPath(pathAttr.value.value)) {
            const fullPath = resolve(mdDirectory, pathAttr.value.value);
            const varName = getImageVarName(fullPath);
            imageImports.set(pathAttr.value.value, { varName, fullPath });
            // 将路径字符串替换为 import 后的变量（如 <Figure path={img_xxx} />）
            pathAttr.value = t.jsxExpressionContainer(t.identifier(varName));
          }
        }
      },
    });

    // 给 AST 头部添加图片的 import 语句（原有逻辑不变）
    const importStatements = Array.from(imageImports.values()).map(
      ({ varName, fullPath }) =>
        t.importDeclaration(
          [t.importDefaultSpecifier(t.identifier(varName))], // import img_xxx from 'xxx'
          t.stringLiteral(fullPath) // 图片的绝对路径
        )
    );
    ast.program.body.unshift(...importStatements);

    // 生成最终可执行的 JS 代码
    const { code } = generate(ast);
    callback(null, code); // 成功：传递生成的代码

  } catch (err) {
    // 错误处理：生成错误提示组件（JSX 会被后续 babel-loader 转译）
    this.emitError(`Loader 执行失败：${err.message}`); // 告诉 Webpack 这是一个错误
    callback(null, `
      import React from 'react';
      // 错误提示组件（样式保留原有逻辑）
      export default () => (
        <div style={{ 
          color: 'red', 
          padding: '16px', 
          border: '1px solid #ffccc7', 
          borderRadius: '4px', 
          background: '#fff5f5',
          fontFamily: 'sans-serif'
        }}>
          <h3 style={{ margin: '0 0 8px', fontSize: '16px' }}>Markdown 加载失败</h3>
          <pre style={{ margin: '0', whiteSpace: 'pre-wrap', wordBreak: 'break-all' }}>
            ${err.message.replace(/`/g, '\\`').replace(/\$/g, '\\$')}
          </pre>
        </div>
      );
    `);
  }
};