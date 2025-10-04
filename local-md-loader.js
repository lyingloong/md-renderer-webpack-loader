const { readFileSync } = require('fs');
const { dirname, resolve } = require('path');
const { pathToFileURL } = require('url');

const submoduleEntry = resolve(__dirname, './md-renderer-react/dist/index.js');

const loadMdRenderer = async () => {
  const fileUrl = pathToFileURL(submoduleEntry).href;
  const submoduleExports = await import(fileUrl);

  const parser = submoduleExports.parser;
  if (!parser) throw new Error('未找到 parser 实例');

  const mdRenderer = submoduleExports.mdRenderer || submoduleExports.default?.mdRenderer || submoduleExports.default;
  if (!mdRenderer) throw new Error('未找到 mdRenderer 实例');
  return { mdRenderer, parser };
};

const getImageVarName = (path) =>
  `img_${Buffer.from(path).toString('base64').replace(/[^a-zA-Z0-9_]/g, '_')}`;

const isLocalPath = (path) => !/^https?:\/\//.test(path) && !path.startsWith('/');

// 从 AST 节点中提取纯文本
const extractTextFromAST = (astNode) => {
  if (typeof astNode === 'string') return astNode;
  if (!astNode) return '';

  // 处理数组类型的 content
  if (Array.isArray(astNode)) {
    return astNode.map(extractTextFromAST).join('');
  }

  // 处理对象类型的 content
  if (astNode.content) {
    return extractTextFromAST(astNode.content);
  }

  // 处理 plain 类型的节点
  if (astNode.type === 'plain' && astNode.content) {
    return astNode.content;
  }

  return '';
};

// 生成标题 ID
const generateHeadingId = (text, index) => {
  const cleanText = text
    .toLowerCase()
    .replace(/[^\w\s-]/g, '') // 移除特殊字符
    .replace(/\s+/g, '-') // 空格替换为连字符
    .replace(/-+/g, '-') // 多个连字符合并为一个
    .replace(/^-|-$/g, ''); // 移除首尾连字符

  return cleanText || `heading-${index}`;
};


module.exports = async function (source) {
  /*
  optional style config
  */
  const options = this.getOptions ? this.getOptions() : this.query || {};
  let cssImport = '';
  if (options.cssPath) {
      cssImport = `import styles from '${options.cssPath}';`;
      console.log(`[local-md-loader] import css: '${options.cssPath}'`);
  }

  const callback = this.async();
  try {
    const { mdRenderer, parser } = await loadMdRenderer();

    const mdFilePath = this.resourcePath;
    this.addDependency(mdFilePath);
    const mdContent = readFileSync(mdFilePath, 'utf8');

    // 获取 AST 数据用于提取标题
    const ast = await parser.parseMarkdownToAST(mdContent);
    const headings = [];

    // 从 AST 中提取标题信息
    ast.forEach((section, index) => {
      if (section.type === 'section' && section.title) {
        // 生成标题 ID
        const titleText = extractTextFromAST(section.title);
        const id = generateHeadingId(titleText, index);

        headings.push({
          id: id,
          text: titleText,
          type: 'heading',
        });
      }
    });

    const codeResult = await mdRenderer.getRenderedCode(mdContent);
    if (!codeResult.success) throw new Error(`[local-md-loader] error: ${codeResult.data}`);

    let html = codeResult.data;
    const imageImports = new Map();

    // 为标题添加 ID 和 CSS 类名
    let headingIndex = 0;
    html = html.replace(/<h2([^>]*)>/g, (match, attrs) => {
      const id = headingIndex < headings.length ? headings[headingIndex].id : '';
      headingIndex++;
      return `<h2${attrs} id="${id}" class="subheading">`;
    });
    html = html.replace(/<p([^>]*)>/g, '<p$1 class="paragraph">');

    // 匹配 <img src="...">
    html = html.replace(/<img\s+[^>]*src="([^"]+)"[^>]*>/g, (match, srcPath) => {
      if (isLocalPath(srcPath)) {
        const fullPath = srcPath.replace(/\\/g, '/');
        const varName = getImageVarName(fullPath);
        if (!imageImports.has(srcPath)) imageImports.set(srcPath, { varName, fullPath });
        return match.replace(srcPath, `\${${varName}}`);
      }
      console.log("[local-md-loader] img match:", match)
      return match;
    });

    // 匹配 <Figure path="..."> 可选
    html = html.replace(/<Figure\s+[^>]*path="([^"]+)"[^>]*\/?>/g, (match, path) => {
      if (isLocalPath(path)) {
        const fullPath = srcPath.replace(/\\/g, '/');
        const varName = getImageVarName(fullPath);
        if (!imageImports.has(path)) imageImports.set(path, { varName, fullPath });
        return match.replace(path, `\${${varName}}`);
      }
      console.log("[local-md-loader] img match:", match)
      return match;
    });

    const importStatements = Array.from(imageImports.values())
      .map(({ varName, fullPath }) => `import ${varName} from '${fullPath}';`)
      .join('\n');

    const finalCode = `
      import React from 'react';
      import 'katex/dist/katex.min.css';
      import 'highlight.js/styles/github.css';
      ${cssImport}

      ${importStatements}
      
      // 导出标题数据
      export const headings = ${JSON.stringify(headings)};
      
      export default () => (
        <div dangerouslySetInnerHTML={{ __html: \`${html}\` }} />
      );
    `;

    callback(null, finalCode);
  } catch (err) {
    console.error('Loader error:', err.stack);
    this.emitError(`Loader 执行失败：${err.message}`);
    callback(null, `
      import React from 'react';
      export default () => (
        <div style={{ 
          color: 'red', 
          padding: '16px', 
          border: '1px solid #ffccc7', 
          borderRadius: '4px', 
          background: '#fff5f5',
          fontFamily: 'sans-serif'
        }}>
          <h3>Markdown 加载失败</h3>
          <pre style={{ whiteSpace: 'pre-wrap' }}>${err.message.replace(/`/g, '\\`').replace(/\$/g, '\\$')}</pre>
        </div>
      );
    `);
  }
};
