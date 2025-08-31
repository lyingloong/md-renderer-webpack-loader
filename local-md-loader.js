const { readFileSync } = require('fs');
const { dirname, resolve } = require('path');

const submoduleEntry = resolve(__dirname, './md-renderer-react/dist/index.js');

const loadMdRenderer = async () => {
  console.log('\n===== 开始加载子模块 =====');
  const submoduleExports = await import(`file://${submoduleEntry}`);
  const mdRenderer = submoduleExports.mdRenderer || submoduleExports.default?.mdRenderer;
  if (!mdRenderer) throw new Error('未找到 mdRenderer 实例');
  return mdRenderer;
};

const getImageVarName = (path) =>
  `img_${Buffer.from(path).toString('base64').replace(/[^a-zA-Z0-9_]/g, '_')}`;

const isLocalPath = (path) => !/^https?:\/\//.test(path) && !path.startsWith('/');

module.exports = async function (source) {
  const callback = this.async();
  try {
    const mdRenderer = await loadMdRenderer();

    const mdFilePath = this.resourcePath;
    this.addDependency(mdFilePath);
    const mdContent = readFileSync(mdFilePath, 'utf8');

    const codeResult = await mdRenderer.getRenderedCode(mdContent);
    if (!codeResult.success) throw new Error(`子模块内部错误: ${codeResult.data}`);

    let html = codeResult.data;
    const imageImports = new Map();

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
      ${importStatements}
      export default () => (
        <div className="react-rendered-content" dangerouslySetInnerHTML={{ __html: \`${html}\` }} />
      );
    `;

    callback(null, finalCode);
  } catch (err) {
    console.error('Loader 处理出错:', err.stack);
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
