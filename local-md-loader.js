const { readFileSync } = require('fs');
const { resolve } = require('path');
const { pathToFileURL } = require('url');

// md-renderer-react submodule path
const submoduleEntry = resolve(__dirname, './md-renderer-react/dist/index.js');

const getImageVarName = (path) =>
  `img_${Buffer.from(path).toString('base64').replace(/[^a-zA-Z0-9_]/g, '_')}`;

const isLocalPath = (path) => !/^https?:\/\//.test(path) && !path.startsWith('/');


module.exports = async function () {
  const options = this.getOptions ? this.getOptions() : this.query || {};
  let cssImport = '';
  if (options.cssPath) {
      cssImport = `import styles from '${options.cssPath}';`;
      console.log(`[local-md-loader] import css: '${options.cssPath}'`);
  }

  const callback = this.async();
  try {
    const mdFilePath = this.resourcePath;
    this.addDependency(mdFilePath);
    const mdContent = readFileSync(mdFilePath, 'utf8');

    const finalCode = `
      import React, { useEffect, useState } from 'react';
      import 'katex/dist/katex.min.css';
      import 'highlight.js/styles/github.css';
      ${cssImport}
      import { mdRenderer } from '${pathToFileURL(submoduleEntry).href}';

      let __ast = null;

      export default function MarkdownComponent() {
        const [element, setElement] = useState(<div>loading...</div>);

        useEffect(() => {
          let canceled = false;

          (async () => {
            try {
              const el = await mdRenderer.renderToElement(\`${mdContent.replace(/`/g, '\\`')}\`);
              if (canceled) return;

              if (!el) {
                console.error("[md-loader] renderToElement 返回 null 或 undefined");
                setElement(<div style={{ color: 'red' }}>Markdown 渲染失败（空元素）</div>);
                return;
              }

              __ast = mdRenderer.lastAST || null;
              setElement(el);
            } catch (err) {
              console.error("[md-loader] 渲染失败:", err);
              setElement(<div style={{ color: 'red' }}>Markdown 渲染失败: {err.message}</div>);
            }
          })();

          return () => {
            canceled = true;
          };
        }, []);

        return element;
      }

      export const getAST = () => __ast;
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
