const { readFileSync } = require('fs');
const { resolve } = require('path');
const { pathToFileURL } = require('url');

// md-renderer-react submodule path
const submoduleEntry = resolve(__dirname, './md-renderer-react/dist/index.js');
console.log('[local-md-loader] 执行环境:', typeof window === 'undefined' ? 'Node (SSR)' : 'Browser');

const isNonEmptyString = (v) => typeof v === 'string' && v.trim().length > 0;
const isNonEmptyArray  = (arr) => Array.isArray(arr) && arr.length > 0;

module.exports = async function () {
  if (this.resourcePath) {
    console.log('[local-md-loader] processing:', this.resourcePath);
  } else {
    console.warn('[local-md-loader] Warning: this.resourcePath is undefined');
  }
  const options = this.getOptions ? this.getOptions() : this.query || {};

  let cssImport = '';
  if (isNonEmptyString(options.cssPath)) {
      cssImport = `import * as styles from '${options.cssPath}';`;
      console.log(`[local-md-loader] import css: '${options.cssPath}'`);
  }

  const callback = this.async();
  try {
    const mdFileName = this.resourcePath ? this.resourcePath.split(/[\\/]/).pop() : 'unknown.md';
    const mdUrl = `static/md/${mdFileName}`; // 浏览器可访问 URL
    console.log(`[local-md-loader] mdUrl: '${mdUrl}'`);

    const finalCode = `
      import React, { useEffect, useState } from 'react';
      import 'katex/dist/katex.min.css';
      import 'highlight.js/styles/github.css';
      ${cssImport}
      import { MDRenderer } from '${pathToFileURL(submoduleEntry).href}';

      const renderer = new MDRenderer({ styles });

      export default function MarkdownComponent() {
        const [element, setElement] = useState(<div>loading...</div>);

        useEffect(() => {
          let canceled = false;

          (async () => {
            try {
              console.log("[frontend md-loader] fetching md from:", '${mdUrl}');
              const res = await fetch('${mdUrl}');
              const mdContent = await res.text();
              console.log("[frontend md-loader] fetched md content:", mdContent);

              const el = await renderer.renderToElement(mdContent);
              if (canceled) return;

              if (!el) {
                console.error("[md-loader] renderToElement 返回 null 或 undefined");
                setElement(<div style={{ color: 'red' }}>Markdown 渲染失败（空元素）</div>);
                return;
              }

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
