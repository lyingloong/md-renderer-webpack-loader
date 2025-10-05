import React from 'react';
import ReactDOM from 'react-dom/client';
import DemoMd, { getAST } from './docs/demo.md';

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <DemoMd />
  </React.StrictMode>
);