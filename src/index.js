import React from 'react';
import ReactDOM from 'react-dom/client';
import DemoMd from './md/demo.md';
import DemoListMd from './md/demo_list.md';

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <DemoListMd />
  </React.StrictMode>
);