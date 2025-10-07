import React from 'react';
import ReactDOM from 'react-dom/client';
import DemoMd from './md/demo.md';
import DemoModelMd from './md/demo_model.md';

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <DemoModelMd />
  </React.StrictMode>
);