# webpack-igemd-loader

A custom Webpack loader to render Markdown (`.md`) files as React components with support for images, KaTeX math, and syntax highlighting. Ideal for projects that include this repository as a submodule.

## Features

- Automatically converts Markdown files into React components.
- Handles local images and generates proper imports.
- Supports KaTeX for mathematical formulas.
- Syntax highlighting via `highlight.js`.
- Compatible with Webpack and modern React setups.

## Usage

### Webpack Configuration

In your main project, configure Webpack to use `local-md-loader.js` for `.md` files. Example:

```
const path = require('path');

module.exports = {
  module: {
    rules: [
      {
        test: /\.md$/,
        use: [
          {
            loader: 'babel-loader',
            options: { presets: ['@babel/preset-react'], babelrc: false, configFile: false },
          },
          {
            loader: path.resolve(__dirname, 'local-md-loader.js'),
          },
        ],
      },
      // other rules...
    ],
  },
  resolve: { extensions: ['.js', '.jsx', '.md'] },
};
```

### Markdown in React

Once configured, you can import Markdown files like this:

```
import MyMarkdown from './docs/example.md';

export default function App() {
  return (
    <div>
      <MyMarkdown />
    </div>
  );
}
```

The loader will automatically:

- Render Markdown to HTML.
- Import local images and update `<img>` and `<Figure>` paths.
- Apply KaTeX and syntax highlighting styles.

## Project Structure

```
webpack-igemd-loader/
├── md-renderer-react/   # Submodule providing the Markdown rendering engine
├── src/                 # Demo / test entry for trying out the loader
├── .gitignore
├── .gitmodules          # Defines the `md-renderer-react` submodule
├── local-md-loader.js   # Custom Webpack loader implementation
├── package.json         # Project scripts and dependencies
├── package-lock.json
├── readme.md            # Documentation
└── webpack.config.js    # Example Webpack configuration for demo
```

## Notes

- **Submodule Path:** The loader dynamically imports the submodule bundle at `md-renderer-react/dist/index.js`. Make sure to run  `npm build:submodule` in this repo or in the submodule before using the loader.
- **Local Images:** Only local image paths are processed. URLs starting with `http(s)://` or `/` are not transformed.
- **Error Handling:** If the loader encounters errors in rendering, it outputs a red error box in the rendered React component for easy debugging.

## Scripts

```
{
  "scripts": {
    "prepare": "cd md-renderer-react && npm install",
    "build:submodule": "cd md-renderer-react && npm run build",
    "prestart": "npm run build:submodule",
    "start": "webpack serve --config webpack.config.js",
    "build": "webpack --config webpack.config.js"
  }
}
```

- `prepare`: Install submodule dependencies.
- `build:submodule`: Build the Markdown renderer submodule.
- `prestart`: Ensures the submodule is built before running the dev server.
- `start`: Start a Webpack dev server with hot reload.
- `build`: Build the main project bundle.

## Dependencies

- React & React-DOM
- Babel (`@babel/preset-env`, `@babel/preset-react`, `babel-loader`)
- CSS handling (`style-loader`, `css-loader`)
- Highlight.js for syntax highlighting
- KaTeX for math rendering

## License

MIT