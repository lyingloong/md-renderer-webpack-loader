const path = require('path');

module.exports = {
  mode: 'development',
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
  },
  module: {
    rules: [
      // 规则1：处理 .md 文件 → 先自定义 Loader 输出 JSX，再 babel-loader 转译 JSX
      {
        test: /\.md$/,
        use: [
          {
            // 第2步执行：转译自定义 Loader 输出的 JSX（必须在自定义 Loader 之前，因为 Loader 从后往前执行）
            loader: 'babel-loader',
            options: {
              presets: ['@babel/preset-react'], // 仅需解析 JSX，无需额外预设
              babelrc: false, // 禁用项目根目录的 babel 配置，避免冲突
              configFile: false, // 禁用单独的 babel 配置文件
            },
          },
          {
            // 第1步执行：自定义 Loader 处理 Markdown，输出含 JSX 的代码
            loader: path.resolve(__dirname, 'local-md-loader.js'),
          },
        ],
      },
      // 规则2：处理图片资源（无问题，保留）
      {
        test: /\.(png|jpe?g|gif|svg)$/i,
        type: 'asset/resource',
        generator: {
          filename: 'images/[name].[hash:6][ext]', // 带 hash 防缓存，合理
        },
      },
      // 规则3：处理主项目自己的 JS/JSX（无问题，保留）
      {
        test: /\.(js|jsx)$/,
        include: path.resolve(__dirname, "src"), // 仅处理主项目 src 目录，合理
        exclude: /node_modules/, // 排除 node_modules，合理
        use: {
          loader: 'babel-loader',
          options: {
            presets: ['@babel/preset-react', '@babel/preset-env'], // 同时处理 JSX 和 ES 语法，合理
          },
        },
      },
    ],
  },
  resolve: {
    extensions: ['.js', '.jsx', '.md'], // 支持省略 .md/.jsx 后缀，合理
    alias: {
      // React 别名：避免主项目与子模块 React 版本冲突，可选但推荐保留
      "react": path.resolve(__dirname, "node_modules/react"),
      "react-dom": path.resolve(__dirname, "node_modules/react-dom")
    }
  },
  devServer: {
    static: path.resolve(__dirname, 'dist'), // 开发服务器静态目录，合理
    port: 3200, // 端口号，合理
    hot: true, // 热更新，合理
  },
};