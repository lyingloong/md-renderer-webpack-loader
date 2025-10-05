const path = require('path');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const CopyPlugin = require('copy-webpack-plugin');

module.exports = {
  mode: 'development',
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
  },
  plugins: [
    new HtmlWebpackPlugin({
      filename: 'index.html',
      template: './src/index.html',
      title: 'Markdown 渲染Demo',
    }),
    new CopyPlugin({
      patterns: [
        {
          from: path.resolve(__dirname, 'src/docs'), // md 源目录
          to: 'docs',                               // 打包到 dist/docs
          globOptions: {
            ignore: ['**/*.DS_Store'],              // 忽略 macOS 文件等
          },
        },
      ],
    }),
  ],
  module: {
    rules: [
      // 规则1：处理 .md 文件 → 自定义 Loader + babel-loader
      {
        test: /\.md$/,
        use: [
          {
            loader: 'babel-loader',
            options: {
              presets: ['@babel/preset-react'],
              babelrc: false,
              configFile: false,
            },
          },
          {
            loader: path.resolve(__dirname, 'local-md-loader.js'),
            options: {
              cssPath: '../styles/style.module.scss'
            }
          },
        ],
      },
      // 规则2：处理图片资源
      {
        test: /\.(png|jpe?g|gif|svg)$/i,
        type: 'asset/resource',
        generator: {
          filename: 'images/[name].[hash:6][ext]',
        },
      },
      // 规则3：处理主项目 JS/JSX
      {
        test: /\.(js|jsx)$/,
        include: path.resolve(__dirname, "src"),
        exclude: /node_modules/,
        use: {
          loader: 'babel-loader',
          options: {
            presets: ['@babel/preset-react', '@babel/preset-env'],
          },
        },
      },
      // 规则4：处理 CSS 文件（包括 KaTeX CSS）
      {
        test: /\.css$/i,
        use: ['style-loader', 'css-loader'],
      },
      // 规则5：处理 SCSS 文件（启用 CSS Modules）
      {
        test: /\.module\.s[ac]ss$/i,   // 只处理 .module.scss
        use: [
          'style-loader',
          {
            loader: 'css-loader',
            options: {
              modules: {
                localIdentName: '[local]_[hash:6]', // 可选：自定义 class 名规则
              },
            },
          },
          'sass-loader',
        ],
      },
      // 处理普通 SCSS 文件，不启用 modules
      {
        test: /\.s[ac]ss$/i,
        exclude: /\.module\.s[ac]ss$/i,
        use: ['style-loader', 'css-loader', 'sass-loader'],
      },
      // 规则6：处理 KaTeX 字体文件
      {
        test: /\.(woff2?|ttf|eot|svg)$/,
        type: 'asset/resource',
        generator: {
          filename: 'fonts/[name].[hash:6][ext]',
        },
      },
    ],
  },
  resolve: {
    extensions: ['.js', '.jsx', '.md'],
    alias: {
      react: path.resolve(__dirname, "node_modules/react"),
      "react-dom": path.resolve(__dirname, "node_modules/react-dom"),
    },
  },
  devServer: {
    static: path.resolve(__dirname, 'dist'),
    port: 3220,
    hot: true,
    open: false,
    // client: {
    //   webSocketURL: 'ws://localhost:8888/ws',
    // },
    webSocketServer: 'ws',
  },
};
