const path = require('path')
const MiniCssExtractPlugin = require('mini-css-extract-plugin')


module.exports = {
    entry: {
        index: './src/js/index.js'
    },
    output: {
        filename: 'js/[name].bundle.js',
        path: __dirname + '/main/static/'
    },
    module: {
        rules: [{
            test: /\.scss$/,
            use: [
                    MiniCssExtractPlugin.loader,
                    {
                        loader: 'css-loader'
                    },
                    {
                        loader: 'sass-loader',
                        options: {
                            sourceMap: true,
                            // options...
                        }
                    }
                ]
        }]
    },
    plugins: [
        new MiniCssExtractPlugin({
            filename: 'css/bundle.css'
        }),
    ]
}