const path = require('path')
const MiniCssExtractPlugin = require('mini-css-extract-plugin')


module.exports = {
    entry: {
        main: './src/js/main.js',
        index: './src/js/index.js'
    },
    output: {
        filename: 'js/[name].bundle.js',
        path: __dirname + '/main/static/'
    },
    module: {
        rules: [
            {
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
            },
            {
                test: /\.js$/,
                exclude: /node_modules/,
                use: {
                    loader: 'babel-loader'
                }
            },
            {
                test: /\.vue$/,
                loader: 'vue-loader'
            }
        ]
    },
    plugins: [
        new MiniCssExtractPlugin({
            filename: 'css/bundle.css'
        }),
    ],
    resolve: {
        alias: {
            'vue$': 'vue/dist/vue.esm.js'
        },
        extensions: ['*', '.js', '.vue', '.json']
    }

}