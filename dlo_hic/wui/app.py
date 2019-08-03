from flask import Flask, g
from flask_login import LoginManager
from flask_bootstrap import Bootstrap
from flask_wtf.csrf import CSRFProtect


def gen_secret_key():
    import string
    from random import choice
    def random_string(length=20):
        letters = string.ascii_lowercase
        return ''.join([choice(letters) for _ in range(length)])
    return random_string()


def create_app(password):
    app = Flask("dlohic")

    from dlo_hic.wui.main.views import main
    app.register_blueprint(main)
    login = LoginManager(app)
    bootstrap = Bootstrap(app)

    @login.user_loader
    def load_user(user_id):
        from dlo_hic.wui.main.user import users
        return users[user_id]

    app.static_folder = 'static'
    app.password = password

    import os
    app.config.update({
        'SECRET_KEY': os.environ.get("SECRET_KEY") or gen_secret_key()
    })

    return app

