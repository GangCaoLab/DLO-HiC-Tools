import logging
log = logging.getLogger(__name__)

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


def create_app(conf):
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

    app.conf = conf
    from os.path import getmtime
    app.conf_modify_time = getmtime(conf['path'])

    from dlo_hic.wui.parse_config import parse
    @app.before_request
    def reload_config():  # update config when server_config is modified
        if app.conf_modify_time != getmtime(conf['path']):
            log.info("Detected change in {}, reloading".format(conf['path']))
            app.conf = parse(conf['path'])

    import os
    app.config.update({
        'SECRET_KEY': os.environ.get("SECRET_KEY") or gen_secret_key()
    })

    return app

