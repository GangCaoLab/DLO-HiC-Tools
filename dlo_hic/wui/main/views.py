from flask import (render_template, Blueprint, url_for,
                   redirect, flash, request,
                   current_app)
from flask_login import current_user, login_required, login_user, logout_user


main = Blueprint('main', __name__,
                 template_folder='templates',
                 static_folder='static',
                 static_url_path='/main/static')


import logging
log = logging.getLogger(__name__)


@main.route("/", methods=['GET', 'POST'])
def index():
    if not current_user.is_authenticated:
        return redirect(url_for('main.login'))
    return render_template("index.html", conf=current_app.conf)


@login_required
@main.route("/gen-pipe-conf", methods=['GET', 'POST'])
def gen_pipe_conf():
    from .form import ConfigForm
    form = ConfigForm()
    if request.method == 'POST':
        log.debug({k:v.data for k, v in form.__dict__.items if 'data' in v.__dict__})


@main.route("/login", methods=['GET', 'POST'])
def login():
    if current_user.is_authenticated:
        return redirect(url_for('main.index'))
    from .form import LoginForm
    form = LoginForm()
    if request.method == 'POST' and form.validate():
        from .user import User
        user = User(current_app.conf['password_hash'])
        if not user.check_password(form.password.data):
            flash("Invalid password.")
            return redirect(url_for('main.login'))
        login_user(user, remember=True)
        return redirect(url_for('main.index'))

    return render_template('login.html', title='Sign In', form=form)


@main.route("/logout")
def logout():
    logout_user()
    return redirect(url_for('main.index'))
