from flask import (render_template, Blueprint, url_for,
                   redirect, flash, request,
                   current_app)
from flask_login import current_user, login_required, login_user

from .form import LoginForm
from .user import User

main = Blueprint('main', __name__,
                 template_folder='templates',
                 static_folder='static',
                 static_url_path='/main/static')


@main.route("/")
def index():
    if not current_user.is_authenticated:
        return redirect(url_for('main.login'))
    return render_template("index.html")


@main.route("/login", methods=['GET', 'POST'])
def login():
    if current_user.is_authenticated:
        return redirect(url_for('main.index'))
    form = LoginForm()
    if request.method == 'POST' and form.validate():
        user = User(current_app.password)
        if not user.check_password(form.password.data):
            flash("Invalid password.")
            return redirect(url_for('main.login'))
        login_user(user, remember=True)
        return redirect(url_for('main.index'))

    return render_template('login.html', title='Sign In', form=form)
