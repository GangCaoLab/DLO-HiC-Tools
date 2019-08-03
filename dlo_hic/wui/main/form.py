from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField
from wtforms.validators import InputRequired


class LoginForm(FlaskForm):
    password = PasswordField('Password', validators=[InputRequired()])
    submit = SubmitField('Login')
