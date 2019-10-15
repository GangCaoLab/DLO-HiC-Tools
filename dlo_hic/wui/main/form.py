from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField, DecimalField, IntegerField
from wtforms.validators import InputRequired


class LoginForm(FlaskForm):
    password = PasswordField('Password', validators=[InputRequired()])
    submit = SubmitField('Login')


class ConfigForm(FlaskForm):
    ncpu = IntegerField('ncpu', validators=[InputRequired])
    working_dir = StringField('working-dir', validators=[InputRequired])
