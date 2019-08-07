from uuid import uuid1
from werkzeug.security import check_password_hash


users = {}


class User():

    def __init__(self, password_hash):
        self.uid = str(uuid1())
        self.password_hash = password_hash
        users[self.uid] = self

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def get_id(self):
        return self.uid

    @property
    def is_authenticated(self):
        return True

    @property
    def is_active(self):
        return True

    @property
    def is_anonymous(self):
        return False