import click

import os

os.environ['FLASK_DEBUG'] = '1'

from dlo_hic.wui.app import create_app

app = create_app("test")
app.run()
