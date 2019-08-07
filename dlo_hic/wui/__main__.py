import sys
import click


# module level loger
import logging
log = logging.getLogger(__name__)


@click.group(name="dlohic_wui")
@click.option("--log-file", help="The log file, default output to stdout.")
@click.option("--debug", is_flag=True)
def wui(log_file, debug):
    from dlo_hic.config import LOGGING_FMT, LOGGING_DATE_FMT
    log = logging.getLogger() # set root logger
    if log_file:
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler(sys.stderr)
    fomatter = logging.Formatter(
        fmt=LOGGING_FMT,
        datefmt=LOGGING_DATE_FMT)
    handler.setFormatter(fomatter)
    log.addHandler(handler)

    if debug:
        log.setLevel(logging.DEBUG)
        import os
        os.environ['FLASK_DEBUG'] = '1'
    else:
        log.setLevel(logging.INFO)


@wui.command(name="server")
@click.option("--config-file", "-c", 
    required=True,
    help="Path to server configuration file.")
def server(config_file):
    """ Run WUI server """
    from dlo_hic.wui.app import create_app
    from dlo_hic.wui.parse_config import parse
    conf = parse(config_file)
    app = create_app(conf)
    app.run(host=conf['host'], port=int(conf['port']))


@wui.command(name="config")
def config():
    """ Generate an example server configuration file. """
    from os.path import abspath, dirname, join, exists
    here = dirname(abspath(__file__))
    fname = "server_config.json"
    example = join(here, fname)
    here_f = join("./", fname)
    if not exists(here_f):
        from subprocess import check_call
        check_call(['cp', example, here_f])
    else:
        log.info("{} aleardy exist.".format(here_f))


@wui.command(name="hash")
@click.argument("PASSWD")
def hash(passwd):
    """ Generate password hash code """
    from werkzeug import generate_password_hash
    hash_ = generate_password_hash(passwd)
    print(hash_)


wui()