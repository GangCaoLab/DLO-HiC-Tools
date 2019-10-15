import sys
if sys.version_info < (3, 5):
    sys.exit("Only suppert >= 3.5 Python Version")

import re
import os
from codecs import open
from setuptools import setup, Extension, find_packages

from Cython.Build import cythonize

here = os.path.abspath(os.path.dirname(__file__))


extensions = [
    Extension('dlo_hic.utils.align',       sources=['dlo_hic/utils/align.pyx']),
    Extension('dlo_hic.utils.fastqio',     sources=['dlo_hic/utils/fastqio.pyx']),
    Extension('dlo_hic.utils.linker_trim', sources=['dlo_hic/utils/linker_trim.pyx']),
]


keywords = [
    "DLO Hi-C",
    "Hi-C",
    "bioinformatics",
]


def get_version():
    with open("dlo_hic/__init__.py") as f:
        for line in f.readlines():
            line = line.strip()
            m = re.match("__version__ ?= ?['\"]([^']+)['\"]", line)
            if m:
                ver = m.group(1)
                return ver
        raise IOError("Version information can not found.")


def get_long_description():
    with open("README.rst") as f:
        desc = f.read()
    return desc


def get_install_requires():
    requirements = []
    with open('requirements.txt') as f:
        for line in f:
            requirements.append(line.strip())
    return requirements


def build_wui():
    from os.path import join, exists, abspath
    wui_dir = abspath("./dlo_hic/wui")
    static_dir = join(wui_dir, 'main/static')
    from shutil import rmtree
    if exists(static_dir):  # clean old static directory
        rmtree(static_dir)
    curr_dir = abspath(os.curdir)
    os.chdir(wui_dir)
    from subprocess import check_call
    npm_dir = join(wui_dir, 'node_modules')
    print(npm_dir, exists(npm_dir))
    if not exists(npm_dir):
        check_call(['npm', 'install'])
    check_call(['npm', 'run', 'build'])
    os.chdir(curr_dir)


def check_npm_installed():
    from subprocess import Popen, PIPE
    try:
        p = Popen(['npm', '--version'], stdout=PIPE)
        out, _ = p.communicate()
        version = out.decode('utf-8').strip()
        print("npm version: "+version)
        p.kill()
        if p.returncode == 0:
            return True
        else:
            return True
    except FileNotFoundError as e:
        print("[Warning] npm is not installed. Don't build WUI.")
        return False


if check_npm_installed():
    try:
        build_wui()
    except Exception as e:
        print(str(e))
        print("[Error] Fail to build WUI.")


setup(
    name='dlo_hic',
    version=get_version(),
    description='Tools for DLO-HiC data analyze',
    long_description=get_long_description(),
    author='nanguage',
    author_email='nanguage@yahoo.com',
    url='https://github.com/Nanguage/DLO-HiC-Tools',
    keywords=keywords,
    packages=find_packages(),
    scripts=['scripts/dlohic', 'scripts/dlohicwui'],
    include_package_data=True,
    zip_safe=False,
    package_dir={'dlo_hic': 'dlo_hic'},
    install_requires=get_install_requires(),
    ext_modules=cythonize(extensions),
    license='GNU GPLv3',
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        "Programming Language :: Cython",
    ]
)
