import sys
if sys.version_info < (3, 4):
    sys.exit("Only suppert >= 3.4 Python Version")

import re
import os
from codecs import open
from setuptools import setup, Extension, find_packages

from Cython.Build import cythonize

here = os.path.abspath(os.path.dirname(__file__))


extensions = [
    Extension('dlo_hic.utils.align', sources=['dlo_hic/utils/align.pyx'])
]


def get_version():
    with open("dlo_hic/__init__.py") as f:
        for line in f.readlines():
            m = re.match("__version__ = '([^']+)'", line)
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


setup(
    name='dlo_hic',
    version=get_version(),
    keywords='HiC',
    description='Tools for DLO-HiC data analyze',
    long_description=get_long_description(),
    author='nanguage',
    author_email='nanguage@yahoo.com',
    url='https://github.com/Nanguage/DLO-HiC-Tools',
    packages=find_packages(),
    scripts=['scripts/dlohic'],
    include_package_data=True,
    package_dir={'dlo_hic': 'dlo_hic'},
    install_requires=get_install_requires(),
    ext_modules=cythonize(extensions),
    license='GNU GPLv3',
    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        "Programming Language :: Cython",
    ]
)
