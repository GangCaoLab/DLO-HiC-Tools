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
    scripts=['scripts/dlohic'],
    include_package_data=True,
    zip_safe=False,
    package_dir={'dlo_hic': 'dlo_hic'},
    install_requires=[],
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
