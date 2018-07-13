import sys
if sys.version_info < (3, 4):
    sys.exit("Only suppert >= 3.4 Python Version")

import os
from codecs import open
from setuptools import setup, Extension, find_packages

from Cython.Build import cythonize

here = os.path.abspath(os.path.dirname(__file__))

# load about information
requires = [
    "numpy >= 1.13.0",
    "scipy >= 0.18.0",
    "matplotlib >= 2.0.0",
    "biopython >= 1.70",
    "Cython >= 0.25.2",
    "click >= 6.7",
    "seaborn >= 0.8.1",
    "pandas >= 0.20.3",
    "pyfaidx >= 0.5.1"
    "cooler >= 0.7.6"
]

extensions = [
    Extension('dlo_hic.utils.align', sources=['dlo_hic/utils/align.pyx'])
]

setup(
    name='dlo_hic',
    version='0.0.1',
    keywords='HiC',
    description='Tools for DLO-HiC data analyze',
    author='nanguage',
    author_email='nanguage@yahoo.com',
    url='https://github.com/Nanguage/DLO-HiC-Tools',
    packages=find_packages(),
    scripts=['scripts/dlohic'],
    package_data={
        '':['LICENSE'],
        'dlo_hic': [],
    },
    include_package_data=True,
    package_dir={'dlo_hic': 'dlo_hic'},
    install_requires=requires,
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
