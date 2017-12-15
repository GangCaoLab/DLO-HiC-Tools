import os
from codecs import open
from setuptools import setup, Extension, find_packages

from Cython.Build import cythonize

here = os.path.abspath(os.path.dirname(__file__))

# load about information
requires = [
    "numpy >= 1.13.3",
    "scipy >= 0.18.1",
    "matplotlib >= 2.0.0",
    "biopython >= 1.70",
    "iced >= 0.4.2",
    "Cython >= 0.25.2",
]

extensions = [
    Extension('dlo_hic.utils.align', sources=['dlo_hic/utils/align.pyx'])
]

setup(
    name='dlo_hic',
    version='0.0.0',
    keywords='HiC',
    description='Tools for DLO-HiC data analyze',
    author='nanguage',
    author_email='nanguage@yahoo.com',
    url='',
    packages=find_packages(),
    package_data={
        '':['LICENSE'],
        'dlo_hic': [],
        },
    package_dir={'dlo_hic': 'dlo_hic'},
    install_requires=requires,
    ext_modules=cythonize(extensions),
    license='MIT License',
    classifiers=[
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        "Programming Language :: Cython",
    ]
)
