#!/usr/bin/env python2

from distutils.core import setup

setup(
    name='parametric_protein_scaffold_design',
    version='0.0.0',
    author='Xingjie Pan',
    author_email='xingjiepan@gmail.com',
    url='https://github.com/Kortemme-Lab/parametric_scaffold_design',
    packages=[
        'parametric_protein_scaffold_design',
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'docopt',
    ],
    extras_require = {
        'weblogo':  ['weblogo'],
        'cylinder_fitting': ['cylinder_fitting']
    },
    entry_points={
        'console_scripts': [
        ],
    },
    description='PyRosetta scripts for parametric protein scaffold design',
    long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 2',
        'Intended Audience :: Science/Research',
    ],
)
