#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os
from io import open

about = {}
# Read version number from deepbgc.__version__.py (see PEP 396)
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'biophi', '__version__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

# Read contents of readme file into string
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Read dependencies from environment.yml
with open(os.path.join(here, 'environment.yml'), encoding='utf-8') as f:
    # We need to parse without yaml parser, since this runs before any packages are installed
    install_requires = []
    started = False
    for line in f:
        if started:
            install_requires.append(line.strip().replace('- ', ''))
        if 'pip:' in line:
            started = True
    if not install_requires:
        raise ValueError(f'Error parsing pip dependencies from environment.yml')

setup(
    name='biophi',
    version=about['__version__'],
    description='BioPhi: Platform for antibody design, humanization and humanness evaluation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='David Prihoda',
    packages=find_packages(include=['biophi.*', 'biophi']),
    author_email='david.prihoda@gmail.com',
    license='MIT',
    python_requires=">=3.8",
    install_requires=install_requires,
    keywords='biophi, antibody humanization, humanness evaluation, sapiens, oasis',
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
    include_package_data=True,
    package_data={'': ['*.js', '*.css', '*.html', '*.png', '*.svg', '*.json']},
    url='https://github.com/Merck/BioPhi',
    entry_points={
        'console_scripts': ['biophi = biophi.common.cli.main:main']
    }
)
