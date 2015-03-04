#!/usr/bin/env python

from setuptools import setup

setup(name='pyasb',
      version='1.0.dev0',
      description='Python - All Sky Brightness pipeline',
      author='Miguel Nievas',
      author_email='miguelnr89[at]gmail[dot]com',
      license='GPLv3',
      packages=['pyasb'],
      install_requires=['numpy', 'astropy', 'ephem', 'matplotlib'],
     )
