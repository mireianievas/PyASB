#!/usr/bin/env python

from setuptools import setup

setup(name='pyasb',
      version='1.0.dev0',
      description='Python - All Sky Brightness pipeline',
      author='Mireia Nievas',
      author_email='mnievas[at]ucm[dot]es',
      license='GPLv3',
      packages=['pyasb'],
      install_requires=['numpy', 'scipy', 'astropy', 'ephem', 'matplotlib'],
     )
