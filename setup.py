#!/usr/bin/env python
 # -*- coding: utf-8 -*-

from distutils.core import setup

setup(name='pyEmto',
      version='1.0',
      description='Program to generate input files for EMTO prgoram',
      author='H. Levämäki, M. Ropo',
      author_email='hpleva@utu.fi, matti.ropo@tut.fi',
      licence =  'License :: Free for non-commercial use',
      classifiers = ['Development Status :: 4 - Beta','Environment :: Console',
                     'Intended Audience :: Science/Research'],
#      url='https://www.python.org/sigs/distutils-sig/',
      packages=['pyemto','pyemto.EOS','pyemto.emtoinputs',
                'pyemto.latticeinputs','pyemto.utilities',
                'pyemto.common','pyemto.emto_parser',
                'pyemto.examples'],
     )
