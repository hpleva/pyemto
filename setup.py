#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

setup(name='pyemto',
      version='1.0.0',
      description='Program to generate input files for EMTO program',
      author='H. Levämäki, M. Ropo',
      author_email='hpleva@utu.fi, matti.ropo@tut.fi',
      license='MIT',
      classifiers = ['Development Status :: 4 - Beta','Environment :: Console',
                     'Intended Audience :: Science/Research'],
      packages=['pyemto','pyemto.docs','pyemto.EOS','pyemto.emtoinputs',
                'pyemto.latticeinputs','pyemto.utilities',
                'pyemto.common','pyemto.emto_parser',
                'pyemto.examples','pyemto.free_energy',
                'pyemto.space_group'],
      package_data={"pyemto.emto_parser": ["*.json"]},
      install_requires=["numpy>=1.10.3","scipy>=0.17.1","matplotlib>=1.5.1"],
      extras_require={"emto_parser": ["pandas>=0.20.3"],
                      "emto_input_generator": ["pymatgen>=4.4.0"]},
     )
