#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("pyemto.c.c_lattice", ["pyemto/c/c_lattice.pyx"],
              include_dirs = [numpy.get_include()],
              #libraries = [...],
              #library_dirs = [...],
              )
]

setup(name='pyemto',
      version='0.9.3',
      description='Program to generate input files for EMTO program',
      author='H. Levämäki, M. Ropo',
      author_email='hpleva@utu.fi, matti.ropo@tut.fi',
      license='MIT',
      classifiers=['Development Status :: 4 - Beta', 'Environment :: Console',
                   'Intended Audience :: Science/Research'],
      packages=['pyemto', 'pyemto.EOS', 'pyemto.emtoinputs',
                'pyemto.latticeinputs', 'pyemto.utilities',
                'pyemto.common', 'pyemto.emto_parser',
                'pyemto.examples', 'pyemto.free_energy',
                'pyemto.space_group', 'pyemto.jij',
                'pyemto.c'],
      package_data={"pyemto.emto_parser": ["*.json", "LICENSE.txt"],
                    "pyemto": ["contributors.txt", "Documentation.html"],
                    "pyemto.c": ["*.py", "*.so", "*.pyc"]},
      install_requires=["numpy>=1.10.3", "scipy>=0.17.1", "matplotlib>=1.5.1"],
      extras_require={"emto_parser": ["pandas>=0.20.3"],
                      "emto_input_generator": ["pymatgen>=4.4.0"]},
      ext_modules = cythonize(extensions),
      scripts=['pyemto/jij/create_jij_input']
      )
