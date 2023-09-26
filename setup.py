#!/usr/bin/env python
## -*- encoding: utf-8 -*-

import os
import re
import sys
from setuptools import setup
from codecs import open  # To open the README file with proper encoding
from setuptools.command.test import test as TestCommand  # for tests
from setuptools.extension import Extension
from sage.env import sage_include_directories
from Cython.Build import cythonize

# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


# For the tests
class SageTest(TestCommand):
    def run_tests(self):
        errno = os.system("sage -t --force-lib pyrforest")
        if errno != 0:
            sys.exit(1)


cythonize_dir = "build"

path = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(path, "pyrforest/lib")
allfiles_in_lib = [
    os.path.relpath(os.path.join(dp, f), path)
    for dp, dn, fn in os.walk(os.path.expanduser(lib_path))
    for f in fn
]

rforest_sources = [
    elt for elt in allfiles_in_lib if elt.endswith(".c") and elt != "test_rforest.c"
]


pyrforest = Extension(
    "pyrforest.rforest",
    language="c",
    sources=[
        "pyrforest/rforest.pyx",
    ]
    + rforest_sources,
    libraries=["gmp", "m"],
    include_dirs=sage_include_directories() + ["pyrforest/lib/"],
    extra_compile_args=[
        "-O3",
        "-fPIC",
        "-fomit-frame-pointer",
        "-funroll-loops",
        "-m64",
        "-std=gnu11",
        "-Wno-sign-compare",
        "-Wno-unused-function",
    ],
)

setup(
    name="pyrforest",
    author="Edgar Costa",
    author_email="edgarc@mit.edu",
    url="https://github.com/edgarcosta/pyrforest",
    license="MIT",
    description="Wrapper for C library rforest to compute remainder forests",
    long_description=readfile("README.md"),  # get the long description from the README
    version=readfile("VERSION"),  # the VERSION file is shared with the documentation
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT",
        "Programming Language :: Python :: 3.7",
    ],  # classifiers list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
    keywords="sagemath rforest",
    setup_requires=[
        "cython",
        "sagemath",
    ],  # currently useless, see https://www.python.org/dev/peps/pep-0518/
    install_requires=["cython", "sagemath", "sphinx"],
    packages=["pyrforest"],
    include_package_data=False,
    ext_modules=cythonize([pyrforest]),
    cmdclass={"test": SageTest}  # adding a special setup command for tests
    # ext_modules = extensions,
    # cmdclass = {'test': SageTest, 'build_ext': Cython.Build.build_ext} # adding a special setup command for tests and build_ext
)
