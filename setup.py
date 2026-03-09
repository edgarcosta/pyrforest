import os
from setuptools import setup
from setuptools.extension import Extension

try:
    from sage.config import get_include_dirs
    sage_include_directories = lambda: [str(p) for p in get_include_dirs()]
except ImportError:
    from sage.env import sage_include_directories

from Cython.Build import cythonize as cython_cythonize

try:
    from sage.misc.package_dir import cython_namespace_package_support
    def cythonize(*args, **kwargs):
        with cython_namespace_package_support():
            return cython_cythonize(*args, **kwargs)
except ImportError:
    cythonize = cython_cythonize

path = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(path, "pyrforest/lib")
rforest_sources = [
    os.path.relpath(os.path.join(dp, f), path)
    for dp, _, fn in os.walk(lib_path)
    for f in fn
    if f.endswith(".c") and f != "test_rforest.c"
]

rforest = Extension(
    "pyrforest.rforest",
    language="c",
    sources=["pyrforest/rforest.pyx"] + rforest_sources,
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
    ext_modules=cythonize([rforest], include_path=sage_include_directories()),
)
