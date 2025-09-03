import os
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.install import install as _install

PKG = "graphy"

def make_louvain():
    print("\n*** Running make on Louvain code...")
    subprocess.call("cd graphy/external/SurpriseMeCPM ; make", shell=True)
    print()

class install(_install):
    def run(self):
        make_louvain()
        super().run()

class build_ext(_build_ext):
    def run(self):
        make_louvain()
        super().run()
    def finalize_options(self):
        super().finalize_options()
        import numpy as np
        self.include_dirs.append(np.get_include())



# Cython modules (if any)
cython_modules = []
try:
    from Cython.Build import cythonize
    ext_modules = cythonize([os.path.join(PKG, s + ".pyx") for s in cython_modules])
except ImportError:
    ext_modules = [Extension(PKG + "." + s, [os.path.join(PKG, s + ".c")]) for s in cython_modules]

setup(
    name=PKG,
    cmdclass={"build_ext": build_ext, "install": install},
    ext_modules=ext_modules,
)
