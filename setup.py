#!/usr/bin/env python

# import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize


sources = ["cjam/_jam_axi.pyx"]
interp = ["src/interp/interp2dpol.c"]
jam = ["src/jam/jam_axi_rms.c", "src/jam/jam_axi_rms_mgeint.c",
    "src/jam/jam_axi_rms_mmt.c", "src/jam/jam_axi_rms_wmmt.c",
    "src/jam/jam_axi_vel.c", "src/jam/jam_axi_vel_losint.c",
    "src/jam/jam_axi_vel_mgeint.c", "src/jam/jam_axi_vel_mmt.c",
    "src/jam/jam_axi_vel_wmmt.c"]
mge = ["src/mge/mge_addbh.c", "src/mge/mge_dens.c", "src/mge/mge_deproject.c",
    "src/mge/mge_qmed.c", "src/mge/mge_read.c", "src/mge/mge_surf.c"]
tools = ["src/tools/maximum.c", "src/tools/median.c", "src/tools/minimum.c",
    "src/tools/range.c", "src/tools/readcol.c", "src/tools/sort_dbl.c",
    "src/tools/where.c"]
sources += interp + jam + mge + tools

ext_modules = Extension("cjam._jam_axi", sources, libraries=["gsl","gslcblas"])

setup(name="cjam",
    ext_modules=cythonize([ext_modules]),
    version="0.1",
    description="Jeans Anisotropic MGE models -- written in C with a Python "
        "wrapper.",
    author="Laura Watkins",
    author_email="lauralwatkins@gmail.com",
    url="https://github.com/lauralwatkins/cjam",
    packages=["cjam"],
    )
