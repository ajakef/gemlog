#!/usr/bin/env python
try:
    from setuptools import setup, Extension
except ImportError:
    raise RuntimeError('setuptools is required')

import os

from distutils.util import convert_path
import sys

# The minimum python version which can be used to run ObsPy
MIN_PYTHON_VERSION = (3, 8)

# Fail fast if the user is on an unsupported version of python.
if sys.version_info < MIN_PYTHON_VERSION:
    msg = ("gemlog requires python version >= {}".format(MIN_PYTHON_VERSION) +
           " you are using python version {}".format(sys.version_info))
    print(msg, file=sys.stderr)
    sys.exit(1)


DESCRIPTION = 'A set of functions for processing Gem datalogger files.'
LONG_DESCRIPTION = """
gemlog is a python package that provides functions for processing the
log files from the Gem infrasound datalogger. The Gem infrasound logger
is a low-cost, lightweight, low-power instrument for recording infrasound
in field campaigns.

Source code: https://github.com/ajakef/gemlog
"""

DISTNAME = 'gemlog'
LICENSE = 'GPL-3.0'
AUTHOR = 'Jake Anderson'
MAINTAINER_EMAIL = 'ajakef@gmail.com'
URL = 'https://github.com/ajakef/gemlog'


# TO-DO: replace this with versioneer
version_dict = {}
version_path = convert_path('gemlog/version.py')
with open(version_path) as version_file:
    exec(version_file.read(), version_dict)

VERSION = version_dict['__version__']

## Dependency notes:
## Main version requirements
# pandas: >= 1.3 (2021-07; earlier versions are incompatible with all gemlog. >=1.4 only works for gemlog >= 1.6.1)
# numpy: >=1.22 (2022-06-22; earlier versions have security issue https://github.com/advisories/GHSA-fpfv-jqm9-f5jm).
# scipy: >=1.4 (2022-11-22; 1.3.0 fails to install now, don't know why)

## Secondary requirements
# python: >=3.8 (2022-06-22; >=3.8 is required by numpy 1.22).
# obspy: >=1.3 (2022-06-22; >=1.3 is required for numpy 1.22 compatibility)

## example range: pandas>=1.3.0,<1.4.0
## example exact: numpy==1.21
INSTALL_REQUIRES = [
    'setuptools>=18.0', # 18.0 needed to handle cython in installation
    'obspy>=1.3',
    'numpy>=1.22', 
    'pandas>=1.3.0', 
    'scipy>=1.4.0', # May 2019
    'matplotlib>=3.2.0', # March 2020
    'cython', # used for reading raw files quickly
    'fpdf' # needed for pdf huddle test reports
]

TESTS_REQUIRE = ['pytest']

EXTRAS_REQUIRE = {
    # 'optional': [...],
    'doc': ['sphinx==5.3.0', 'sphinx-rtd-theme==1.1.1'],
    'test': TESTS_REQUIRE
}
EXTRAS_REQUIRE['all'] = sorted(set(sum(EXTRAS_REQUIRE.values(), [])))

CLASSIFIERS = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Physics'
]

ENTRY_POINTS = {
    'console_scripts': [
        'gem2ms = gemlog.gemconvert:main',
        'gemconvert = gemlog.gemconvert:main',
        'gemconvert_single = gemlog.gemconvert_single:main',
        'gem_cat = gemlog.gem_cat:main',
        'change_mseed_times = gemlog.change_mseed_times:main',
        'verify_huddle_test = gemlog.huddle_test:main',
        'waveform_calc_lags = gemlog.xcorr:xcorr_all_terminal',
        'waveform_calc_directions = gemlog.xcorr:calculate_direction_terminal',
        'gem_make_inventory = gemlog.gem_network:summarize_gps_terminal'
    ]
}

setuptools_kwargs = {
    'zip_safe': False,
    'scripts': [],
    'include_package_data': True,
}

PACKAGES = ['gemlog']

setup(name=DISTNAME,
      setup_requires=['setuptools>=18.0','cython'],
      version=VERSION,
      packages=PACKAGES,
      install_requires=INSTALL_REQUIRES,
      extras_require=EXTRAS_REQUIRE,
      tests_require=TESTS_REQUIRE,
      ext_modules=[Extension('gemlog.parsers', sources=['gemlog/parsers.pyx'])],
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      maintainer_email=MAINTAINER_EMAIL,
      license=LICENSE,
      url=URL,
      classifiers=CLASSIFIERS,
      entry_points=ENTRY_POINTS,
      **setuptools_kwargs)
