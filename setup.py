#! /usr/bin/env python
import pdb
try:
    import setuptools  # @UnusedImport # NOQA
except ImportError:
    pass

try:
    import numpy  # @UnusedImport # NOQA
except ImportError:
    msg = ("No module named numpy. "
           "Please install numpy first, it is needed before installing ObsPy.")
    raise ImportError(msg)

#import fnmatch
import glob
import inspect
import os
import sys
import platform
from distutils.util import change_root

from setuptools import setup, find_packages
#from numpy.distutils.core import DistutilsSetupError, setup
#from numpy.distutils.ccompiler import get_default_compiler
#from numpy.distutils.command.build import build
#from numpy.distutils.command.install import install
#from numpy.distutils.exec_command import exec_command, find_executable
#from numpy.distutils.misc_util import Configuration

from rpy2.robjects.packages import importr
utils = importr('utils')
#utils.install_packages('https://cran.r-project.org/src/contrib/gemlog_0.36.tar.gz')
utils.install_packages('gemlog', repos='https://cloud.r-project.org')
# The minimum python version which can be used to run ObsPy
MIN_PYTHON_VERSION = (3, 6)

# Fail fast if the user is on an unsupported version of python.
if sys.version_info < MIN_PYTHON_VERSION:
    msg = ("gemlog requires python version >= {}".format(MIN_PYTHON_VERSION) +
           " you are using python version {}".format(sys.version_info))
    print(msg, file=sys.stderr)
    sys.exit(1)

# Directory of the current file in the (hopefully) most reliable way
# possible, according to krischer
SETUP_DIRECTORY = os.path.dirname(os.path.abspath(inspect.getfile(
    inspect.currentframe())))
INSTALL_REQUIRES = [
    'obspy',
    'rpy2',
    'numpy>=1.15.0',
    'scipy>=1.0.0',
    'matplotlib>=3.2.0',
    'lxml',
    'setuptools',
    'sqlalchemy',
    'decorator',
    'requests']

EXTRAS_REQUIRE = []
ENTRY_POINTS = {
    'console_scripts': [
        'gem2ms = gemlog.gem2ms:main'
        #'obspy-runtests = obspy.scripts.runtests:main',
    ]
}

KEYWORDS = ['']
#pdb.set_trace()
DOCSTRING = ['', '', '', '']
classifiers=[
    'Development Status :: 5 - Production/Stable',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: GNU Library or ' +
    'Lesser General Public License (LGPL)',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Physics'
]

def setupPackage():
    # setup package
    setup(
        name='gemlog',
        #version=get_git_version(),
        version = '0.0.2',
        #description=DOCSTRING[1],
        #long_description="\n".join(DOCSTRING[3:]),
        #url="",
        #author='',
        #author_email='',
        #license='',
        #platforms='OS Independent',
        #classifiers = classifiers,
        #keywords=KEYWORDS,
        #packages=['gemlog'],
        packages=find_packages(),
        #namespace_packages=[],
        #zip_safe=False,
        #python_requires=f'>={MIN_PYTHON_VERSION[0]}.{MIN_PYTHON_VERSION[1]}',
        #install_requires=INSTALL_REQUIRES,
        #extras_require=EXTRAS_REQUIRE,
        #features=add_features(),
        # this is needed for "easy_install obspy==dev"
        #download_url=("https://github.com/obspy/obspy/zipball/master"
        #              "#egg=obspy=dev"),
        #include_package_data=True,
        entry_points=ENTRY_POINTS#,
        #ext_package='obspy.lib',
        #cmdclass={
        #    'build_man': Help2ManBuild,
        #    'install_man': Help2ManInstall
        #},
        #configuration=configuration
    )


if __name__ == '__main__':
    # clean --all does not remove extensions automatically
    if 'clean' in sys.argv and '--all' in sys.argv:
        import shutil
        # delete complete build directory
        path = os.path.join(SETUP_DIRECTORY, 'build')
        try:
            shutil.rmtree(path)
        except Exception:
            pass
        # delete all shared libs from lib directory
        path = os.path.join(SETUP_DIRECTORY, 'gemlog', 'lib')
        for filename in glob.glob(path + os.sep + '*.pyd'):
            try:
                os.remove(filename)
            except Exception:
                pass
        for filename in glob.glob(path + os.sep + '*.so'):
            try:
                os.remove(filename)
            except Exception:
                pass
    else:
        setupPackage()

