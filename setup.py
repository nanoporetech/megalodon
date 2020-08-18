import os
import re
import sys

from setuptools import setup
from setuptools.extension import Extension


__pkg_name__ = 'megalodon'

# Get the version number from _version.py, and exe_path
verstrline = open(
    os.path.join(os.path.dirname(__file__), 'megalodon', '_version.py'),
    'r').readlines()[-1]
vsre = r"^MEGALODON_VERSION = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError(
        'Unable to find version string in "megalodon/_version.py".')


install_requires = [
    "h5py >= 2.2.1",
    "numpy >= 1.9.0",
    "scipy >= 1.1.0",
    "Cython >= 0.25.2",
    "mappy >= 2.16",
    "pysam >= 0.15",
    "ont_fast5_api >= 1.1",
    "tqdm",
    "ont-pyguppy-client-lib",
    "sklearn",
    "seaborn"
]


#  Build extensions
try:
    import numpy as np
    include_dirs = [np.get_include()]
except ImportError:
    sys.stderr.write(
        '*' * 60 + '\nINSTALLATION ERROR:\n'
        '\tNeed to install numpy before megalodon installation.\n' +
        '\tThis is required in order to get maximum efficincy from ' +
        'cython code optimizations.\n' +
        'To install run:\n$ pip install numpy\n' +
        '*' * 60 + '\n')
    sys.exit(1)

extra_compile_args = ['-std=c99']
if sys.platform == 'darwin':
    extra_compile_args.append('-mmacosx-version-min=10.9')
    print('Using macOS clang args')
extensions = [
    Extension(str("megalodon.decode"),
              [str("megalodon/_decode.pyx")],
              include_dirs=include_dirs,
              extra_compile_args=extra_compile_args,
              language="c"),
]
extensions[0].cython_directives = {"embedsignature": True}


def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name=__pkg_name__,
    version=__version__,
    description='Nanopore base calling augmentation',
    long_description=readme(),
    keywords='nanopore sequencing basecalling mapping methylation variants',
    author='Marcus Stoiber',
    maintainer='Marcus Stoiber',
    maintainer_email='marcus.stoiber@nanoporetech.com',
    url='https://github.com/nanoporetech/megalodon',
    license='mpl-2.0',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: GPU',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    packages=['megalodon', 'megalodon_extras'],
    package_data={__pkg_name__: ['model_data/*/*', ]},
    exclude_package_data={'': ['*.hdf', '*.c', '*.h']},
    ext_modules=extensions,
    install_requires=install_requires,
    zip_safe=False,
    entry_points={
        'console_scripts': [
            '{0} = {0}.__main__:_main'.format(__pkg_name__),
            '{0}_extras = {0}_extras.__main__:_main'.format(__pkg_name__),
        ]
    },
)
