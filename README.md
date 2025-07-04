# FreeSASA

[![DOI](https://zenodo.org/badge/18467/mittinatten/freesasa.svg)](https://zenodo.org/badge/latestdoi/18467/mittinatten/freesasa)

FreeSASA is a C library and C++ command line tool for calculating Solvent
Accessible Surface Area (SASA) of biomolecules. It is designed to be
simple to use with defaults, but allows customization of all
parameters of the calculation and provides a few different tools to
analyze the results. Python bindings are provided separately (see below).

By default Lee & Richards' algorithm is used, but Shrake & Rupley's is
also available. Both can be parameterized to arbitrary precision, and
for high resolution versions of the algorithms, the calculations give
identical results.

FreeSASA assigns a radius and a class to each atom. The atomic radii
are by default the _ProtOr_ radii defined by Tsai et
al. ([JMB 1999, 290: 253](http://www.ncbi.nlm.nih.gov/pubmed/10388571))
for standard amino acids and nucleic acids, and the van der Waals
radius of the element for other atoms. Each atom is also classified as
either polar or apolar.

Users can provide their own atomic radii and classifications via
configuration files. The input format for configuration files is
described in the
[online documentation](http://freesasa.github.io/doxygen/Config-file.html),
and the `share/` directory contains some sample configurations,
including one for the NACCESS parameters
([Hubbard & Thornton 1993](http://www.bioinf.manchester.ac.uk/naccess/)).

Version 2.0 adds some new features and breaks a few parts of the
interface from 1.x (mainly the API), see CHANGELOG.md for detailed
information.

## Building and installing

After cloning the repository, add git submodules

    git submodule init
    git submodule update

FreeSASA can be compiled and installed using the following

    autoreconf -i # only necessary if you're cloning git repo
    ./configure
    make && make install

NB: If the source is downloaded from the git repository the
configure-script needs to be set up first using `autoreconf -i`. Users
who don't have autotools installed, can download a tarball that
includes the autogenerated scripts from http://freesasa.github.io/ or
from the latest
[GitHub-release](https://github.com/mittinatten/freesasa/releases).

If you are upgrading from a pre 2.1.0 build you might need to call
`make clean` before building.

The above commands build and install the command line tool `freesasa`
(built in `src/`), the commands

    freesasa -h

and, if installed,

    man freesasa

give an overview of options. To run a calculation from PDB-file input
using the defaults, simply type

    freesasa <pdb-file>

In addition, `make install` installs the header `freesasa.h` and the
library `libfreesasa`.

The configuration can be changed with these options:

- `--disable-json` build without support for JSON output.
- `--disable-xml` build without support for XML output.
- `--disable-threads` build without multithreaded calculations
- `--enable-doxygen` activates building of Doxygen documentation

For developers:

- `--enable-check` enables unit-testing using the Check framework
- `--enable-gcov` adds compiler flags for measuring coverage of tests
  using gcov
- `--enable-parser-generator` rebuild parser/lexer source from
  Flex/Bison sources (the autogenerated code is included in the
  repository, so no need to do this if you are not going to change
  the parser).

### Installing using package managers

You can install binaries that have already been built using the package managers.

With homebrew on MacOS (amd64, arm64) or Linux (amd64):

    brew install brewsci/bio/freesasa

With conda, mamba or pixi on Linux (amd64, arm64, ppc64le) or MacOS (amd64, arm64):

    # conda
    conda install -c conda-forge freesasa-c

    # mamba
    mamba install -c conda-forge freesasa-c

    # pixi
    pixi add freesasa-c

## Python module

The Python bindings are available from PyPi and can be installed using

    pip install freesasa

This module is found in a separate repository
https://github.com/freesasa/freesasa-python
The PyPi module has binaries for Mac OS X and Windows, for a number
of Python versions and [separate documentation](http://freesasa.github.io/python/).

## Documentation

Enabling Doxygen builds a [full reference
manual](http://freesasa.github.io/doxygen/), documenting both CLI and
API in the folder `doc/html/doxygen/`, also available on the web at
http://freesasa.github.io/.

After building the package, calling

    freesasa -h

explains how the commandline tool can be used.

## Compatibility and dependencies

The library has been tested successfully with several versions of GNU
C/C++ Compiler and Clang/LLVM. It can be built using only
standard C and GNU libraries, in addition to
[Gemmi](https://github.com/project-gemmi/gemmi) (as a git submodule).
The standard build depends on
[json-c](https://github.com/json-c/json-c) and
[libxml2](http://xmlsoft.org/). These can be disabled by configuring
with `--disable-json` and `--disable-xml` respectively.

Developers who want to do testing need to install the Check unit
testing framework. Building the full reference manual requires Doxygen
(version > 1.8.8). Changing the selection parser and lexer requires Flex and
Bison. These build options, which add extra dependencies, are disabled
by default to simplify installation for users only interested in the
command line tool and and/or C Library.

The C API can be built using MSVC, see:
https://github.com/mittinatten/freesasa/issues/22#issuecomment-374661526
(no project files provided though), but probably not the command line tool.
It should be relatively straightforwad to build the command line tool for
Windows using MinGW or Cygwin, but this hasn't been tested (let me know if
you've got it to work).

### Prerequisites for Ubuntu

The following command will install all dependencies needed, some of which most users will already have,
for a minimal build of FreeSASA in Ubuntu (verified for version 16, 18 and 20).

    apt-get update
    apt-get install git build-essential autoconf libc++-dev libc++abi-dev

For a fully featured build, with ability to run unit tests, these additional packages are needed

    apt-get install check libjson-c-dev libxml2-dev libxml2-utils

### Prerequisites for RHEL / CentOS

_Thanks to [@jvlehtonen](https://github.com/jvlehtonen) for providing the following_:

GNU Build System and git install on RHEL 7 and 8 -based systems (including CentOS 7, AlmaLinux, Rocky Linux, etc) as "yum group":
   
    sudo yum install @development

The additional dependencies can be installed on RHEL-based systems with:

    sudo yum install json-c-devel libxml2-devel

The system GCC in RHEL 7 (and hence CentOS 7) does not support required C++ standard.Later version of GCC are available via Software Collections (SCL). RHEL 7 has its method to enable SCL. CentOS 7 enables SCL repo with:

    sudo yum install centos-release-scl-rh

"Toolchain" from Developer Toolset is sufficient for building with GCC. For example:

    sudo yum install devtoolset-9-toolchain

In order to use programs from SCL one has to enable environment:

    scl enable devtoolset-9 bash

The bash that starts above does have gcc version 9 on PATH. The built freesasa binary does not require the environment.

## Citing FreeSASA

FreeSASA can be cited using the following publication

- Simon Mitternacht (2016) FreeSASA: An open source C library for
  solvent accessible surface area calculations. _F1000Research_
  5:189. (doi:
  [10.12688/f1000research.7931.1](http://dx.doi.org/10.12688/f1000research.7931.1))

The [DOI numbers from Zenodo](https://zenodo.org/badge/latestdoi/18467/mittinatten/freesasa)
can be used to cite a specific version of FreeSASA.
