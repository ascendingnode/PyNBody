PyNBody
=======

Python library for few-body numerical integrations, with a small NAIF CSPICE wrapper

## Requirements

* Python2 2.7 or greater, or Python3 3.4 or greater
* Numpy
* GCC 3.4 or greater (needs C++11 support)
* NAIF CSPICE N0065 or greater; see below

### Download CSPICE 

To avoid having to edit the setup.py file:

1. Download the appropriate version of CSPICE [from the NAIF website](http://naif.jpl.nasa.gov/naif/toolkit_C.html) into the PyNBody directory
2. Extract the tarball: `tar xzvf cspice.tar.Z`

## Installation

To build and install locally:

python setup.py build_ext install --user

To build and install globally:

sudo python setup.py build_ext install
