PyNBody
=======

Python library for few-body numerical integrations, with a small NAIF CSPICE wrapper

## Requirements

* Python2 2.7 or greater, or Python3 3.4 or greater
* Numpy
* GCC 3.4 or greater (needs C++11 support)
* NAIF CSPICE N0065 or greater, installed such that cspice/SpiceUsr.h is in the include path and -lcspice works (see below)

### Installing CSPICE globally 

To avoid having to edit the setup.py file:

1. [Obtain CSPICE from the NAIF website](http://naif.jpl.nasa.gov/naif/toolkit_C.html)
2. Copy "cspice/include/SpiceUsr.h" to "/usr/include/cspice/SpiceUsr.h"
3. Copy "cspice/lib/cspice.a" to "/usr/lib/libcspice.a"
