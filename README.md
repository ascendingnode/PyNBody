PyNBody
=======

Python library for few-body numerical integrations

The integrator is based on [ACM Algorithm 670 by Brankin, Gladwell, Dormand, Prince, and Seward](http://dl.acm.org/citation.cfm?id=69650), converted from FORTAN 77 to C++.

## Requirements

* Python2 2.7 or greater, or Python3 3.4 or greater
* Numpy
* Cython 0.17 or greater

## Installation

To build and install locally:

```bash
python setup.py build_ext install --user
```

To build and install globally:

```bash
sudo python setup.py build_ext install
```
