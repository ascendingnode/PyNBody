# distutils: language = c++
# distutils: extra_compile_args = -std=c++11 -O3

# Numpy and math
import numpy as np
#cimport numpy as np
import math
import scipy.optimize as opt

# C++ Standard Library vector class
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

####################################################

# Import the parts of the C++ Conic class we need
cdef extern from "conic.hpp":
    cdef cppclass _Conic "Conic":
        double rp,e,i,O,w,M0,t0,mu
        _Conic() except +
        void setup_elements(vector[double])
        void setup_state(double,double,vector[double])
        vector[double] elements(double)
        vector[double] state(double)

# Python Conic class
cdef class Conic:

    # Contain a copy of the C++ Conic class
    cdef _Conic *thisptr
    def __cinit__(self):
        self.thisptr = new _Conic()
    def __del__(self):
        del self.thisptr

    # Directly get the orbital elements
    property rp:
        def __get__(self): return self.thisptr.rp
    property e:
        def __get__(self): return self.thisptr.e
    property i:
        def __get__(self): return self.thisptr.i
    property O:
        def __get__(self): return self.thisptr.O
    property w:
        def __get__(self): return self.thisptr.w
    property M0:
        def __get__(self): return self.thisptr.M0
    property t0:
        def __get__(self): return self.thisptr.t0
    property mu:
        def __get__(self): return self.thisptr.mu

    # Setup the Conic class
    def setup_elements(self, elts):
        self.thisptr.setup_elements(elts)
    def setup_state(self, t0,mu,sv):
        self.thisptr.setup_state(t0,mu,sv)
    def setup_rv(self, t0,mu,r,v):
        sv = list(r)+list(v)
        self.thisptr.setup_state(t0,mu,sv)

    # Extract elements or state from the class
    def elements(self, t):
        return np.array(self.thisptr.elements(t))
    def state(self, t):
        return np.array(self.thisptr.state(t))
    def rv(self, t):
        sv = self.state(t)
        return sv[:3],sv[3:]

####################################################

# Import the parts of the C++ NBody class we need
cdef extern from "nbody.hpp":
    cdef cppclass _NBody "NBody":
        double T, maxdist
        unsigned nobj
        int failed
        vector[double] u,maxe
        vector[string] name
        bool verbose
        _NBody() except +
        unsigned add_object_state(string,double,double,vector[double])
        unsigned add_object(string,double,double,vector[double],vector[double])
        vector[double] get_position(int)
        vector[double] get_velocity(int)
        vector[double] get_state(int)
        int lookup(string)
        void calc_bary(double &,vector[double] &,vector[double] &)
        _Conic orbit(int,int)
        void evolve_self(double)
        void move2bary()
    void rkn_evolve(_NBody &,double)

# Python nbody class
cdef class NBody:
    
    # Contain a copy of the C++ Conic class
    cdef _NBody *thisptr
    def __cinit__(self):
        self.thisptr = new _NBody()
    def __del__(self):
        del self.thisptr

    # Add an object to the system
    def add_object(self, name,double u0,double R0,vector[double] r0,vector[double] v0):
        name = name.encode('UTF-8')
        return self.thisptr.add_object(name, u0,R0,r0,v0)
    def add_object_state(self, name,double u0,double R0,vector[double] state):
        name = name.encode('UTF-8')
        return self.thisptr.add_object_state(name, u0,R0,state)

    # Get object state
    def get_position(self, int i):
        return np.array(self.thisptr.get_position(i))
    def get_velocity(self, int i):
        return np.array(self.thisptr.get_velocity(i))
    def get_state(self, int i):
        return np.array(self.thisptr.get_state(i))

    # Get a relative orbit; j<0 means barycentric
    def orbit(self, int i,int j=-1):
        o = Conic()
        o.thisptr[0] = self.thisptr.orbit(i,j)
        return o

    def lookup(self, string s):
        return self.thisptr.lookup(s)

    def calc_bary(self):
        cdef double ub = 0
        cdef vector[double] rb,vb
        rb.resize(3); vb.resize(3)
        self.thisptr.calc_bary(ub,rb,vb)
        return ub,np.array(rb),np.array(vb)

    def u(self, i):
        return self.thisptr.u[i]

    def set_u(self, i,double u):
        self.thisptr.u[i] = u

    def name(self, i):
        return self.thisptr.name[i]

    def set_maxe(self, i,e):
        self.thisptr.maxe[i] = e

    # Get/set verbosity
    property verbose:
        def __get__(self): return self.thisptr.verbose
        def __set__(self,v): 
            if v: self.thisptr.verbose = 1
            else: self.thisptr.verbose = 0

    # Get/set time
    property T:
        def __get__(self): return self.thisptr.T
        def __set__(self,double v): self.thisptr.T = v

    # Get/set time
    property maxdist:
        def __get__(self): return self.thisptr.maxdist
        def __set__(self,double v): self.thisptr.maxdist = v

    # Get number of objects
    property nobj:
        def __get__(self): return self.thisptr.nobj

    # Get failure state
    property failed:
        def __get__(self): return self.thisptr.failed

    # Make a copy
    def copy(self):
        nbs = NBody()
        nbs.thisptr[0] = self.thisptr[0]
        return nbs

    # Evolve system
    def evolve(self, double tgoal):
        self.thisptr.evolve_self(tgoal)
        #cdef _NBody nb = self.thisptr[0]
        #rkn_evolve(nb,tgoal)
        #self.thisptr[0] = nb

    # Move the barycenter to the center
    def move2bary(self):
        self.thisptr.move2bary()
