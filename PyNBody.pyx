# distutils: language = c++
# distutils: extra_compile_args = -O3

####################################################

# Numpy and math
import numpy as np

# C++ Standard Library vector class
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool

####################################################

# Import the parts of the C++ NBody class we need
cdef extern from "nbody.hpp":
    cdef cppclass _NBody "NBody":
        double t , maxdist
        unsigned nobj
        vector[double] u#,maxe
        vector[string] name
        bool verbose
        _NBody() except +
        unsigned add_object_state(string,double,double,vector[double])
        unsigned add_object(string,double,double,vector[double],vector[double])
        void set_oblate(string,double,double,double,double)
        void add_emax(int,double,int)
        vector[double] get_position(int)
        vector[double] get_velocity(int)
        vector[double] get_state(int)
        double get_eccentricity(int,int)
        int lookup(string)
        void calc_bary(double &,vector[double] &,vector[double] &)
        int evolve_rkn(double,double)
        void move2bary()
        vector[double] momentum()
        double hamiltonian()
        bool is_sane()

####################################################

# Python NBody class
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

    # Define one object to be oblate
    def set_oblate(self, name,double J20,double J40,double ra0,double dec0):
        name = name.encode('UTF-8')
        self.thisptr.set_oblate(name,J20,J40,ra0,dec0)

    def add_emax(self, int i,double e,j=None):
        if j==None: self.thisptr.add_emax(i,e,-1)
        else: self.thisptr.add_emax(i,e,j)

    # Get object state
    def get_position(self, int i):
        return np.array(self.thisptr.get_position(i))
    def get_velocity(self, int i):
        return np.array(self.thisptr.get_velocity(i))
    def get_state(self, int i):
        return np.array(self.thisptr.get_state(i))

    def get_eccentricity(self, int i,j=None):
        if j==None: return self.thisptr.get_eccentricity(i,-1)
        else: return self.thisptr.get_eccentricity(i,j)

    # Lookup the index of an object given a name
    def lookup(self, string s):
        return self.thisptr.lookup(s)

    # Calculate the barycenter of an object
    def calc_bary(self):
        cdef double ub = 0
        cdef vector[double] rb,vb
        rb.resize(3); vb.resize(3)
        self.thisptr.calc_bary(ub,rb,vb)
        return ub,np.array(rb),np.array(vb)

    # Return the GM of an object
    def u(self, i):
        return self.thisptr.u[i]

    # Set the GM of an object
    def set_u(self, i,double u):
        self.thisptr.u[i] = u

    # Get the name of an index
    def name(self, i):
        return self.thisptr.name[i]

    # Set the maximum allowed eccentricity
    #def set_maxe(self, i,e):
    #    self.thisptr.maxe[i] = e

    # Get/set verbosity
    property verbose:
        def __get__(self): return self.thisptr.verbose
        def __set__(self,v): 
            if v: self.thisptr.verbose = 1
            else: self.thisptr.verbose = 0

    # Get/set time
    property t:
        def __get__(self): return self.thisptr.t
        def __set__(self,double v): self.thisptr.t = v

    # Get/set time
    property maxdist:
        def __get__(self): return self.thisptr.maxdist
        def __set__(self,double v): self.thisptr.maxdist = v

    # Get number of objects
    property nobj:
        def __get__(self): return self.thisptr.nobj

    # Make a copy
    def copy(self):
        nbs = NBody()
        nbs.thisptr[0] = self.thisptr[0]
        return nbs

    # Evolve system
    def evolve(self, double tgoal, precision=None):
        if precision==None: 
            return self.thisptr.evolve_rkn(tgoal,1e-12)
        else:
            return self.thisptr.evolve_rkn(tgoal,float(precision))

    # Move the barycenter to the center
    def move2bary(self):
        self.thisptr.move2bary()

    def momentum(self):
        return np.array(self.thisptr.momentum())

    def hamiltonian(self):
        return float(self.thisptr.hamiltonian())

    # Check that the system is sane
    def is_sane(self):
        return self.thisptr.is_sane()
