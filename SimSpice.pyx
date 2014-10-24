# distutils: libraries = cspice

# Import numpy 
cimport numpy as np
import numpy as np

# Import NAIF SPICE functions
cdef extern from "cspice/SpiceUsr.h":
    void furnsh_c(char *)
    void unload_c(char *)
    void kclear_c()
    void spkgps_c(int,double,char *,int,double *,double *)
    void spkgeo_c(int,double,char *,int,double *,double *)
    void spkpos_c(char *,double,char *,char *,char *,double *,double *)
    void str2et_c(char *,double *)
    void timout_c(double,char *,int,char *)
    void tpictr_c(char *,int,int,char *,int *,char *)
    void recrad_c(double *,double *,double *,double *)
    void radrec_c(double,double,double,double *)

def furnsh(kernel):
    k2 = kernel.encode('UTF-8')
    cdef char* k3 = k2
    furnsh_c(k3)

def unload(kernel):
    k2 = kernel.encode('UTF-8')
    cdef char* k3 = k2
    unload_c(k3)

def kclear():
    kclear_c()

def spkgps(int targ,double et,ref,int obs):
    ref = ref.encode('UTF-8')
    cdef char* ref2 = ref
    cdef double lt        
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(3)
    spkgps_c(targ,et,ref2,obs,&state[0],&lt)
    return state,lt

def spkgeo(int targ,double et,ref,int obs):
    ref = ref.encode('UTF-8')
    cdef char* ref2 = ref
    cdef double lt
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(6)
    spkgeo_c(targ,et,ref2,obs,&state[0],&lt)
    return state,lt

def spkpos(targ,double et,ref,abcorr,obs):
    targ = targ.encode('UTF-8')
    ref = ref.encode('UTF-8')
    abcorr = abcorr.encode('UTF-8')
    obs = obs.encode('UTF-8')
    cdef char* targ2 = targ
    cdef char* ref2 = ref
    cdef char* abcorr2 = abcorr
    cdef char* obs2 = obs
    cdef double lt
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(3)
    spkpos_c(targ2,et,ref2,abcorr2,obs2,&state[0],&lt)
    return state,lt

def str2et(s):
    s2 = s.encode('UTF-8')
    cdef char* s3 = s2
    cdef double et
    str2et_c(s3,&et)
    return et

def timout(double et,picture):
    picture = picture.encode('UTF-8')
    cdef char* pic2 = picture
    cdef char buff[256]
    timout_c(et,pic2,256,buff)
    return buff

def tpictr(example):
    example = example.encode('UTF-8')
    cdef char* example2 = example
    cdef char buff1[256]
    cdef char buff2[256]
    cdef int ok = 0
    tpictr_c(example2,256,256,buff1,&ok,buff2)
    return buff1,ok,buff2

def recrad(rectan):
    cdef double ran,ra,dec
    cdef np.ndarray[double, ndim=1, mode="c"] state = rectan
    recrad_c(&state[0],&ran,&ra,&dec)
    return ran,ra,dec

def radrec(double ran,double ra,double dec):
    cdef np.ndarray[double, ndim=1, mode="c"] state = np.zeros(3)
    radrec_c(ran,ra,dec,&state[0])
    return state
