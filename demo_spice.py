from __future__ import print_function
from PyNBody import NBody
import SimSpice as spice
import numpy as np

# Create the Uranus system from SPICE
def make_uranus(et):
    
    # GMs for massive objects
    gm_ura = 5.793951322279009E+06 # 799
    gm_ari = 8.346344431770477E+01 # 701
    gm_umb = 8.509338094489388E+01 # 702
    gm_tit = 2.269437003741248E+02 # 703
    gm_obe = 2.053234302535623E+02 # 704
    gm_mir = 4.319516899232100E+00 # 705

    # Radii for objects
    R_ura = 2.555900000000000E+04 # 799

    # Basic obliquity parameters
    J702   =  3.510685384697763E-03 # Uranus J2
    J704   = -3.416639735448987E-05 # Uranus J4
    ZACPL7 =  7.730990252631723E+01 # RA of Uranus pole at epoch
    ZDEPL7 =  1.517245819840212E+01 # Dec of Uranus pole at epoch
    
    nb = NBody()
    nb.t = et

    nb.add_object_state("Uranus",gm_ura,R_ura,np.zeros(6))
    nb.set_oblate("Uranus",J702,J704,ZACPL7,ZDEPL7)

    state = spice.spkgeo(705,et,"j2000",799)[0]
    nb.add_object_state("Miranda",gm_mir,0,state)

    state = spice.spkgeo(701,et,"j2000",799)[0]
    nb.add_object_state("Ariel",gm_ari,0,state)

    state = spice.spkgeo(702,et,"j2000",799)[0]
    nb.add_object_state("Umbriel",gm_umb,0,state)

    state = spice.spkgeo(703,et,"j2000",799)[0]
    nb.add_object_state("Titania",gm_tit,0,state)

    state = spice.spkgeo(704,et,"j2000",799)[0]
    nb.add_object_state("Oberon",gm_obe,0,state)

    nb.move2bary()

    return nb

# Measure RMS error of an integrated NBody system
def print_diff(et0,nb1):
    et = nb1.t
    nb0 = make_uranus(et)
    dt = (et-et0)/(24.*3600.)
    print("Runtime: {:.3e} years, {:.3e} days\n".format(dt/365.25,dt))
    drt,dvt = (0,0)
    for i in range(1,nb1.nobj):
        r0 = nb0.get_position(i)-nb0.get_position(0)
        v0 = nb0.get_velocity(i)-nb0.get_velocity(0)
        r1 = nb1.get_position(i)-nb1.get_position(0)
        v1 = nb1.get_velocity(i)-nb1.get_velocity(0)
        dr = np.sqrt(np.mean(np.square(r0-r1)))
        dv = np.sqrt(np.mean(np.square(v0-v1)))
        drt += np.sum(np.square(r0-r1))
        dvt += np.sum(np.square(v0-v1))
        print("{:<10} {:.8e} {:.8e}".format(nb0.name(i),dr,dv))
    drt = np.sqrt(drt/(3.0*nb1.nobj))
    dvt = np.sqrt(dvt/(3.0*nb1.nobj))
    print("\n{:<10} {:.8e} {:.8e}".format("Total",drt,dvt))

# Integrate the Uranus system
def run_uranus(et1,et2):
    nb0 = make_uranus(et1)
    nb1 = nb0.copy()
    ifail = nb1.evolve(et2)
    print("ifail: {:.0f}".format(ifail))
    print_diff(et1,nb1)

if __name__=='__main__':

    # Update to reflect where you keep these kernels
    kdir = '../../kernels/'
    lsk = kdir+"naif0010.tls"
    spk = kdir+"ura111.bsp"

    spice.furnsh(lsk)
    spice.furnsh(spk)

    # Pull out the temporal coverage of the kernel
    et0,ete = spice.spkcov(spk,799)[0]

    print("Running forward for length of kernel:")
    run_uranus(et0,ete)

    print("\nRunning backward for length of kernel:")
    run_uranus(ete,et0)

    spice.kclear()
