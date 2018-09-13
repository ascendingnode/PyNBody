from __future__ import print_function
from PyNBody import NBody
import numpy as np

def run_test(e):

    # Parameters
    norb = 100.
    P = 2.*np.pi

    print('Testing a two-body system with e = {:.3f} for {:.0f} orbits'.format(e,norb))

    # Set up the NBody system
    rp = 1.+e
    vp = np.sqrt((1.-e)/(1.+e))
    nb1 = NBody()
    nb1.t = 0
    nb1.verbose = True
    nb1.add_object("Big"  ,0.9,0.1,[0,0,0],[0,0,0])
    nb1.add_object("Small",0.1,0.1,[rp,0,0],[0,vp,0])
    nb1.maxdist = np.inf
    nb1.move2bary()

    # Use .copy() to get a deep copy of the object
    nb0 = nb1.copy()

    # Save the initial energy and momentum
    L0,H0 = (nb1.momentum(),nb1.hamiltonian())

    # Integrate the system for 100 orbits
    ifail = nb1.evolve(P*norb)

    # Find the difference in position and velocity
    dr = (nb1.get_position(0)-nb0.get_position(0)) + (nb1.get_position(1)-nb0.get_position(1))
    dr = np.sqrt(np.mean(np.square(dr)))
    dv = (nb1.get_velocity(0)-nb0.get_velocity(0)) + (nb1.get_velocity(1)-nb0.get_velocity(1))
    dv = np.sqrt(np.mean(np.square(dv)))

    # Find the difference in energy and momentum
    L1,H1 = (nb1.momentum(),nb1.hamiltonian())
    dH = abs((H1-H0))
    dL = np.linalg.norm(L1-L0)

    # Output
    print('\tFailure state:              {:.0f}'.format(ifail))
    print('\tNumber of orbits completed: {:.3f}'.format(nb1.t/P))
    print('\tAbsolute error in energy:   {:e}'.format(dH))
    print('\tAbsolute error in momentum: {:e}'.format(dL))
    if ifail==0:
        print('\tAbsolute error in position  {:e}'.format(dr))
        print('\tAbsolute error in velocity  {:e}'.format(dv))

if __name__=='__main__':

    # Run an eccentric but stable system
    run_test(0.5)

    print()

    # Run a system so eccentric that it impacts before reaching periapse
    run_test(0.9)
