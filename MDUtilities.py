"""
CMod Project B: auxiliary MD methods
"""

import random
import numpy as np


def setInitialPositions(rho, particles):
    # Determine number of particles
    nAtoms = len(particles)

    # Set box dimensions
    boxSize = (nAtoms / rho) ** (1. / 3.)

    # Number or particles in each direction
    nDim = int(float((nAtoms - 1) / 4.0) ** (1. / 3.)) + 1

    # Give warning if fcc lattice will not be fully occupied
    if 4 * nDim ** 3 != nAtoms:
        print "Atoms will not fill a fcc lattice completely.\n"

    # Separation between particles
    delta = boxSize / nDim

    # Set particle positions
    iAtom = 0
    for ix in range(nDim):
        for iy in range(nDim):
            for iz in range(nDim):
                if iAtom < nAtoms:
                    pos = np.array([ix * delta, iy * delta, iz * delta])
                    particles[iAtom].position = pos
                    iAtom += 1
                if iAtom < nAtoms:
                    pos = np.array([(ix + 0.5) * delta, (iy + 0.5) * delta, iz * delta])
                    particles[iAtom].position = pos
                    iAtom += 1
                if iAtom < nAtoms:
                    pos = np.array([(ix + 0.5) * delta, iy * delta, (iz + 0.5) * delta])
                    particles[iAtom].position = pos
                    iAtom += 1
                if iAtom < nAtoms:
                    pos = np.array([ix * delta, (iy + 0.5) * delta, (iz + 0.5) * delta])
                    particles[iAtom].position = pos
                    iAtom += 1

    # Some output
    print "{0:d} atoms placed on a face-centered cubic lattice.\n".format(nAtoms)
    print "Box dimensions: {0:f} {0:f} {0:f}\n".format(boxSize)

    # Return the box size as Vector3D object
    return np.array([boxSize, boxSize, boxSize])


def setInitialVelocities(temp, particles):
    # Determine number of particles
    nAtoms = len(particles)

    # Zero the accumulators
    xv0 = 0.0
    yv0 = 0.0
    zv0 = 0.0
    vsq = 0.0

    # Loop over particles, set velocities
    for i in range(nAtoms):
        # Random inital velocities
        xvt = random.random() - 0.5
        yvt = random.random() - 0.5
        zvt = random.random() - 0.5

        particles[i].velocity = np.array([xvt, yvt, zvt])

        # Add to total velocity
        xv0 += xvt
        yv0 += yvt
        zv0 += zvt
        vsq += xvt ** 2 + yvt ** 2 + zvt ** 2

    # Centre-of-mass motion
    xv0 /= nAtoms
    yv0 /= nAtoms
    zv0 /= nAtoms

    # Boltzmann factor
    kB = (3 * nAtoms * temp / vsq) ** (1. / 2.)

    # Zero the probe accumulators
    xv0Tot = 0.0
    yv0Tot = 0.0
    zv0Tot = 0.0
    v0sq = 0.0

    # Rescale all velocities
    for i in range(nAtoms):
        vtemp = particles[i].velocity
        xvt = kB * (vtemp[0] - xv0)
        yvt = kB * (vtemp[1] - yv0)
        zvt = kB * (vtemp[2] - zv0)

        particles[i].velocity = np.array([xvt, yvt, zvt])

        xv0Tot += xvt
        yv0Tot += yvt
        zv0Tot += zvt
        v0sq += xvt ** 2 + yvt ** 2 + zvt ** 2

    # Output
    print "Temperature = {0:f}\n".format(temp)
    print "Centre-of-mass velocity = {0:f} {1:f} {2:f}\n".format(xv0Tot / nAtoms, yv0Tot / nAtoms, zv0Tot / nAtoms)
