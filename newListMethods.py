from Particle3D import Particle3D
import numpy as np
from math import pi

dt = 0.01 # Adjust accordingly. numstep = 1000, dt = 0.001. numstep = 10000, dt = 0.0001 and so on.


def traj_output(outfile, listObjects):
    for j in range(len(listObjects)):
        outfile.write("{}\n".format(Particle3D.__str__(listObjects[j])))

def normaliser(histData, binData, density, timestep, nParticles):
    histData = [histData[i] * 1 / (4 * pi * density * (binData[i] ** 2) * (nParticles - 1) * 0.2 *
               (timestep - 5)) for i in range(len(histData)) if i != 0]
    histData.insert(0, 0) # Inserts 0 at beginning. Prevents a divide by zero error due to 1/r**2 term.
    return histData

def mirror_image(listObjects, boxDim):
    """
    This is basically the same as your code. But now I will only call
    positions rather than an entirely new object.
    """
    mirrorPos = [] # No need to call for this in the main code.

    for i in range(len(listObjects)):
        pos = listObjects[i].position
        #if np.any(pos) > (boxDim - r_cut) or np.any(pos) < r_cut:
        if pos[0] > boxDim/2:
            newPos = pos - np.array([boxDim, 0, 0])
            mirrorPos.append(newPos)
        if pos[0] < boxDim/2:
            newPos = pos + np.array([boxDim, 0, 0])
            mirrorPos.append(newPos)
        if pos[1] > boxDim/2:
            newPos = pos - np.array([0, boxDim, 0])
            mirrorPos.append(newPos)
        if pos[1] < boxDim/2:
            newPos = pos + np.array([0, boxDim, 0])
            mirrorPos.append(newPos)
        if pos[2] > boxDim/2:
            newPos = pos - np.array([0, 0, boxDim])
            mirrorPos.append(newPos)
        if pos[2] < boxDim/2:
            newPos = pos + np.array([0, 0, boxDim])
            mirrorPos.append(newPos)

    return mirrorPos

def new_force(listObjects, r_cut, boxSide, rangeParticles, mirrorImages):
    """ Using triangular matrices!"""
    force_sum = []
    force_vectors = []
    allPD = []
    for i in rangeParticles:
        force_upper = np.zeros((len(listObjects) - 1, 3)) # Upper triangle
        force_lower = np.zeros((len(listObjects) - 1, 3)) # Lower triangle.

        for j in rangeParticles:
            if j > i:
                # Upper triangle value generation.
                separation = Particle3D.vector_split(listObjects[i], listObjects[j])
                sepMag = np.linalg.norm(separation)
                if sepMag <= boxSide/2:
                    allPD += [sepMag, sepMag] # To fill box fully
                force = Particle3D.inter_force(separation)
                force_upper[j - 1] = force

            if j < i:
                # Take values from upper triangle and multiply by - 1.
                force_lower[j] = -1 * force_vectors[j][i - 1]
        force_full = np.add(force_upper, force_lower) # Evaluate full matrix.
        force_vectors.append(force_full)

        # Including minimum image convention
        mirrorForces = []
        for k in range(len(mirrorImages)):
            separation = listObjects[i].position - mirrorImages[k]
            sepMag = np.linalg.norm(separation)
            if sepMag <= boxSide/2:
                allPD.append(sepMag)
            if sepMag < r_cut:
                force = Particle3D.inter_force(separation)
                mirrorForces.append(force)

        force_sum.append(sum(force_full) + sum(mirrorForces))

    return np.asarray(force_sum), allPD

def velUpdate(listObjects, sum_forces, rangeParticles) :
    for i in range(len(listObjects)) :
        Particle3D.leap_velocity(listObjects[i], dt, sum_forces[i])


def posUpdatePBC(listObjects, sum_forces, boxSide):
    """Merged posUpdate into pbc."""
    for i in range(len(listObjects)):
        Particle3D.leap_pos2nd(listObjects[i], dt, sum_forces[i])
        np.putmask(listObjects[i].position, listObjects[i].position > boxSide, listObjects[i].position - boxSide)
        np.putmask(listObjects[i].position, listObjects[i].position < 0, listObjects[i].position + boxSide)