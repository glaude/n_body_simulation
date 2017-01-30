# Method testing
import MDUtilities as md
from Particle3D import Particle3D
import numpy as np



def force_sum(listObjects):
    forces = []
    r_cut = 1.5
    for i in range(len(listObjects)):
        for j in range(len(listObjects)):
            if i != j:
                forces.append(Particle3D.inter_force(r_cut, listObjects[i], listObjects[j]))
    return sum(forces)

def traj_output(outfile, listObjects, nParticles):
    for j in range(nParticles):
        outfile.write("{}\n".format(Particle3D.__str__(listObjects[j])))

def velUpdate(listObjects, nParticles, force): # Also just making do with a symplectic Euler Integration for this.
    for l in range(nParticles):
        dt = 0.1
        listObjects[l].leap_velocity(dt, force)

def posUpdate(listObjects, nParticles, force): # Note: trying to make do with symplectic Euler integration for this test.
    for k in range(nParticles):
        dt = 0.1
        listObjects[k].leap_pos2nd(dt, force)

