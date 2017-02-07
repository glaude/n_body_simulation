from Particle3D import Particle3D
import numpy as np

"""
Pull requests before committing to main code!

In order to successfully simulate the program, the minimum requirement for completion here is the following:
(i) Force method
(ii) Velocity update
(iii) Position update
(iv) Periodic boundary conditions
(v) Method to continuously write positions to file.

So best to split the project into 2 phases, with phase 1 being the above. Then phase 2 shall be the following:
(i) Energy (***)
(ii) Mean squared displacement
(iii) Radial distribution function (***)
(iv) Whatever we may have missed

*** = Methods we can use to ensure that the simulation is running correctly, of course apart from visual tracking of particles.
"""

def traj_output(outfile, listObjects):
    for j in range(len(listObjects)):
        outfile.write("{}\n".format(Particle3D.__str__(listObjects[j])))

def force_sum(listObjects, r_cut): # Issue: 3rd particle not evaluated?
    sum_forces = []
    for i in range(len(listObjects)):
        list_forces = []
        for j in range(len(listObjects)):
            if i != j:
                list_forces.append(Particle3D.inter_force(r_cut, listObjects[i], listObjects[j]))
            else:
                list_forces.append(np.zeros(3))
        sum_forces.append(sum(list_forces))
    return sum_forces
