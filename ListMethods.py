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
        
def all_pos(listObjects, rangeParticles):
    """
    Stores the positions of all Particle3D instances in an array for mean square displacement calculations.
    :param listObjects: list containing Particle3D instances
    :param rangeParticles: The 'range' of the number of Particle3D instances
    :return: Returns an array of positions of all Particle3D instances.
    """
    allPos = np.array([listObjects[i].position for i in rangeParticles])
    return allPos

def msd(finalPos, initPos, rangeParticles):
    """
    Calculates the mean square displacement of all particles at time t
    :param finalPos: Position of all Particle3D instances at time t
    :param initPos: Initial position of all Particle3D instances (ie. at t = 0)
    :param rangeParticles:
    :return: Returns the mean square displacement at time t.
    """
    posVar = finalPos - initPos
    x = len(rangeParticles)
    posVarSq = [np.linalg.norm(posVar[i]) ** 2 for i in rangeParticles]
    msd = sum(posVarSq) / x
    return msd

def force_sum(listObjects, r_cut): 
    """
    The method currently does not take into account minimum image convention.
    I will need the minimum image convention method to update this.
    Lastly, if you are ok with making the new method regarding pair distances calculations,
    this can be significantly shortened.
    """
    sum_forces = [] # Contains list of force experienced by each particle in listObjects.
    for i in range(len(listObjects)):
        list_forces = [] # Pair forces stored here.
        for j in range(len(listObjects)):
            if i != j:
                separation = Particle3D.vector_split(listObjects[i], listObjects[j]) # Needed by pair_force as an argument.
                list_forces.append(Particle3D.pair_force(r_cut, separation)) # As above.
            else:
                list_forces.append(np.zeros(3))
        sum_forces.append(sum(list_forces)) # sum(list_forces) is the amount of force experienced by 1 particle. 
                                            # This calculation is repeated over listObjects.
    return sum_forces


def pbc(listObjects, boxSide, rangeParticles):
    """
    I just had to. Browsing through numpy documentation and saw this amazing method called putmask.
    This method subjects particle positions to periodic boundary conditions.
    :param listObjects: list of Particle3D instances.
    :param boxSide: Length of simulation box.
    :param rangeParticles: a 'range' over nParticles
    """
    for m in rangeParticles:
        """
        How it works: np.putmask('what you want to be updated', 'condition for update', 'the update operation')
        """
        pos = listObjects[m].position
        np.putmask(pos, pos > boxSide, pos - boxSide)
        np.putmask(pos, pos < 0, pos + boxSide)
