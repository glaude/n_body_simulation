from Particle3D import Particle3D
import numpy as np
from math import pi

dt = 0.001 # Adjust accordingly. numstep = 1000, dt = 0.001. numstep = 10000, dt = 0.0001 and so on.

def traj_output(outfile, listObjects, rangeParticles):
    """
    A method which writes positions of all objects in listObjects to a file.
    :param outfile: Trajectory file
    :param listObjects: List containing Particle3D instances
    :param rangeParticles: Range of length of listObjects
    :return:
    """
    for j in rangeParticles:
        outfile.write("{}\n".format(Particle3D.__str__(listObjects[j])))

def normaliser(histData, binData, density, timestep, nParticles):
    """
    Normalises the histogram generated in order to yield a proper radial distribution function.
    :param histData: All pair distances binned as a histogram
    :param binData: A list containing all the 'bins'.
    :param density: Reduced density
    :param timestep: Length of simulation
    :param nParticles: Number of particles
    :return: Normalised histogram data
    """
    histData = [histData[i] * 1 / (4 * pi * density * (binData[i] ** 2) * (nParticles) * 0.05 *
               (timestep - 200)) for i in range(len(histData)) if i != 0]
    histData.insert(0, 0) # Inserts 0 at beginning. Prevents a divide by zero error due to 1/r**2 term.
    return histData

def msd(initParticles, currentParticles, boxDim, nParticles, rangeParticles):
    """
    Calculates the mean squared displacement averaged over all particles
    :param initParticles: Reference particles at initial position
    :param currentParticles: Particles with position at time t
    :param boxDim: box dimensions
    :param nParticles: no. of particles
    :param rangeParticles: range of length of listObjects
    :return: Mean squared displacements averaged over all particles
    """

    posDif = np.asarray([mirror_image(initParticles[i], currentParticles[i], boxDim) for i in rangeParticles])
    posDifMag = np.linalg.norm(posDif, axis = 1)
    msd = np.sum(posDifMag ** 2) * 1/nParticles
    return msd

def mirror_image(p_1, p_2,boxDim):
    """
    Evaluates the separation of a particle and its mirror image, provided that the separation vector between two
    particles p_1 and p_2 are above half of the box dimensions.
    :param p_1: Particle 1
    :param p_2: Particle 2
    :param boxDim: Box dimension
    :return: The separation vector between particle 1 and the mirror image of particle 2.
    """
    separation = Particle3D.vector_split(p_1, p_2)
    pos = p_2.position
    if abs(separation[0]) > boxDim/2:
        if pos[0] > boxDim/2:
            separation[0] = separation[0] + boxDim
        else:
            separation[0] = separation[0] - boxDim
    if abs(separation[1]) > boxDim/2:
        if pos[1] > boxDim/2:
            separation[1] = separation[1] + boxDim
        else:
            separation[1] = separation[1] - boxDim
    if abs(separation[2]) > boxDim/2:
        if pos[2] > boxDim/2:
            separation[2] = separation[2] + boxDim
        else:
            separation[2] = separation[2] - boxDim
    return separation

def new_force(listObjects, r_cut, boxSide, rangeParticles):
    """
    Evaluates force, pair displacements and stores these in a list. The way the force calculation works is by evaluating
    a 'triangular matrix' (ie. abusing symmetry). Do check the 'Readme.txt' file provided for a visual representation of
    the triangular matrices and how this will be made to work. Also, a cutoff radius is implemented as in force
    calculations, pair distances scale as 1/r ** 14 and 1/r ** 8, meaning that high pair distance values will result
    in force values tending to 0. Thus, a cutoff radius is implemented to make sure that forces are not calculated past
    this value.
    :param listObjects: list of Particle3D instances
    :param r_cut: Cutoff radius
    :param boxSide: box dimensions
    :param rangeParticles: range of length of listObjects
    :return: list of forces and list of pair displacements.
    """
    force_sum = [] # List of forces experienced by each particle
    force_vectors = [] # List of force vectors for each particle (check 'Readme.txt')
    allPD = [] # List of pair displacement magnitudes, to be used in evaluating RDF.
    for i in rangeParticles:
        force_upper = np.zeros((len(listObjects) - 1, 3)) # Upper triangle matrix
        force_lower = np.zeros((len(listObjects) - 1, 3)) # Lower triangle matrix

        for j in rangeParticles:
            if j > i:
                # Upper triangle matrix value generation.
                new_separation = mirror_image(listObjects[i], listObjects[j], boxSide)
                sepMag = np.linalg.norm(new_separation)
                allPD += [sepMag, sepMag] # Separation magnitude appended twice as appending once does not represent a
                                          # full box.
                if sepMag < r_cut: # Imposing cutoff radius condition to minimise number of calculations needed.
                    force = Particle3D.inter_force(new_separation, sepMag)
                    force_upper[j - 1] = force
                else:
                    continue

            if j < i:
                # Takes values from upper triangle and multiplies them by - 1.
                force_lower[j] = -1 * force_vectors[j][i - 1]

        force_full = np.add(force_upper, force_lower) # Evaluate full force matrix by adding upper and lower matrices.
        force_vectors.append(force_full)
        force_sum.append(sum(force_full))

    return np.asarray(force_sum), allPD

def velUpdate(listObjects, sum_forces, rangeParticles) :
    """
    Updates the velocity of all Particle3D instances in listObjects
    :param listObjects: list containing Particle3D instances
    :param sum_forces: Forces experienced by each particle
    :param rangeParticles: range of the length of listObjects
    """
    for i in rangeParticles :
        Particle3D.leap_velocity(listObjects[i], dt, sum_forces[i])



def posUpdatePBC(listObjects, sum_forces, boxSide, rangeParticles):
    """
    Updates the position of the particles considering periodic boundary conditions (pbc).
    np.putmask() accounts for pbc, where conditions are such that if the position of the particle is outside of the
    box, the particle shall be placed back inside the box by adding or subtracting a box length to its position
    depending on whether the position of the particle is above the box length or below 0 in any coordinate (x, y, z).
    :param listObjects: list containing Particle3D instances
    :param sum_forces: forces experienced by each particle
    :param boxSide: box dimensions
    :param rangeParticles: range of length of listObjects
    """
    for i in rangeParticles:
        Particle3D.leap_pos2nd(listObjects[i], dt, sum_forces[i])
        np.putmask(listObjects[i].position, listObjects[i].position > boxSide, listObjects[i].position - boxSide)
        np.putmask(listObjects[i].position, listObjects[i].position < 0, listObjects[i].position + boxSide)
