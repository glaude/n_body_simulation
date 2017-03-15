"""Slightly different notation from yours. I had ListMethods import as method. Apart from that, nothing else."""

import newListMethods as method
import MDUtilities as md
from Particle3D import Particle3D
import numpy as np
import time
from copy import deepcopy


def MDEngine(listObjects, r_cut, boxSide, rangeParticles, density, nParticles, outfileTraj, outfileRDF, outfileMSD):
    """
    This is the 'engine' of the program. Basically, where everything from force, energy and observables are calculated.
    :param listObjects: a list containing all Particle3D instances.
    :param r_cut: Cutoff radius, in which any pair distance past this value, is not calculated.
    :param boxSide: Dimensions of the box.
    :param rangeParticles: The range of the length of the number of particles.
    :param density: Reduced density
    :param nParticles: No. of particles.
    :param outfileTraj: File containing list of trajectories (extension must be .xyz)
    :param outfileRDF: File containing radial distribution function data.
    :param outfileMSD: File containing mean squared displacement data.
    """
    numstep = 10000 # Simulation time
    force, allPD = method.new_force(listObjects, r_cut, boxSide, rangeParticles) # Calculates initial force and
                                                                                 # pair displacements.
    initParticles = deepcopy(listObjects) # Generates a reference point for mean squared displacement calculations.

    # Set up bins and histogram for rdf.
    histPlot = np.zeros(200) # Where the binned pair displacements shall be stored.
    binRange = [0] # Contains bin widths.
    n = 0
    while n < 10.0: # Defines the range of the bins
        binRange.append(0.05 + n)
        n += 0.05 # Bin size.

    # Starting the engine
    for i in range(numstep):

        # Time skipping
        if i % 5 == 0:  # Timestep skipping. Best to set to 5
            outfileTraj.write("{}\nPoint = {}\n".format(nParticles, i))
            method.traj_output(outfileTraj, listObjects, rangeParticles)

        # Position, velocity, force updates
        method.posUpdatePBC(listObjects, force, boxSide, rangeParticles)
        force_new, allPD_new = method.new_force(listObjects, r_cut, boxSide, rangeParticles)
        force_avg = 0.5 * (force + force_new)
        method.velUpdate(listObjects, force_avg, rangeParticles)

        # Calculate Radial Distribution Function (RDF)
        if i >= 200:  # Only record after equilibration.
            histData = np.histogram(allPD, bins=binRange)
            histPlot = np.add(histData[0], histPlot)
        # Calculation of Mean Squared Displacement and writing to file
        if i % 5 == 0:
            outfileMSD.write(
                "{0}  {1}\n".format(i, method.msd(initParticles, listObjects, boxSide, nParticles, rangeParticles)))

        # Update force and allPD
        force, allPD = force_new, allPD_new

    # RDF Normalisation and writing to file
    normHist = method.normaliser(histPlot, binRange, density, numstep, nParticles)
    for i in range(len(histPlot)):
        outfileRDF.write("{0}  {1}\n".format(binRange[i], normHist[i]))


def main():
    # Setup of lists, parameters, etc. etc.
    listObjects = []  # Contains list of Particle3D instances
    nParticles = int(raw_input("Input number of particles: "))  # Note: 108 atoms can fill a lattice fully.
    listObjects.extend([Particle3D(i, 0, 0, 0, 0, 0, 0, 1) for i in range(nParticles)])
    rangeParticles = range(nParticles)
    density = float(raw_input("Define density: "))
    temp = float(raw_input("Define temperature: "))
    box = md.setInitialPositions(density, listObjects)  # Generates box dimension and set initial position
    md.setInitialVelocities(temp, listObjects)
    r_cut = float(raw_input("Define cutoff radius: "))

    # User input of outfile names.
    outfileNameTraj = raw_input("Name of trajectory output file: ")  # This has to be .xyz
    outfileNameRDF = raw_input("Name of rdf file: ")
    outfileNameMSD = raw_input("Name of msd file: ")
    outfile_traj = open(outfileNameTraj, "w")
    outfileRDF = open(outfileNameRDF, "w")
    outfileMSD = open(outfileNameMSD, "w")

    # Starting the Engine!
    print ("Simulation started at", time.ctime())  # Measuring how long simulation is.
    MDEngine(listObjects, r_cut, box[0], rangeParticles, density, nParticles, outfile_traj, outfileRDF, outfileMSD)
    outfile_traj.close()
    outfileRDF.close()
    outfileMSD.close()


main()
print("Simulation ended on", time.ctime())  # Measurement, again.
