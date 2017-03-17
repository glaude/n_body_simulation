"""
This is the main program that enacts the simulation of particles interacting through a Lennard-Jones pair potential.

Written by Gabriel Laude (s1565983) and Mikolaj Roguski (s14*****).

Average simulation times for 10k time steps:
Solid (108 particles) : ~12 mins
Liquid/Gas (100 particles): 9 - 10 mins (Gas ~ 9min 12secs, Liquid ~ 10mins)
"""

import newListMethods as method
import MDUtilities as md
from Particle3D import Particle3D
import numpy as np
import time
from copy import deepcopy


def MDEngine(listObjects, r_cut, boxSide, rangeParticles, density, nParticles, timeskip, dt, outfileTraj, outfileRDF,
             outfileMSD):
    """
    This is the 'engine' of the program. Basically, where everything from force, energy and observables are calculated.
    :param listObjects: a list containing all Particle3D instances.
    :param r_cut: Cutoff radius, in which any pair distance past this value, is not calculated.
    :param boxSide: Dimensions of the box.
    :param rangeParticles: The range of the number of particles.
    :param density: Reduced density
    :param nParticles: No. of particles.
    :param timeskip: No. of time steps to skip for rdf recording (it should only be recorded after equilibration).
    :param dt = time step size
    :param outfileTraj: File containing list of trajectories (extension must be .xyz)
    :param outfileRDF: File containing radial distribution function data.
    :param outfileMSD: File containing mean squared displacement data.
    """
    numstep = 10000  # Simulation time
    force, allPD = method.new_force(listObjects, r_cut, boxSide, rangeParticles)  # Calculates initial force and
    # pair displacements (allPD).
    initParticles = deepcopy(listObjects)  # Generates a reference point for mean squared displacement calculations.

    # Set up bins and histogram for rdf.
    histPlot = np.zeros(200)  # Where the binned pair displacements shall be stored.
    binRange = [0]  # Contains bin widths.
    n = 0
    while n < 10.0:  # Defines the range of the bins
        binRange.append(0.05 + n)
        n += 0.05  # Bin size.

    # Starting the engine
    for i in range(numstep):

        # Time skipping
        if i % 5 == 0:  # Timestep skipping. Best to set to 5
            outfileTraj.write("{}\nPoint = {}\n".format(nParticles, i))
            method.traj_output(outfileTraj, listObjects, rangeParticles)

        # Position, velocity, force updates
        method.posUpdatePBC(listObjects, force, boxSide, rangeParticles, dt)
        force_new, allPD_new = method.new_force(listObjects, r_cut, boxSide, rangeParticles)
        force_avg = 0.5 * (force + force_new)
        method.velUpdate(listObjects, force_avg, rangeParticles, dt)

        # Calculate Radial Distribution Function (RDF)
        if i >= 200:  # Only record after equilibration.
            histData = np.histogram(allPD, bins=binRange)  # Binning of current pair displacements
            histPlot = np.add(histData[0], histPlot)  # Addition of current pair displacements to current overall total

        # Calculation of Mean Squared Displacement and writing to file
        if i % 5 == 0:
            outfileMSD.write(
                "{0}  {1:.5f}\n".format(i, method.msd(initParticles, listObjects, boxSide, nParticles, rangeParticles)))

        # Update force and allPD
        force, allPD = force_new, allPD_new

    # RDF Normalisation and writing to file
    normHist = method.normaliser(histPlot, binRange, density, numstep, nParticles, timeskip)
    for i in range(len(histPlot)):
        outfileRDF.write("{0}  {1}\n".format(binRange[i], normHist[i]))


def main():
    """
    User input regarding variables such as number of particles, density, temperature, etc. as well as file output names
    are taken here. These are then taken as arguments to start 'MDEngine' which will evaluate trajectories, observables,
    etc. over a given time (defined by numstep).
    """
    # Setup of lists, parameters, etc. etc.
    listObjects = []  # Contains list of Particle3D instances
    nParticles = int(raw_input("Input number of particles: "))  # Note: 108 atoms can fill a lattice fully.
    listObjects.extend([Particle3D(i, 0, 0, 0, 0, 0, 0, 1) for i in range(nParticles)])  # Particle generator loop
    rangeParticles = range(nParticles)
    density = float(raw_input("Define density: "))
    temp = float(raw_input("Define temperature: "))
    box = md.setInitialPositions(density, listObjects)  # Generates box dimension and sets initial positions
    md.setInitialVelocities(temp, listObjects)
    r_cut = float(raw_input("Define cutoff radius: "))
    dt = float(raw_input("Define time step size (suggested: 0.001): "))
    timeskip = int(raw_input("Define number of time steps to skip for RDF calculation: "))

    # User input of outfile names.
    outfileNameTraj = raw_input("Name of trajectory output file: ")  # File extension must be .xyz
    outfileNameRDF = raw_input("Name of rdf file: ")
    outfileNameMSD = raw_input("Name of msd file: ")
    outfileTraj = open(outfileNameTraj, "w")
    outfileRDF = open(outfileNameRDF, "w")
    outfileMSD = open(outfileNameMSD, "w")

    # Starting the Engine!
    print ("Simulation started at {}".format(time.ctime()))  # Record time at start of simulation.
    MDEngine(listObjects, r_cut, box[0], rangeParticles, density, nParticles, timeskip, dt, outfileTraj, outfileRDF,
             outfileMSD)

    # Close output files
    outfileTraj.close()
    outfileRDF.close()
    outfileMSD.close()


main()
print("Simulation ended on {}".format(time.ctime()))  # Record time at end of simulation.
