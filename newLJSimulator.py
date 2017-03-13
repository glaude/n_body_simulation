"""Slightly different notation from yours. I had ListMethods import as method. Apart from that, nothing else."""

import newListMethods as method
import MDUtilities as md
from Particle3D import Particle3D
import numpy as np
import time

def MDEngine(listObjects, r_cut, boxSide, rangeParticles, density, outfileTraj, outfileRDF):
    numstep = 100
    mirrorImages = method.mirror_image(listObjects, boxSide)
    force, allPD = method.new_force(listObjects, r_cut, boxSide, rangeParticles, mirrorImages)


    # Set up bins and histogram for rdf.
    histPlot = np.zeros(25)
    binRange = [0]
    n = 0
    while n < 5.0:
        binRange.append(0.2 + n)
        n += 0.2

    # Engine start!!!
    for i in range(numstep):

        # Time skipping
        if i % 5 == 0: # Timestep skipping. Best to set to 5
            outfileTraj.write("{}\nPoint = {}\n".format(len(listObjects), i))
            method.traj_output(outfileTraj, listObjects)

        # Calculate RDF
        if i >= 5: # Only record after equilibration. Adjust accordingly for larger time steps. Ie. 1000, i>=50.
            histData = np.histogram(allPD, bins = binRange)
            histPlot = histData[0] + histPlot

        # Pos, velocity, force updates
        method.posUpdatePBC(listObjects, force, boxSide)
        newMirror = method.mirror_image(listObjects, boxSide)
        force_new, allPD_new = method.new_force(listObjects, r_cut, boxSide, rangeParticles, newMirror)
        force_avg = 0.5 * (force + force_new)
        method.velUpdate(listObjects, force_avg, rangeParticles)
        force, allPD = force_new, allPD_new


    # RDF Normalisation and writing
    normHist = method.normaliser(histPlot, binRange, density, numstep, len(rangeParticles))
    for i in range(len(histPlot)):
        outfileRDF.write("{0}  {1}\n".format(binRange[i], normHist[i]))



def main():
    # Setup of lists, parameters, etc. etc.
    listObjects = [] # Contains list of Particle3D instances
    nParticles = int(raw_input("Input number of particles: ")) # Note: 108 atoms can fill a lattice fully.
    listObjects.extend([Particle3D(i, 0, 0, 0, 0, 0, 0, 1) for i in range(nParticles)])
    rangeParticles = range(nParticles)
    density = float(raw_input("Define density: "))
    temp = float(raw_input("Define temperature: "))
    box = md.setInitialPositions(density, listObjects) # Generates box dimension and set initial position
    md.setInitialVelocities(temp, listObjects)
    r_cut = float(raw_input("Define cutoff radius: "))

    # User input of outfile names. Deactivate if not using.
    outfileNameTraj = raw_input("Name of trajectory output file: ")  # This has to be .xyz
    outfileNameRDF = raw_input("Name of rdf file: ")
    outfile_traj = open(outfileNameTraj, "w")
    outfileRDF = open(outfileNameRDF, "w")

    # Starting the Engine!
    print ("Simulation started at" ,time.ctime()) # Measuring how long simulation is.
    MDEngine(listObjects, r_cut, box[0], rangeParticles, density, outfile_traj, outfileRDF)
    outfile_traj.close()
    outfileRDF.close()


main()
print("Simulation ended on", time.ctime()) # Measurement, again.
