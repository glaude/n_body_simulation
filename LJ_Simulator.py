"""
Pull requests before committing to main branch!

This can be executed to generate initial positions as of now.
"""
import ListMethods as method
import MDUtilities as md
from Particle3D import Particle3D
import time

def loop():
    """Nothing so far."""
    numstep = 300
    for k in range(numstep):
        """Nothing yet. All updates to position, etc. go in here"""

def main(): #
    """
    Currently only asks for trajectory output file. The rest can be added later.
    :return:
    """
    # Setup of lists, parameters, etc. etc.
    listObjects = [] # Contains list of Particle3D instances
    listNumber = [] # Particle3D instance behaviour over time.
    nParticles = int(raw_input("Input number of particles: ")) # Note: 108 atoms can fill a lattice fully.
    listObjects.extend([Particle3D(i, 0, 0, 0, 0, 0, 0, 1) for i in range(nParticles)]) # Generates nParticles at the edge of the lattice.
    density = float(raw_input("Define density: "))
    temp = float(raw_input("Define temperature: "))
    r_cut = float(raw_input("Define cutoff radius: "))
    outfileNameTraj = raw_input("Name of output file: ") # This has to be .xyz
    outfile_traj = open(outfileNameTraj, "w")


    # Set initial position, velocity and box dimensions.
    box = md.setInitialPositions(density, listObjects) # Generates box dimension, sets initial positions.
    md.setInitialVelocities(temp, listObjects)
    
    method.traj_output(outfile_traj, listObjects) # Example call to write to file. Use this line to confirm crystal structure (fcc)
                                                  # The line above must be inside loop() and not in the main function.
    # Generate values and write to file
    # loop(instances, listNumber, outfile_traj, box, (more shit here) This code will be executed once the loop function
    # meets the minimum progress required to execute a simulation.
    outfile_traj.close()
    
    # Useful for testing: recording start and end time
    print ("Simulation started at" ,time.ctime())

main()
print("Simulation ended on", time.ctime()) # Notifies when the simulation has ended.

