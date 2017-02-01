"""
Pull requests before committing to main branch!

Do not attempt to execute if loop() is empty!
"""
import ListMethods as method
import MDUtilities as md
from Particle3D import Particle3D

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
    nParticles = int(raw_input("Input number of particles: "))
    file_name = raw_input("Name of parameters file: ") # Parameters.txt
    file_handle = open(file_name, "r")
    p_1 = Particle3D.parameter_reader(file_handle) # Opens parameter file containing one particle
    density = float(raw_input("Define density: "))
    temp = float(raw_input("Define temperature: "))
    outfileNameTraj = raw_input("Name of output file: ") # This has to be .xyz
    outfile_traj = open(outfileNameTraj, "w")
    for i in range(nParticles): # This loop generates 100 particles at origin.
        listObjects.append(
            Particle3D(i, p_1.position[0], p_1.position[1], p_1.position[2], p_1.velocity[0], p_1.velocity[1],
                       p_1.velocity[2], p_1.mass))

    # Box dimension
    box = md.setInitialPositions(density, listObjects) # Generates box dimension
    # Set initial positions and velocity
    md.setInitialPositions(density, listObjects)
    md.setInitialVelocities(temp, listObjects)

    # Generate values and write to file
    # loop(instances, listNumber, outfile_traj, box) This code will be executed once the loop function
    # meets the minimum progress required to execute a simulation.
    outfile_traj.close()