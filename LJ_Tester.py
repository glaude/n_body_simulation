import MDUtilities as md
from Particle3D import Particle3D
import Method_Test as method

def loop(listObjects, listNumber, outfile):
    numstep = 300
    force = method.force_sum(listObjects) * 1 / len(listObjects)

    # Update everything with this loop!
    for k in range(numstep):
        # Output file
        listNumber.append(listObjects)
        outfile.write("{}\nPoint = {}\n".format(len(listObjects), k))
        method.traj_output(outfile, listNumber[k], len(listObjects))

        # Update all shits
        method.posUpdate(listObjects, len(listObjects), force)
        force_new = method.force_sum(listObjects) * 1 / len(listObjects)
        force_avg = 0.5 * (force + force_new)
        method.velUpdate(listObjects, len(listObjects), force_avg)
        force = force_new

def main():
    # Setup of lists, parameters, etc. etc.
    instances = []
    listNumber = []
    nParticles = int(raw_input("Input number of particles: "))
    file_name = raw_input("Name of parameters file: ")
    file_handle = open(file_name, "r")
    p_1 = Particle3D.parameter_reader(file_handle)
    density = float(raw_input("Define density: "))
    temp = float(raw_input("Define temperature: "))
    outfileName = raw_input("Name of output file: ")
    outfile = open(outfileName, "w")
    for i in range(nParticles):
        instances.append(
            Particle3D(i, p_1.position[0], p_1.position[1], p_1.position[2], p_1.velocity[0], p_1.velocity[1],
                       p_1.velocity[2], p_1.mass))

    # Set initial positions and velocity
    md.setInitialPositions(density, instances)
    md.setInitialVelocities(temp, instances)

    # Generate values and write to file
    loop(instances, listNumber, outfile)


main()


