"""
Describes a 3D particle, methods to update particle position, velocity, and calculation of energy values.
"""
import numpy as np


class Particle3D(object):
    def __init__(self, label, x_pos, y_pos, z_pos, vel_x, vel_y, vel_z, mass):
        """
        Initialises a 3D particle instance.
        :param label: particle instance label.
        :param x_pos: x-coordinates of particle position
        :param y_pos: y-coordinates of particle position
        :param z_pos: z-coordinates of particle position
        :param vel_x: velocity of particle in the x-direction
        :param vel_y: velocity of particle in the y-direction
        :param vel_z: velocity of particle in the z-direction
        :param mass: mass of particle
        """
        self.label = label
        self.position = np.array([x_pos, y_pos, z_pos])
        self.velocity = np.array([vel_x, vel_y, vel_z])
        self.mass = mass

    def __str__(self):
        """
        Provides output of 3D particle parameters as a string
        :return: String output of position and mass of particle.
        """
        return ("{0} {1:.5f} {2:.5f} {3:.5f}".format(self.label, self.position[0], self.position[1], self.position[2]))

    def kinetic_energy(self):
        """
        Calculates the kinetic energy of a particle.
        :return: Kinetic energy value.
        """
        return 0.5 * self.mass * (np.linalg.norm(self.velocity) ** 2)

    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        Updates velocity of particle for a given force vector and time step.
        :param dt: Time step
        :param force: Force experienced by particle due to gravitational interactions.
        """
        self.velocity = self.velocity + (dt * force) / self.mass  # a = dv/dt = (1/m)F

    def leap_pos1st(self, dt):
        """
        Updates particle position for a given time step at first order.
        :param dt: Time step
        """
        self.position = self.position + (dt * self.velocity)

    def leap_pos2nd(self, dt, force):
        """
        Updates particle position for a given time step at 2nd order.
        :param dt: Time step.
        :param force: Force experienced by particle due to gravitational interactions.
        """
        self.position = self.position + (dt * self.velocity) + (0.5 * (dt ** 2) * (force / self.mass))


    # Vector Separation
    @staticmethod
    def vector_split(p_1, p_2):
        """
        Static method which calculates the vector separation between two particles.
        :param p_1: Particle 3D instance
        :param p_2: Another particle 3D instance.
        :return: The relative vector separation between 2 Particle3D instances.
        """
        vector_separation = p_1.position - p_2.position
        return vector_separation

    @staticmethod
    def inter_force(separation, sepMag):
        """
        Static method which calculates the force between two particles.
        :param separation: separation vector between two particles.
        :param sepMag: magnitude of separation vector.
        :return: force vector between two particles.
        """
        force = np.multiply(48 * ((1 / (sepMag ** 14)) - (1 / (2 * sepMag ** 8))), separation)
        return force
    
    @staticmethod
    def potential_energy(r_cut, sepMag):
        """
        Static method which calculates the potential energy between two particles.
        :param r_cut: Cutoff radius
        :param sepMag: separation vector magnitude
        :return: Potential energy if the separation vector magnitude is below r_cut, otherwise 'zero'.
        """
        if sepMag < r_cut:
            pEnergy = 4 * (1/(sepMag ** 12) - 1/(sepMag **6))
            return pEnergy
        else:
            return 0

    @staticmethod
    def parameter_reader(file_handle):
        """
        Static method which generates a Particle3D instance from file input.
        :param file_handle:
        :return: Label, position, velocity vectors and particle mass values read from file input.
        """
        output = file_handle.readline()
        tokens = output.split()
        label = str(tokens[0])
        x_pos = float(tokens[1])
        y_pos = float(tokens[2])
        z_pos = float(tokens[3])
        vel_x = float(tokens[4])
        vel_y = float(tokens[5])
        vel_z = float(tokens[6])
        mass_part = float(tokens[7])

        return Particle3D(label, x_pos, y_pos, z_pos, vel_x, vel_y, vel_z, mass_part)