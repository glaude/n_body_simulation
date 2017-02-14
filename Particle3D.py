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
        return ("{} {} {} {}".format(self.label, self.position[0], self.position[1], self.position[2]))

    def kinetic_energy(self):
        """
        Calculates the kinetic energy of a particle.
        :return: Kinetic energy value.
        """
        return 0.5 * self.mass * (np.linalg.norm(self.velocity) ** 2)

    #  Potential Energy
    def potential_energy(self):
        """
        Calculates potential energy of a particle.
        :return: Potential energy value.
        """
        G = 1  # As defined in exercise 3, section 2.2.
        return G * self.mass * (-1 / np.linalg.norm(self.position))

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

    # Parameter Reader
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

    # Vector Separation
    @staticmethod
    def vector_split(p_1, p_2):
        """
        Calculates the vector separation between two particles.
        :param p_1: Particle 3D instance
        :param p_2: Another particle 3D instance.
        :return: The relative vector separation between 2 Particle3D instances.
        """
        vector_separation = p_1.position - p_2.position # changed from p_2.position - p_1.position
        return vector_separation

    @staticmethod
    def inter_force(r_cut, p_1, p_2):
        vect_mag = np.linalg.norm(Particle3D.vector_split(p_1, p_2))
        if vect_mag <= r_cut:
            force = np.multiply(48 * ((1 / (vect_mag ** 14)) - 1 / (2 * vect_mag ** 8)),
                                Particle3D.vector_split(p_1, p_2)) # Edited a mistake. Changed the '+' to a '-'.
            return force
        else:
            force = np.zeros(3) # changed from 0 to np.zeros(3)
            return force
