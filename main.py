

# main

from functions import *
from objects import *
from collision import *
from geometry import *
from Lattice import *
from math import sqrt
import random
import numpy as np


def set_environment(walls, discs, l, radius):
    # place the discs first in the corners of a hexagone
    discs.append(obstacle([-l/2., - l*sqrt(3.)/2.], radius, 10))
    discs.append(obstacle([+l/2., - l*sqrt(3.)/2.], radius, 11))
    discs.append(obstacle([l, 0.], radius, 12))

    discs.append(obstacle([l/2., l*sqrt(3.)/2.], radius, 13))
    discs.append(obstacle([-l/2., l*sqrt(3.)/2.], radius, 14))
    discs.append(obstacle([-l, 0.], radius, 15))

    walls.append(wall(discs[0], discs[1], -10))
    walls.append(wall(discs[1], discs[2], -11))
    walls.append(wall(discs[2], discs[3], -12))

    walls.append(wall(discs[3], discs[4], -13))
    walls.append(wall(discs[4], discs[5], -14))
    walls.append(wall(discs[5], discs[0], -15))

    walls[0].set_periodic_wall_id(-13)
    walls[1].set_periodic_wall_id(-14)
    walls[2].set_periodic_wall_id(-15)

    walls[3].set_periodic_wall_id(-10)
    walls[4].set_periodic_wall_id(-11)
    walls[5].set_periodic_wall_id(-12)


def print_walls(walls):
    for w in walls:
        print w.get_position()[0], w.get_position()[1], w.get_id(), w.get_periodic_wall_id(), w.get_parameters()


def print_obstacles(obs):
    for o in obs:
        print o.get_position()[0], o.get_position()[1], o.get_radius(), o.get_id()


def print_particles(parts):
    for particle in parts:
        print particle.get_position(),particle.get_speed(), particle.get_radius(), particle.get_id()


def set_particles(length, radius, radius_a, radius_b, discs):
    # Generate two particles which are not inside any object
    # the total  conserved momentum is momentum P = (1,1)
    # the labels are A = 1 or B = -1
    # the id is an integer (100 or 200 in this case)
    particles = []
    pol = []
    for vertex in discs:
        pol.append(vertex.get_position())

    angle_speed = random.random() * 2. * pi
    # relative speed is unitary
    # urel = (2V*cos(), 2*V*sin()-
    # ||urel|| = sqrt(4*V^2cos^2+4*V^2*sin^2) = 2*V
    V = 1. / 2.

    # set first disc
    while True:
        # set random position inside the outer circle of the hexagon
        l = random.random() * length
        angle = random.random() * 2. * pi
        # check if it is inside the polygon
        if is_in_polygon(pol, (l*cos(angle), l*sin(angle))):
            particles.append(particle([l*cos(angle), l*sin(angle)], [V*cos(angle_speed), V*sin(angle_speed)], radius_a, 1, 100))
            break

    # set second disc
    while True:
        # set random position inside the outer circle of the hexagon
        l = random.random() * length
        angle = random.random() * 2. * pi
        # check if it is inside the polygon
        if is_in_polygon(pol, (l*cos(angle), l*sin(angle))):
            # check if it is not superposed with the previous disc
            if not are_superposed_in_all_cell(particles[0].get_position(),radius_a,[l*cos(angle), l*sin(angle)],radius_b, length):
                particles.append(particle([l*cos(angle),l*sin(angle)], [-V*cos(angle_speed),-V*sin(angle_speed)], radius_b, -1, 200))
                break

    particles[0].set_initial_relative_pos_coll(particles[1])
    particles[1].set_initial_relative_pos_coll(particles[0])

    return particles


def update_histogram(List_of_color, Histo_A, N, dt):
    for p in List_of_color:
        i = p[0][0]
        i = i + N
        Histo_A[i] = Histo_A[i] + dt


def main():
    random.seed(987)
    current_time = 0.

    l = 1. / sqrt(3)
    radius_of_obstacles = 0.
    Max_Time = 8000

    count = 0
    walls = []
    discs = []
    pol = []
    # SET THE UNIT CELL: A Hexagonal cell with periodic boundary conditions
    # The are no solid walls.
    # Solid fixed disc can be added to the vertex

    set_environment(walls, discs, l, radius_of_obstacles)
    for vertex in discs:
        pol.append(vertex.get_position())

    # SET INITIAL CONDITION FOR THE DISCS
    radius = 0.235
    radius_A = radius
    radius_B = radius
    particles = set_particles(l, radius_of_obstacles, radius_A, radius_B, discs)

    # SET THE LATTICE. It does correctly work for a one-dimensional lattice
    # for a 2-D lattice, consistency must be checked
    # First argument: Dimension of he lattice,
    # Second argument: Probability of adding a colored disc from the boundaries
    #                  (A from left and B from right)
    # Third argument: Probability of reaction
    N = 0
    Hexagonal_Lattice = Lattice((N, 0), (0., 0.), 1 - 0.00151)
    # Set number of colored disc in the initial configuration
    # (randomly distributed)
    Hexagonal_Lattice.set_initial_conditions(0, 0)
    Hexagonal_Lattice.get_density_distribution("0")

    # ---------------MAIN LOOP ---------------- #
    while current_time <= Max_Time:
        # Collide() iterates the system to the next event
        # (either collisions or a disc passing across a boundary )
        # it returns the time to the next event
        dt = Collide(particles, discs, walls, Hexagonal_Lattice, current_time)

        current_time = current_time + dt

        # ------------COMPUTE OBSERVABLES--------------#

        # ------------END COMPUTE OBSERVABLES----------#

        count += 1
    # ------------END MAIN LOOP ---------------- #

if __name__ == "__main__":
    main()
