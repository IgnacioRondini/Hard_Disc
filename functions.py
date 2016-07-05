
# functions

# some miscelaneous functions
from math import sqrt


def scalar_product(a, b):
    if len(a) != len(b):
        return -900000000

    product = 0
    for i in range(len(a)):
        product = product + a[i] * b[i]

    return product


def are_superposed_in_all_cell(center_a, radius_a, center_b, radius_b, L):

    R = (radius_a + radius_b)

    H = L * sqrt(3.) / 2.
    Displacement = ((0., 0.), (3.*L/2., H), (0., 2.*H), (-3.*L/2., H), (-3.*L/2., -H), (0., -2*H), (3.*L/2., -H))
    for vect in Displacement:
        dy = center_a[1] + vect[1] - center_b[1]
    #    dy = dy.real
        dx = center_a[0] + vect[0] - center_b[0]
    #    dx = dx.real
        if dx**2 + dy**2 <= R**2:
            return True

    return False


def are_superposed_simple(center_a, radius_a, center_b, radius_b):
    R = (radius_a + radius_b).real
    dy = center_a[1] - center_b[1]
    #  dy = dy.real
    dx = center_a[0] - center_b[0]
    #  dx = dx.real
    if dx**2 + dy**2 <= R**2:
        return True
    return False


def get_energy(Particles):
    energy = 0
    for p in Particles:
        for speed in p.get_speed():
            energy = energy + speed**2

    return energy


def get_lattice_deplacement(wall_id):
    if wall_id == -10:
        return [-1, +1]
    if wall_id == -11:
        return [0, 1]
    if wall_id == -12:
        return [1, 0]
    if wall_id == -13:
        return [1, -1]
    if wall_id == -14:
        return [0, -1]
    if wall_id == -15:
        return [-1, 0]
