# geometry.py

from math import *

# Some methods needed for the initial conditions
# of the disc.
def is_point_in_triangle(triangle, point):
    v1 = [triangle[0][0], triangle[0][1]]
    v2 = [triangle[1][0], triangle[1][1]]
    v3 = [triangle[2][0], triangle[2][1]]
    p = [point[0], point[1]]

    # translate so V1 = (0,0)
    v2[0] = v2[0] - v1[0]
    v2[1] = v2[1] - v1[1]
    v3[0] = v3[0] - v1[0]
    v3[1] = v3[1] - v1[1]
    p[0] = p[0] - v1[0]
    p[1] = p[1] - v1[1]
    v1[0] = 0.
    v1[1] = 0.

    d = v2[0]*v3[1] - v3[0]*v2[1]
    wa = p[0]*(v2[1] - v3[1]) + p[1]*(v3[0] - v2[0]) + v2[0]*v3[1] - v3[0]*v2[1]
    wa = (wa / d)  # .real
    wb = p[0]*v3[1] - p[1]*v3[0]
    wb = (wb / d) # .real
    wc = p[1]*v2[0]-p[0]*v2[1]
    wc = (wc / d)  # .real
    if 0 <= min(wa, wb, wc) and max(wa, wb, wc) <= 1:
        return True
    return False


def rotate(point, angle):
    p_x = cos(angle)*point[0] - sin(angle)*point[1]
    p_y = sin(angle)*point[0] + cos(angle)*point[1]
    return[p_x, p_y]


def is_in_polygon(polygon, point):
    # polygon must be centered in zero
    for i in range(len(polygon) - 1):

        v2 = polygon[i]
        v3 = polygon[i + 1]

        if is_point_in_triangle([[0, 0], v2, v3], point):
            return True

    return is_point_in_triangle([[0, 0], polygon[-1], polygon[0]], point)
