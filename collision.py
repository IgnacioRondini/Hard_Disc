
# collision.py Compute the times collision and others
from objects import *
from math import sqrt
from functions import *


def get_collision_times_with_walls(particle, w):
    # It gives the time untill the collision
    # between the disc and the wall. It checks for instantaneous
    # collision that must be avoided
    if w.get_id() != particle.get_last_collision():
        y0 = (particle.get_position())[1]
        x0 = (particle.get_position())[0]
        a = (w.get_parameters())[0]
        b = (w.get_parameters())[1]

        vx = (particle.get_speed())[0]
        vy = (particle.get_speed())[1]

        t = (y0 - a*x0 - b) / (a*vx - vy) #.real

        if t > 0:
            return t
        else:
            return float('inf')
    else:
        return float('inf')


def get_collision_times_with_obstacles(particle, obs):
    # It gives the time untill the collision
    # between the disc and one of the obstacle. It checks for instantaneous
    # collision that must be avoided
    if obs.get_radius() <= 0:
        return float('inf')

    if obs.get_id() != particle.get_last_collision():
        x_bar = ((particle.get_position())[0] - (obs.get_position())[0])#.real
        y_bar = ((particle.get_position())[1] - (obs.get_position())[1])#.real
        vx = (particle.get_speed())[0]#.real
        vy = (particle.get_speed())[1]#.real
        rho = obs.get_radius()#.real
        if ((x_bar*vx + y_bar*vy)**2 -(vx**2+vy**2)*(x_bar**2+y_bar**2 - rho**2)) >= 0:
            if vx**2 + vy**2 == 0:
                return float('inf')
            new_time = -(x_bar*vx + y_bar*vy)-sqrt( abs((x_bar*vx + y_bar*vy)**2 -(vx**2+vy**2)*(x_bar**2+y_bar**2 - rho**2)))
            new_time = new_time/(vx**2 + vy**2)
            if new_time >= 0:
                return new_time
            else:
                return float('inf')
        else:
            return float('inf')

    return float('inf')


def get_collision_times(particle, collision):
    # Create a list with all the times until the next collision
    # for each wall and obstacle
    times = []
    for coll in collision:
        if coll.get_id() > 0:
            times.append(abs(get_collision_times_with_obstacles(particle, coll)))
        else:
            times.append(abs(get_collision_times_with_walls(particle, coll)))
    return times


def get_collision_AB_times(A, B, L):
    # Computes the time until the next collision between the two discs
    # Because we have boundary conditions, we need to check for collision
    # through the boundaries
    times = []
    x_bar = A.get_position()[0] - B.get_position()[0]
    y_bar = A.get_position()[1] - B.get_position()[1]
    vx_bar = A.get_speed()[0] - B.get_speed()[0]
    vy_bar = A.get_speed()[1] - B.get_speed()[1]

    H = L * sqrt(3.) / 2.

    Displacement = ((0., 0.),(3.*L/2., H),(0., 2.*H), (-3.*L/2., H), (-3.*L/2., -H), (0., -2*H), (3.*L/2., -H))
    for vect in Displacement:
        index = Displacement.index(vect)

        if index != A.get_AB_index() or index != B.get_AB_index():

            x_bar = A.get_position()[0] + vect[0] - B.get_position()[0]
            y_bar = A.get_position()[1] + vect[1] - B.get_position()[1]
            if (x_bar**2 + y_bar**2) <= ((A.get_radius() + B.get_radius())**2):
                print "disks are overlapping"
                exit()

            else:
                delta = (x_bar*vx_bar + y_bar*vy_bar)**2
                delta = delta -(vx_bar**2 + vy_bar**2)*(x_bar**2 + y_bar**2 - (A.get_radius() + B.get_radius())**2)
                if (x_bar*vx_bar + y_bar*vy_bar) >= 0:
                    times.append(float('inf'))
                elif delta < 0:
                    times.append(float('inf'))
                else:
                 # delta > 0
                    new_time = -(x_bar*vx_bar + y_bar*vy_bar + sqrt(abs(delta))) / (vx_bar**2 + vy_bar**2)
                    if new_time > 0.:
                        times.append(new_time)
                    else:
                        times.append(float('inf'))
        else:
            times.append(float('inf'))
    return times


def Collide(particles, obstacles, walls, lattice, current_time):
    # particles = list of particle
    # concat obstacles and walls
    objects = obstacles + walls

    times_A = get_collision_times(particles[0], objects)
    times_B = get_collision_times(particles[1], objects)
    times_AB = get_collision_AB_times(particles[0],particles[1], walls[0].get_length())

    t_ca = min(times_A)
    t_cb = min(times_B)
    t_AB = min(times_AB)

    L = walls[0].get_length()
    H = walls[0].get_length() * sqrt(3.) / 2.

    Displacement = ((0., 0.),(3.*L/2., H),(0., 2.*H), (-3.*L/2., H), (-3.*L/2., -H), (0., -2*H), (3.*L/2., -H))
    Displacement_lattice = ((0, 0), (1, 0), (1, -1), (0, -1), (-1, 0), (-1, +1), (0, 1))

    for v in Displacement:
        x_bar = particles[0].get_position()[0] + v[0] - particles[1].get_position()[0]
        y_bar = particles[0].get_position()[1] + v[1] - particles[1].get_position()[1]

    if(t_AB < t_ca and t_AB < t_cb):
        # it's a collision between particles
        # advance to collision position

        particles[0].free_flight(t_AB)
        particles[1].free_flight(t_AB)

        index = times_AB.index(t_AB)
        # change speed
        particles[0].Collide_with_particle(particles[1], Displacement[index])

        particles[0].set_last_collision(0)
        particles[1].set_last_collision(0)
        particles[0].set_last_AB_index(index)
        particles[1].set_last_AB_index(index)

        lattice.update_collision(Displacement_lattice[index], particles[0].get_energy_exchanged(), current_time + t_AB)
        return t_AB

    else:
        if t_ca < t_cb:
            # it is a collision of 100 = A-particle
            index = times_A.index(t_ca)
            collision_id = (objects[index]).get_id()

            particles[0].free_flight(t_ca)
            particles[1].free_flight(t_ca)

            particles[0].set_last_AB_index(-1)
            if collision_id < 0:
                # if it was    a wall
                particles[0].set_last_collision(objects[index].get_periodic_wall_id())

                vect_index = Displacement_lattice.index(objects[index].get_lattice_deplacement())
                vect = Displacement[vect_index]

                particles[0].pass_trough_boundary(vect, current_time, objects[index].get_lattice_deplacement())

                lattice.update_move(particles[0].get_id(), objects[index].get_lattice_deplacement())

            elif collision_id > 0:
                # if it was an obstacle
                particles[0].set_last_collision(collision_id)
                normal = objects[index].get_normal(particles[0].get_position())
                particles[0].change_speed_collision(normal)
            else:
                print "no collision ERROR"
                exit()
            return t_ca
        else:
            if t_ca <= t_cb:
                print "error in times of collision "
                exit()
            # it is a collision of 200 = B-particle
            index = times_B.index(t_cb)
            collision_id = (objects[index]).get_id()
            particles[0].free_flight(t_cb)

            particles[1].free_flight(t_cb)

            particles[1].set_last_AB_index(-1)
            if collision_id < 0:

                vect_index = Displacement_lattice.index(objects[index].get_lattice_deplacement())
                vect = Displacement[vect_index]

                particles[1].set_last_collision(objects[index].get_periodic_wall_id())
                particles[1].pass_trough_boundary(vect, current_time, objects[index].get_lattice_deplacement())
                lattice.update_move(particles[1].get_id(), objects[index].get_lattice_deplacement())

            elif collision_id > 0:
                # if it was an obstacle
                particles[1].set_last_collision(collision_id)
                normal = objects[index].get_normal(particles[1].get_position())
                particles[1].change_speed_collision(normal)
                print "error collision with obstacle "
                exit()
            else:
                print "no collision ERROR"
                exit()
            return t_cb
