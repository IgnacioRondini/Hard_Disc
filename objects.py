
# from numpy import *
from math import *
from functions import scalar_product


# Class particle models a hard disc of radius r and mass = 1
class particle:
    def __init__(self, _position, _speed, _radius, _color, _id):
        self.position = [_position[0], _position[1]]
        self.initial_position = [_position[0], _position[1]]
        self.speed = _speed
        self.radius = _radius
        self.last_collision = 0
        self.id = _id
        self.color = _color
        self.last_ABcollision_index = -1
        self.real_position = [_position[0], _position[1]]
        self.number_of_collision = 0
        self.number_of_transitions = 0
        self.mass = 1.
        self.lattice_pos = [0, 0]
        self.supper_lattice_pos = [0, 0]


        name = "histogram_of_jumping_times" + str(_id)
        name = name+ ".dat"
        self.outfile_hist_jump_times = open(name, "w")
        name = str(_id) + "_pos.dat"
        self.outfile_pos = open(name, "w")
        # position_of_last_collision[0] = real position of the last pos x
        #
        self.position_of_last_collision = [_position[0], _position[1]]
        self.free_path = 0.

        self.relative_position_of_last_collision = [_position[0], _position[1]]
        self.relative_free_path = 0.

        self.time_of_last_entry = 0
        self.cage_time = 0

        self.energy_exchange = 0

    def get_initial_position(self):
        return self.initial_position

    def get_position(self):
        return self.position

    def get_speed(self):
        return self.speed

    def get_radius(self):
        return self.radius

    def get_id(self):
        return self.id

    def get_number_of_transition(self):
        return self.number_of_transitions

    def get_AB_index(self):
        return self.last_ABcollision_index

    def get_angle_speed(self):
        return atan2(self.speed[1], self.speed[0])

    def get_lattice_position(self):
        return self.lattice_pos

    def get_kinetic_energy(self):
        v_square = self.speed[0]**2 + self.speed[1]**2
        return self.mass*(v_square)/2.

    def get_cage_time(self, current_time):
        if self.number_of_transitions <= 0:
            return current_time
        return self.cage_time/self.number_of_transitions

    def get_energy_exchanged(self):
        return self.energy_exchange

    def get_fraction_energy_exch(self):
        return self.energy_exchange/self.get_kinetic_energy()

    def close_files(self):
        self.outfile_hist_jump_times.close()
        self.outfile_pos.close()

    def set_initial_relative_pos_coll(self, particle_B):
        dx = self.real_position[0] - particle_B.get_real_position()[0]
        dy = self.real_position[1] - particle_B.get_real_position()[1]

        self.relative_position_of_last_collision = [dx, dy]

    def set_speed(self, new_speed):
        for i in range(len(new_speed)):
            self.speed[i] = new_speed[i]

    def set_relative_last_collision(self, pos):
        self.relative_position_of_last_collision[0] = pos[0]
        self.relative_position_of_last_collision[1] = pos[1]

    def add_collision(self):
        self.number_of_collision += 1

    def get_number_of_collision(self):
        return self.number_of_collision

    def get_mass(self):
        return self.mass

    def change_free_path_0(self):
        dx = self.real_position[0] - self.position_of_last_collision[0]
        dy = self.real_position[1] - self.position_of_last_collision[1]
        self.free_path += (sqrt(dx**2 + dy**2))
        self.position_of_last_collision[0] = self.real_position[0]
        self.position_of_last_collision[1] = self.real_position[1]

    def change_free_path(self, particle_B):
        dx = self.real_position[0] - self.position_of_last_collision[0]
        dy = self.real_position[1] - self.position_of_last_collision[1]
        self.free_path += (sqrt(dx**2 + dy**2))
        self.position_of_last_collision[0] = self.real_position[0]
        self.position_of_last_collision[1] = self.real_position[1]

        dx = self.real_position[0] - particle_B.get_real_position()[0]
        dy = self.real_position[1] - particle_B.get_real_position()[1]
        dx = dx - self.relative_position_of_last_collision[0]
        dy = dy - self.relative_position_of_last_collision[1]

        self.relative_free_path += sqrt(dx * dx + dy * dy)

        dx = self.real_position[0] - particle_B.get_real_position()[0]
        dy = self.real_position[1] - particle_B.get_real_position()[1]
        self.relative_position_of_last_collision[0] = dx
        self.relative_position_of_last_collision[1] = dy
        particle_B.set_relative_last_collision([dx, dy])

    def get_diffusion_coefficient(self, particle_B, current_time):
        dx = self.real_position[0]  # - particle_B.get_real_position()[0]
        dy = self.real_position[1]  # - particle_B.get_real_position()[1]
        dx0 = self.initial_position[0]  # - particle_B.get_initial_position()[0]
        dy0 = self.initial_position[1]  # - particle_B.get_initial_position()[1]
        dx = dx - dx0
        dy = dy - dy0
        Dx = pow(dx, 2) / (2. * current_time)
        Dy = pow(dy, 2) / (2. * current_time)
        D = pow(dx, 2) + pow(dy, 2)
        D = D / (4. * current_time)

        return [Dx, Dy, D]

    def get_mean_free_path(self):
        if self.number_of_collision > 0:
            return self.free_path / self.number_of_collision
        else:
            delta = self.real_position[0]**2 + self.real_position[1]**2
            return (sqrt(delta))  # .real

    def get_relative_mean_free_path(self):
        if self.number_of_collision > 0:
            return self.relative_free_path / self.number_of_collision
        else:
            delta = 4 * (self.real_position[0]**2 + self.real_position[1]**2)
            return (sqrt(delta))  # .real

    def set_last_AB_index(self, _index):
        self.last_ABcollision_index = _index

    def change_speed_collision(self, _normal):
        norm = _normal[0]**2 +_normal[1]**2
        vx = self.speed[0] - 2*scalar_product(self.speed, _normal) * _normal[0] / norm
        vy = self.speed[1] - 2*scalar_product(self.speed, _normal) * _normal[1] / norm
        self.speed[0] = vx
        self.speed[1] = vy

    def Collide_with_particle(self, particle_B, vect):

        dx = self.position[0] + vect[0] - particle_B.get_position()[0]
        dy = self.position[1] + vect[1] - particle_B.get_position()[1]
        # dx = dx.real
        # dy = dy.real
        sigma = (self.radius + particle_B.get_radius())
        # sigma = sigma.real
        v_rel = [self.speed[0] - particle_B.speed[0], self.speed[1] - particle_B.speed[1]]

        momentum_involved_A = scalar_product((self.speed[0],self.speed[1]),(dx,dy))/sqrt(scalar_product((dx,dy),(dx,dy)))
        speed_B = particle_B.get_speed()
        momentum_involved_B = scalar_product(speed_B, (dx, dy)) / sqrt(scalar_product((dx, dy), (dx, dy)))

        self.energy_exchange = (momentum_involved_A**2)/2 + (momentum_involved_B**2)/2

        J = 2*self.mass*particle_B.get_mass()*(scalar_product(v_rel, (dx, dy))) / (sigma*(self.mass+particle_B.get_mass()))

        Jx = J * dx / sigma
        Jy = J * dy / sigma

        # change particle 1
        vx = self.speed[0] - Jx / self.mass
        vy = self.speed[1] - Jy / self.mass

        self.speed[0] = vx
        self.speed[1] = vy

        # change particle 2
        vx = particle_B.speed[0] + Jx / particle_B.get_mass()
        vy = particle_B.speed[1] + Jy / particle_B.get_mass()

        particle_B.set_speed([vx, vy])

        # change number of collisions
        self.number_of_collision += 1
        particle_B.add_collision()

        # change free path
        self.change_free_path(particle_B)
        particle_B.change_free_path_0()

    def get_last_collision(self):
        return self.last_collision

    def get_real_position(self):
        return self.real_position

    def set_last_collision(self, a):
        self.last_collision = a

    def free_flight(self, _time_c):
        x_old = self.position[0]
        y_old = self.position[1]
        self.position[0] = x_old + _time_c * self.speed[0]
        self.position[1] = y_old + _time_c * self.speed[1]

        x_old = self.real_position[0]
        y_old = self.real_position[1]
        self.real_position[0] = x_old + _time_c * self.speed[0]
        self.real_position[1] = y_old + _time_c * self.speed[1]

    def pass_trough_boundary(self, vector, current_time, lattice_deplacement):
        self.outfile_pos.write("%f %f %f \n" % (current_time, self.position[0], self.position[1]))
        self.position[0] = self.position[0] - (vector[0])
        self.position[1] = self.position[1] - (vector[1])
        self.last_ABcollision_index = -1

        if lattice_deplacement[1] == 0:
            if lattice_deplacement[0] == 1:
                self.lattice_pos[0] += 1
            elif lattice_deplacement[0] == -1:
                self.lattice_pos[0] += -1

        self.number_of_transitions += 1
        self.cage_time += current_time - self.time_of_last_entry
        self.outfile_hist_jump_times.write("%f \n" % (current_time - self.time_of_last_entry))
        self.time_of_last_entry = current_time
        self.outfile_pos.write("%f %f %f \n" % (current_time, self.position[0],self.position[1]))


# Boundaries of the unit cell
# no solid walls
class wall:
    def __init__(self, point_1, point_2, _id):
        # accept any object with the get_position method, for instance, discs
        dy = (point_1.get_position())[1] - (point_2.get_position())[1]
        dx = (point_1.get_position())[0] - (point_2.get_position())[0]

        _angle = atan2(dy, dx)

        m = tan(_angle)
        n = point_1.get_position()[1] - m*point_1.get_position()[0]
        l = sqrt(dx**2 + dy**2)

        self.position = [point_1.get_position()[0], point_1.get_position()[1]]
        self.angle = _angle
        self.length = l
        self.id = _id
        self.parameters = [m, n]
        self.periodic_wall_id = 0

    def get_position(self):
        return self.position

    def get_angle(self):
        return self.angle

    def get_length(self):
        return self.length

    def get_id(self):
        return self.id

    def get_position_at_distance(self, distance):
        if distance > self.length:
            return -1
        else:
            return (self.position[0] + distance*cos(self.angle), self.position[1] + distance*sin(self.angle))

    def get_parameters(self):
        return self.parameters

    def get_normal(self, point_1):
        if self.parameters[0] == 0:
            return (1., 0.)
        else:
            angle = atan(1. / self.parameters[0])
            return (cos(angle), sin(angle))

    def get_periodic_wall_id(self):
        return self.periodic_wall_id

    def set_periodic_wall_id(self, _id):
        self.periodic_wall_id = _id

    def get_lattice_deplacement(self):
        if self.id == -10:
            return (-1, 1)

        if self.id == -11:
            return (0, 1)

        if self.id == -12:
            return (1, 0)

        if self.id == -13:
            return (1, -1)

        if self.id == -14:
            return (0, -1)

        if self.id == -15:
            return (-1, 0)

# Solid disc inserted in the vertex of the unit cell
class obstacle:
    # position = the position of the center of the disc
    def __init__(self, _position, _radius, _id):
        self.position = _position
        self.radius = _radius
        self.id = _id

    def get_position(self):
        return self.position

    def get_radius(self):
        return self.radius

    def get_normal(self, point):
        dy = point[1] - self.position[1]
        dx = point[0] - self.position[0]

        alpha = atan2(dy, dx)
        return (cos(alpha), sin(alpha))

    def get_id(self):
        return self.id
