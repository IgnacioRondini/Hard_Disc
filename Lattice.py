
from random import randint, random

# Implement a 2D Lattice with a hexagonal unit cell

class Lattice:

    def __init__(self, _size, _prob_boundaries, _prob_reaction):
        # size = (N,M)
        self.size = _size

        # A[i] = [(i,j),id]
        self.color_A = []
        self.color_B = []
        self.reaction_histogram = [0] * (_size[0]*2 + 1)
        self.prob_boundaries = _prob_boundaries  # (p,q)
        self.prob_reaction = _prob_reaction
        self.number_of_reactions = 0
        self.time_of_the_last_reaction = 0
        self.just_react = False

    def get_size(self):
        return self.size

    def get_area(self):
        # to implement
        return -1

    def get_A_particles(self):
        return self.color_A

    def get_B_particles(self):
        return self.color_B

    def get_number_of_reactios(self):
        return self.number_of_reactions

    def get_number_of_colored_particles(self):
        return [len(self.color_A), len(self.color_B)]

    def get_time_last_reaction(self):
        return self.time_of_the_last_reaction

    def get_reaction_histogram(self):
        return self.reaction_histogram

    def has_just_react(self):
        return self.just_react

    def reset_just_react(self):
        self.just_react = False

    def set_last_time_reaction(self, current_time):
        self.time_of_the_last_reaction = current_time

    def add_particle(self, color, position, _id):
        # A color = 1
        # B color = -1
        if color == 1:
            self.color_A.append([position, _id])
        if color == -1:
            self.color_B.append([position, _id])

    def set_initial_conditions(self, Na, Nb):
        N = self.size[0]
        M = self.size[1]
        if Na + Nb > 2 * (2*N + 1) * (2*M + 1):
            print "error, too many colored particles"
            exit()

        count_a = 0
        count_b = 0
        while count_a < Na:
            if random() <= 0.5:
                _id = 100
            else:
                _id = 200
            i = randint(-N, -1)
            j = randint(-M, M)

            if [[i, j], _id] not in self.color_A:
                self.color_A.append([[i, j], _id])
                count_a += 1

        while count_b < Nb:
            if random() <= 0.5:
                _id = 100
            else:
                _id = 200
            i = randint(1, N)
            j = randint(-M, M)

            if [[i, j], _id] not in self.color_A:
                if [[i, j], _id] not in self.color_B:
                    self.color_B.append([[i, j], _id])
                    count_b += 1

    def update_move(self, _id, lattice_vector):

        N = self.size[0]
        M = self.size[1]

        for i in xrange(len(self.color_A) - 1, -1, -1):
            if self.color_A[i][1] == _id:

                if lattice_vector[1] == 0:
                    if lattice_vector[0] == 1:
                        if self.color_A[i][0][0] == self.size[0]:

                            del self.color_A[i]
                        else:
                            self.color_A[i][0][0] = self.color_A[i][0][0] +1

                    elif lattice_vector[0] == -1:
                        if self.color_A[i][0][0] == -self.size[0]:
                            del self.color_A[i]
                        else:
                            self.color_A[i][0][0] = self.color_A[i][0][0] -1

        for i in xrange(len(self.color_B) - 1, -1, -1):

            if self.color_B[i][1] == _id:
                if lattice_vector[1] == 0:
                    if lattice_vector[0] == 1:
                        if self.color_B[i][0][0] == self.size[0]:
                            del self.color_B[i]

                        else:
                            self.color_B[i][0][0] = self.color_B[i][0][0] + 1

                    elif lattice_vector[0] == -1:
                        if self.color_B[i][0][0] == -self.size[0]:
                            del self.color_B[i]

                        else:
                            self.color_B[i][0][0] = self.color_B[i][0][0] - 1

        # add the particles if it was a deplacement with a probability p
        if lattice_vector[0] == 1 and lattice_vector[1] == 0:
            if random() <= self.prob_boundaries[0]:
                # add a left particle (unidimensional j = 0)
                self.add_particle(1, [-self.size[0], 0], _id)

        # add the particles if it was a deplacement
        if lattice_vector[0] == -1 and lattice_vector[1] == 0:
            if random() <= self.prob_boundaries[1]:

                # add a right particle
                self.add_particle(-1, [self.size[0], 0], _id)
#

    def update_collision(self, lattice_vector, energy_exch, current_time):

        N = self.size[0]
        M = self.size[1]

        for a in self.color_A:
            for b in self.color_B:
                # check if they are different particles
                # in the unit cell

                if a[1] != b[1]:
                    new_index_a = [0, 0]
                    new_index_a[0] = a[0][0]
                    new_index_a[1] = a[0][1]
                    new_index_b = [0, 0]
                    new_index_b[0] = b[0][0]
                    new_index_b[1] = b[0][1]

                    if lattice_vector[1] == 0:

                        new_index_a[0] += lattice_vector[0]
                        new_index_a[1] += lattice_vector[1]

                        new_index_b[0] += lattice_vector[0]
                        new_index_b[1] += lattice_vector[1]

                    if new_index_a == b[0] or new_index_b == a[0]:
                        # start collide colors ------------------------------

                        condition = False
                        if energy_exch < 0:
                            print "ERROR energy"
                            exit()
                        if random() >= self.prob_reaction:

                            condition = True

                        if condition is True:
                            # erase a,b
                            self.number_of_reactions += 1
                            self.just_react = True
                            index = self.size[0] + a[0][0]
                            self.reaction_histogram[index] += 1
                            self.color_A.remove(a)
                            self.color_B.remove(b)

                            break

                        if condition is True:
                            print "condition error"
                            exit()
                        # end collide colors --------------------------------
#

    def get_density_distribution(self, _string):
        outfile_A = open(_string + "A.dat", "w")
        outfile_B = open(_string + "B.dat", "w")

        for a in self.color_A:
            outfile_A.write(" %f %f %f \n" % (a[0][0], a[0][1], a[1]))
        for b in self.color_B:
            outfile_B.write(" %f %f %f \n" % (b[0][0], b[0][1], b[1]))

        outfile_A.close()
        outfile_B.close()

    def get_min_max(self):
        # return the rigthmost A particle and the leftmost B particle' position
        max_pos = -self.size[0]
        min_pos = self.size[0]

        for a in self.color_A:
            if max_pos < a[0][0]:
                max_pos = a[0][0]
        for b in self.color_B:
            if min_pos > b[0][0]:
                min_pos = b[0][0]
        return (max_pos, min_pos)

    def print_histo_pos(self, str):
        outfile = open(str, "w")

        for i in range(len(self.reaction_histogram)):
            outfile.write("%f %f \n" % (-self.size[0] + i, self.reaction_histogram[i]))

        outfile.close()
