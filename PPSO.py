from mpi4py import MPI
import numpy as np
import random
import time                                   

n_iterations = 100
n = 10000
c1 = 0.1
c2 = 0.1
W = 0.2

# Get our MPI communicator, our rank, and the world size.
comm = MPI.COMM_WORLD
world_rank = comm.Get_rank()
world_size = comm.Get_size()

class Particle():
    def __init__(self):
        self.position = np.array([(-1) ** (bool(random.getrandbits(1))) * \
                        random.random() * 50, (-1) ** (bool(random.getrandbits(1))) * \
                        random.random() * 50])
        self.pbest_position = self.position
        self.pbest_value = float('inf')
        self.velocity = np.array([0.0])
    
    def move(self):
        self.position = self.position + self.velocity

class Space():
    def __init__(self, n_particles):
        self.n_particles = n_particles
        self.particles = []
        self.gbest_position = np.array([random.random() * 50, random.random() * 50])
        self.gbest_value = float('inf')

    def fitness(self, particle):
        return (particle.position[0] - 20) ** 2 + (particle.position[1] - 20) ** 2 + 1

    def update_pbest(self):
        for particle in self.particles:
            fitness = self.fitness(particle)
            if(fitness < particle.pbest_value):
                particle.pbest_value = fitness
                particle.pbest_position = particle.position

    def update_gbest(self):
        for particle in self.particles:
            fitness = self.fitness(particle)
            if(fitness < self.gbest_value):
                self.gbest_value = fitness
                self.gbest_position = particle.position
                #print(self.gbest_value)    

    def move_particles(self):
        for particle in self.particles:
            global W
            new_velocity = W * particle.velocity + \
                            c1 * random.random() * (particle.pbest_position - particle.position) + \
                            c2 * random.random() * (self.gbest_position - particle.position)
            particle.velocity = new_velocity
            particle.move()
    
if __name__ == '__main__':
    iteration = 0
    global_space = Space(n)
    avg = int(n / world_size)
    space = Space(avg)
    particles_swarm = [Particle() for _ in range(avg)]
    space.particles = particles_swarm
    start = MPI.Wtime()

    while(iteration < n_iterations):
        space.update_pbest()
        iteration+=1
    if(world_rank != 0):
        comm.send(space.particles, dest = 0, tag = 0)
    iteration = 0

    if(world_rank == 0):
        for src in range(1, world_size):
            local = comm.recv(source = src, tag = 0)
            global_space.particles = global_space.particles + local
        while(iteration < n_iterations):
           global_space.update_gbest()
            global_space.move_particles()
            iteration +=1
        print("Executing time: ", MPI.Wtime() - start)
