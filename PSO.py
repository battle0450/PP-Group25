from mpi4py import MPI
import numpy as np
import random
import time
import math                    
import random                                       

n_iterations = 100
n_particles = 10
c1 = 0.1
c2 = 0.1
W = 0.2

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
    space = Space(n_particles)
    particles_swarm = [Particle() for _ in range(space.n_particles)]
    space.particles = particles_swarm
    start = time.time()
    while(iteration < n_iterations):
        space.update_pbest()
        space.update_gbest()
        space.move_particles()
        iteration +=1
    print("Executing time: ", time.time() - start)
    
