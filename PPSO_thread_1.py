#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 00:06:58 2020
@author: jifbvhqp
"""

from mpi4py import MPI
import numpy as np
import random                                       
import threading
n_iterations = 100
n_particles = 10000
c1 = 0.1
c2 = 0.1
W = 0.2
numThread = 4

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

    def update_pbest(self,start,end):
        '''
        for particle in self.particles:
            fitness = self.fitness(particle)
            if(fitness < particle.pbest_value):
                particle.pbest_value = fitness
                particle.pbest_position = particle.position
        '''
        for i in range(start,end+1):
            fitness = self.fitness(self.particles[i])
            if fitness < self.particles[i].pbest_value:
                self.particles[i].pbest_value = fitness
                self.particles[i].pbest_position = self.particles[i].position
    
    def update_gbest(self):
        for particle in self.particles:
            fitness = self.fitness(particle)
            if(fitness < self.gbest_value):
                self.gbest_value = fitness
                self.gbest_position = particle.position
                #print(self.gbest_value)

    def move_particles(self,start,end):
        '''
        for particle in self.particles:
            global W
            new_velocity = W * particle.velocity + \
                            c1 * random.random() * (particle.pbest_position - particle.position) + \
                            c2 * random.random() * (self.gbest_position - particle.position)
            particle.velocity = new_velocity
            particle.move()
        '''
        for i in range(start,end+1):
            global W
            self.particles[i].velocity = W * self.particles[i].velocity + \
                            c1 * random.random() * (self.particles[i].pbest_position - self.particles[i].position) + \
                            c2 * random.random() * (self.gbest_position - self.particles[i].position)
            self.particles[i].move()
        
        
def job(space,startParticleIndex,endParticleIndex):
   # print(startParticleIndex," ",endParticleIndex)
   space.update_pbest(startParticleIndex,endParticleIndex)
   
    
def job2(space,startParticleIndex,endParticleIndex):
    space.move_particles(startParticleIndex,endParticleIndex)
 
if __name__ == '__main__':
    start = MPI.Wtime()
    iteration = 0
    space = Space(n_particles)
    space.particles = [Particle() for _ in range(n_particles)]
    onePerThread = n_particles//numThread
    particles_frame = []

        
    while(iteration < n_iterations):
        index = 0
        threads = []
        threads_2 = []
        for i in range(numThread):
            threads.append(threading.Thread(target = job, args = (space,index,index+(onePerThread-1),))) 
            threads[i].start()
            index += onePerThread

        for i in range(numThread):
            threads[i].join()
        index = 0
        
        space.update_gbest()
        
        for i in range(numThread):
            threads_2.append(threading.Thread(target = job2, args = (space,index,index+(onePerThread-1),)))
            threads_2[i].start()
            index += onePerThread
        for i in range(numThread):
            threads_2[i].join()
            
        iteration +=1
        
        
        
    #while(iteration < n_iterations):
    #    space.update_pbest()
    #    space.update_gbest()
    #    space.move_particles()
    #    iteration +=1
    print("Executing time: ", MPI.Wtime() - start)  
