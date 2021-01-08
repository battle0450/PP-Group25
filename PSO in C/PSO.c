#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>
#include <math.h>
#define PI 3.14159265358979323846
int nDimensions, nIterations;
double x_min = -10.f;
double x_max = 50.f;
double e = 2.71828182845904523536;

double fitness(double x[], int dim)
{
     
    double c = 2 * PI;
    double b = 0.2;
    double a = 20;
    double sum1 = 0;
    double sum2 = 0;
    int i;
    for (i=0; i < dim; i++) {
        sum1 = sum1 + pow(x[i], 2);
        sum2 = sum2 + cos(c * x[i]);
    }
    double term1 = -a * exp(-b * sqrt(sum1 / dim));
    double term2 = -exp(sum2 / dim);
    return term1 + term2 + a + e;
}

int main(int argc, char *argv[]) {
    int i,j;
    double nParticles;
    //Argument handling START
    nDimensions = 5;
    nParticles = 100;
    nIterations = 500;
    int seed = 1;

    int size,myrank,distributed_particles=nParticles;
    double result[(int)distributed_particles];
    int step;
    double a,b;
    double c1, c2, rho1, rho2, w, fit;
    c1 = c2 = 0.1f;
    w = 0.2f;
    
    //Random number generator initialization

    double positions[(int)distributed_particles][(int)nDimensions];
    double velocities[(int)distributed_particles][(int)nDimensions];
    double pBestPositions[(int)distributed_particles][(int)nDimensions];    
    double pBestFitness[(int)distributed_particles];
    double gBestPosition[(int)nDimensions];
    double gBestFitness = RAND_MAX;
    int min;
    double start = clock();
    //particle initialization
    //#pragma omp parallel for private(a,b)  reduction(min:gBestFitness)
    for (i=0; i<distributed_particles; i++) {        
        for (j=0; j<(int)nDimensions; j++) 
        {
            a = x_min + (x_max - x_min) *  (double)rand_r(&seed)/RAND_MAX;
            b = x_min + (x_max - x_min) *  (double)rand_r(&seed)/RAND_MAX;
            positions[i][j] = a;
            pBestPositions[i][j] = a;
            velocities[i][j] = (a-b) / 2.;
        }
        pBestFitness[i] = fitness(positions[i], nDimensions);
        if (pBestFitness[i] < gBestFitness) {
            gBestFitness = pBestFitness[i];
            memmove((void *)gBestPosition, (void *)&positions[i], sizeof(double) * nDimensions);
            
        } 
    }

    //actual calculation
    for (step=0; step<nIterations; step++) {
        //#pragma omp parallel num_threads(4) shared(min)
        {

        //#pragma omp for private(a,b) 
        for (i=0; i<distributed_particles; i++) {
             
            for (j=0; j<nDimensions; j++) {
                // calculate stochastic coefficients
                rho1 = c1 * (double)rand_r(&seed)/RAND_MAX;
                rho2 = c2 * (double)rand_r(&seed)/RAND_MAX;
                // update velocity
                velocities[i][j] = w * velocities[i][j] + \
                rho1 * (pBestPositions[i][j] - positions[i][j]) +  \
                rho2 * (gBestPosition[j] - positions[i][j]);
                // update position
                positions[i][j] += velocities[i][j];

                if (positions[i][j] < x_min) {
                    positions[i][j] = x_min;
                    velocities[i][j] = 0;
                } else if (positions[i][j] > x_max) {
                    positions[i][j] = x_max;
                    velocities[i][j] = 0;
                }

            }

            // update particle fitness
            fit = fitness(positions[i], nDimensions);
            // update personal best position?
            if (fit < pBestFitness[i]) {
                pBestFitness[i] = fit;
                // copy contents of positions[i] to pos_b[i]
                memmove((void *)&pBestPositions[i], (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
            // update gbest??
           
            
                
        }

            //#pragma omp for reduction(min:gBestFitness)
            for(i=0;i<(int)nParticles;i++)
            if (pBestFitness[i] < gBestFitness) {
                        // update best fitness
                        gBestFitness = pBestFitness[i];
                        // copy particle pos to gbest vector
                        }
            
            //#pragma omp  for  
             for(i=0;i<(int)nParticles;i++)
             {if (gBestFitness==pBestFitness[i])
                min=i;  
               
            }
        }
         memmove((void *)gBestPosition, (void *)&pBestPositions[min],sizeof(double) * nDimensions);
            
    }
    double finish = clock();
    double overhead = (double)(finish - start)/CLOCKS_PER_SEC;
    printf("Result: %f\n", gBestFitness);
    printf("\n Time in Seconds (T) : %lf\n", overhead);
    
    

}
