//run : rm serial; gcc serialpso.c -lm -o serial; ./serial 100000 3 100
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>

// Time parameters
struct timeval TimeValue_Start;
struct timezone TimeZone_Start;
struct timeval TimeValue_Final;
struct timezone TimeZone_Final;
long time_start, time_end;
double time_overhead;

// Benchmark Function
double ackley(double x[], double nDimensions) {
    double c = 2*M_PI;
    double b = 0.2;
    double a = 20;
    double sum1 = 0;
    double sum2 = 0;
    int i;
    for (i=0; i<nDimensions; i++) {
        sum1 = sum1 + pow(x[i],2);
        sum2 = sum2 + cos(c*x[i]);
    }
    double term1 = -a * exp(-b*sqrt(sum1/nDimensions));
    double term2 = -exp(sum2/nDimensions);
    return term1 + term2 + a + M_E;
}

int main(int argc, char *argv[]) {
    
    // PSO parameters
    int nParticles, nDimensions, nIterations;
    double x_min = -50;
    double x_max = 50;
    double a,b;
    double c1, c2, rho1, rho2, w, fit;
    c1 = c2 = 1.496;
    w = 0.7298;
    
    nParticles = atoi(argv[1]);
    nDimensions = atoi(argv[2]);
    nIterations = atoi(argv[3]);

    gettimeofday(&TimeValue_Start, &TimeZone_Start);  
    
    // random number generator initialization
    unsigned int seed = time(NULL);

    // PSO space initailization
    double positions[(int)nParticles][(int)nDimensions];
    double velocities[(int)nParticles][(int)nDimensions];
    double pBestPositions[(int)nParticles][(int)nDimensions];    
    double pBestFitness[(int)nParticles];
    double gBestPosition[(int)nDimensions];
    double gBestFitness = DBL_MAX;

    // particle initialization
    for (int i=0; i<nParticles; i++) {
        for (int j=0; j<nDimensions; j++) {
            a = x_min + (x_max - x_min) *  ((double)rand_r(&seed)/RAND_MAX);
            b = x_min + (x_max - x_min) *  ((double)rand_r(&seed)/RAND_MAX);
            positions[i][j] = a;
            pBestPositions[i][j] = a;
            velocities[i][j] = (a-b) / 2.;
        }
        pBestFitness[i] = ackley(positions[i],nDimensions);
        if (pBestFitness[i] < gBestFitness) {
            memmove((void *)gBestPosition, (void *)&positions[i], sizeof(double) * nDimensions);
            gBestFitness = pBestFitness[i];
        } 
    }

    //actual calculation
    for (int step=0; step<nIterations; step++) {
        for (int i=0; i<nParticles; i++) {
            for (int j=0; j<nDimensions; j++) {
                // calculate stochastic coefficients
                rho1 = c1 * ((double)rand_r(&seed)/RAND_MAX);
                rho2 = c2 * ((double)rand_r(&seed)/RAND_MAX);
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
            fit = ackley(positions[i], nDimensions);
            // update personal best position?
            if (fit < pBestFitness[i]) {
                pBestFitness[i] = fit;
                // copy contents of positions[i] to pos_b[i]
                memmove((void *)&pBestPositions[i], (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
            // update gbest??
            if (fit < gBestFitness) {
                // update best fitness
                gBestFitness = fit;
                // copy particle pos to gbest vector
                memmove((void *)gBestPosition, (void *)&positions[i],
                    sizeof(double) * nDimensions);
            }
        }
    }
    
    // print pso parameters
    printf("nParticles  : %d\n", nParticles);
    printf("nDimensions : %d\n", nDimensions);
    printf("nIterations : %d\n", nIterations);
    printf("Best Fitness : %f\n", gBestFitness);
    printf("Best Position: ");
    for(int i=0; i < nDimensions; i++)  printf("%lf ", gBestPosition[i]);
    gettimeofday(&TimeValue_Final, &TimeZone_Final);
    time_start = TimeValue_Start.tv_sec * 1000000 + TimeValue_Start.tv_usec;
    time_end = TimeValue_Final.tv_sec * 1000000 + TimeValue_Final.tv_usec;
    time_overhead = (time_end - time_start)/1000000.0;
    printf("\nTime in Seconds (T) : %lf\n",time_overhead);
}
