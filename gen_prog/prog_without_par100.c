#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N  100
double maxeps = 0.1e-7;
int itmax = 100;
int i,j,k;
double eps;

double A [N][N][N];


void init() { 
    for (i = 0; i <= N - 1; i++) {
        for (j = 0; j <= N - 1; j++) {
	        for (k = 0; k <= N - 1; k++) {
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1) {
                    A[i][j][k] = 0.;
                } else {
                    A[i][j][k] = (4. + i + j + k);
                }
            }
        }
    }
}


void relax()
{

    for (i = 1; i <= N - 2; i++) { 
        for (j = 1; j <= N - 2; j++) {
	        for (k = 1; k <= N - 2; k++) {
                double e;
                e = A[i][j][k];
                A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i][j-1][k] + A[i][j+1][k] + A[i][j][k-1] + A[i][j][k+1]) / 6.;
                eps=Max(eps, fabs(e - A[i][j][k]));
            }    
        }
    }
}


void verify()
{ 
	double s;

	s = 0.;
	for (i = 0; i <= N - 1; i++) {
	    for (j = 0; j <= N - 1; j++) {
	        for (k = 0; k <= N - 1; k++) {
                s = s + A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
            }
        }
    }
	/* printf("  S = %f\n",s); */

}


void printa() {
    for (i = 0; i <= N - 1; i++) {
        for (j = 0; j <= N - 1; j++) {
	        for (k = 0; k <= N - 1; k++) {
                printf("%f ", A[i][j][k]);
            }
        }
    }
    return;
}


int main(int an, char **as) {
    double start = omp_get_wtime();
    int it;
	init();
	for (it=1; it<=itmax; it++) {
		eps = 0.;
		relax();
		if (eps < maxeps) {
            break;
        }
	}
	verify();
    double end = omp_get_wtime();
    printf("Time = %f\n", end - start);
	return 0;
}


