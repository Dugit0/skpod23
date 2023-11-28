#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N  500
double maxeps = 0.1e-7;
int itmax = 100;
int i,j,k;
double eps;

double A [N][N][N];

void init() { 
#pragma omp parallel for collapse(3) default(none) private(i, j, k) shared(A)
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
#pragma omp parallel for ordered(3) default(none) private(i, j, k) shared(A) reduction(max:eps)
    for (i = 1; i <= N - 2; i++) { 
        for (j = 1; j <= N - 2; j++) {
	        for (k = 1; k <= N - 2; k++) {
                double e;
                e = A[i][j][k];
#pragma omp ordered depend(sink: i - 1, j, k) depend(sink: i, j - 1, k) depend(sink: i, j, k - 1)
                A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i][j-1][k] + A[i][j+1][k] + A[i][j][k-1] + A[i][j][k+1]) / 6.;
#pragma omp ordered depend(source)
                eps=Max(eps, fabs(e - A[i][j][k]));
            }    
        }
    }
}


void verify()
{ 
	double s;

	s = 0.;
#pragma omp parallel for collapse(3) default(none) private(i, j, k) shared(A) reduction(+:s)
	for (i = 0; i <= N - 1; i++) {
	    for (j = 0; j <= N - 1; j++) {
	        for (k = 0; k <= N - 1; k++) {
                s += A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
            }
        }
    }
	/* printf("S = %f\n",s); */
}


int main(int an, char **as) {
    double start = omp_get_wtime();
    int it;
	init();
	for (it=1; it<=itmax; it++) {
		eps = 0.;
		relax();
		/* printf( "it=%4i   eps=%f\n", it,eps); */
        
		if (eps < maxeps) {
            break;
        }
	}
	verify();
    double end = omp_get_wtime();
    printf("Time = %f\n", end - start);
	return 0;
}


