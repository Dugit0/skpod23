#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define  Max(a,b) ((a)>(b)?(a):(b))

#define  N  200 
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
    for (int H = 3; H < 3 * N - 5; H++) {
#pragma omp parallel for collapse(1) default(none) private(i, j, k) shared(A, H) reduction(max:eps)
        for (int v = 0; v <= H; v++) {
            /* for (int u = 0; u <= H; u++) { */
            /*     if (v < u) { */
            /*         continue; */
            /*     } */
            for (int u = 0; u <= v; u++) {
                i = u;
                j = v - u;
                k = H - v;
                /* printf("%d %d %d\n", i, j, k); */
                if (i == 0 || j == 0 || k == 0 || i >= N - 1 || j >= N - 1 || k >= N - 1) {
                    continue;
                }
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
#pragma omp parallel for collapse(3) default(none) private(i, j, k) shared(A) reduction(+:s)
	for (i = 0; i <= N - 1; i++) {
	    for (j = 0; j <= N - 1; j++) {
	        for (k = 0; k <= N - 1; k++) {
                s = s + A[i][j][k]*(i+1)*(j+1)*(k+1)/(N*N*N);
            }
        }
    }
	/* printf("  S = %f\n",s); */

}


int main(int an, char **as) {
	

    double start = omp_get_wtime();
	int it;
	init();
    /* int hyp_size = 3 * N - 2; */
    /* for (int x = 0; x < hyp_size; x++) { */
    /*     printf("%d : ", x); */
    /*     for (int y = 0; y < hyper[x]->size; y++) { */
    /*         printf("(%d %d %d) ", (hyper[x]->arr)[y].i, (hyper[x]->arr)[y].j, (hyper[x]->arr)[y].k); */
    /*     } */
    /*     printf("\n"); */
    /* } */
	
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


