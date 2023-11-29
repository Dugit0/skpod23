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

struct Ind {
    int i;
    int j;
    int k;
};


struct IndArr {
    struct Ind *arr;
    int size;
    int cap;
};

struct IndArr *init_indarr(void) {
    struct IndArr *array = calloc(1, sizeof(*array));
    array->cap = 16;
    array->size = 0;
    array->arr = calloc(array->cap, sizeof(struct Ind));
    return array;
}

void push(struct IndArr *arr, struct Ind x) {
    if (arr->size == arr->cap) {
        arr->cap *= 2;
        arr->arr = realloc(arr->arr, arr->cap * sizeof(struct Ind));
    }
    (arr->arr)[arr->size] = x;
    (arr->size)++;
}

void free_arr(struct IndArr *arr) {
    free(arr->arr);
    free(arr);
}


struct IndArr **get_hyper(void) {
    int hyp_size = 3 * N - 2;
    struct IndArr **res = calloc(hyp_size, sizeof(struct IndArr *));
    for (i = 0; i < hyp_size; i++) {
        res[i] = init_indarr();
    }
#pragma omp parallel for collapse(3) default(none) private(i, j, k) shared(res)
    for (i = 1; i <= N - 2; i++) {
        for (j = 1; j <= N - 2; j++) {
	        for (k = 1; k <= N - 2; k++) {
                int ind = i + j + k;
                struct Ind trind = {i, j, k};
#pragma omp critical
                push(res[ind], trind);
            }
        }
    }
    return res;
}




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


void relax(struct IndArr **hyper)
{
    for (int t = 3; t < 3 * N - 5; t++) {
        /* printf("t = %d\n", t); */
        for (int q = 0; q < hyper[t]->size; q++) {
            i = (hyper[t]->arr)[q].i;
            j = (hyper[t]->arr)[q].j;
            k = (hyper[t]->arr)[q].k;
            /* printf("%d %d %d\n", i, j, k); */
            double e;
            e = A[i][j][k];
            A[i][j][k] = (A[i-1][j][k] + A[i+1][j][k] + A[i][j-1][k] + A[i][j+1][k] + A[i][j][k-1] + A[i][j][k+1]) / 6.;
            eps=Max(eps, fabs(e - A[i][j][k]));
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
	printf("  S = %f\n",s);

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
    struct IndArr **hyper = get_hyper();

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
		relax(hyper);
		if (eps < maxeps) {
            break;
        }
	}
	verify();
    double end = omp_get_wtime();
    printf("Time = %f\n", end - start);
    

    for (i = 0; i < 3 * N - 2; i++) {
        free_arr(hyper[i]);
    }
    free(hyper);
	return 0;
}


