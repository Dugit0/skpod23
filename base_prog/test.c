#include <stdio.h>

int main(void) {

    int N;
    int num_threads;
    scanf("%d %d", &N, &num_threads);
    int m = 0;
    int M = 0;
    for (int cur_m = 0; cur_m <= N; cur_m++) {
        int cur_M = (N + cur_m - 1) / cur_m;
        int next_M = (N + cur_m) / (cur_m + 1);
        if (next_M*next_M*next_M > num_threads && cur_M*cur_M*cur_M <= num_threads) {
            m = cur_m;
            M = next_M;
        }
    }
    printf("N = %d, num_threads = %d\nm = %d, M = %d\n", N, num_threads, m, M);
    return 0;
}
