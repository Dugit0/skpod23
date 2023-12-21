#include <stdlib.h>
#include <stdio.h>

int main(void) {
    int ***arr = calloc(1, sizeof(*arr));
    arr[0] = calloc(1, sizeof(*(arr[0])));
    arr[0][0] = calloc(1, sizeof(*(arr[0][0])));
    arr[0][0][0] = 123;
    printf("%p\n%p\n%p\n%d", arr, arr[0], arr[0][0], arr[0][0][0]);
    return 0;
}
