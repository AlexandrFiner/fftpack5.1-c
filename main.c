#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

#define lensav(n) (2*n + (int)(log2(n)) + 6)

int cfftmi_(int *n, float *wsave, int *lensav, int *ier);

int main(int argc, char* argv[]) {

    int n = 2;          // Длина сигнала
    float wsave = 0.1;  // Частота
    int ier = 5;
    int lensv = lensav(n);

    cfftmi_(&n, &wsave, &lensv, &ier);
    printf("Hello world %d\n", ier);
    return EXIT_SUCCESS;
}
