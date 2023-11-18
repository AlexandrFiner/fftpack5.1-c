#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

#define lensav(n) (2*n + (int)(log2(n)) + 6)
//
//int main(int argc, char* argv[]) {
//
//    int n = 1;
//    float wsave = 10;
//    int ier;
//
//    cfftmi_(n, wsave, lensav(n), ier);
//    printf("Hello world %d\n", ier);
//    return EXIT_SUCCESS;
//}

int cfftmi_(int *n, float *wsave, int *lensav, int *ier);

int main(int argc, char* argv[]) {

    int n = 2;
    float wsave = 0.1;
    int ier = 5;
    int lensv = lensav(n);

    cfftmi_(&n, &wsave, &lensv, &ier);
    printf("Hello world %d\n", ier);
    return EXIT_SUCCESS;
}
