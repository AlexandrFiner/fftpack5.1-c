#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

//#define lensav(n) (2*n + (int)(log2(n)) + 6)
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

int main(int argc, char* argv[]) {

    int n = 1;
    float wsave = 10;
    int ier;

//    cfftmi_(n, wsave, lensav(n), ier);
    printf("Hello world %d\n", ier);
    return EXIT_SUCCESS;
}