#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

extern int cfft1i_(integer *, real *, integer *, integer *);
extern int cfft1b_(integer *, integer *, complex *, integer *
					lenc, real *, integer *, real *, integer *,
					integer *);

// Переименование структуры complex в вашем коде

int main(int argc, char* argv[]) {
    integer I, N, LENSAV, IER, LENWRK;
    N = 1000;
    LENSAV = 2013;
    LENWRK = 2 * N;

    float RR, RI;

    static real WSAVE[2013];
	static real WORK[2013];

    printf("Program cfft1i and related messages:\n");
    cfft1i_(&N, WSAVE, &LENSAV, &IER);

    RR = 20.0;
    RI = 10.0;

    complex C = {RR, RI};

	int ONE = 1;
	// cfft1b_(&N, &ONE, &C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
	// for(int i = 0; i < 2013; i++)
	//	printf("%f", WORK[i]);

    return 0;
}
