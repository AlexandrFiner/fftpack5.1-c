#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"

extern int cfft1i_(integer *, real *, integer *, integer *);
extern int cfft1b_(integer *, integer *, complex *, integer *
					lenc, real *, integer *, real *, integer *,
					integer *);
extern int cfft1f_(integer *, integer *, complex *, integer *,
					real *, integer *, real *, integer *,
					integer *);

// Переименование структуры complex в вашем коде

int main(int argc, char* argv[]) {
    integer I, N, LENSAV, IER, LENWRK;
    N = 1000;
    LENSAV = 2013;
    LENWRK = 2 * N;

    float RR, RI;

    static real WSAVE[2013];
	static real WORK[2000];

    printf("Program cfft1i and related messages:\n");
    cfft1i_(&N, WSAVE, &LENSAV, &IER);


    // double realPart = 2.0;
    // double imagPart = 3.0;
	// complex z = {20.0, 30.0};

    RR = 20.0;
    RI = 10.0;

	static integer INC = 1;
	complex C[1000];
	for(int i = 0; i < 1000; i++) {
		C[i].r = RR;
		C[i].i = RI;
	}
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);

	for(int i = 0; i < 2013; i++) {
		printf("%f\n", WSAVE[i]);
	}

	// for(int i = 0; i < N; i++) {
	//	printf("%f", C[i].r);
	//	printf("%f", C[i].i);
	// }


    return 0;
}
