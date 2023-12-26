#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include <time.h>

extern int cfft1i_(integer *, real *, integer *, integer *);
extern int cfft1b_(integer *, integer *, complex *, integer *
					lenc, real *, integer *, real *, integer *,
					integer *);
extern int cfft1f_(integer *, integer *, complex *, integer *,
					real *, integer *, real *, integer *,
					integer *);

// Переименование структуры complex в вашем коде

int main(int argc, char* argv[]) {
	FILE *file = fopen("result.txt", "w");

    // Проверить, был ли файл успешно открыт
    if (file == NULL) {
        fprintf(stderr, "Ошибка открытия файла.\n");
        return 1;  // Вернуть код ошибки
    }

    integer I, N, LENSAV, IER, LENWRK;
    N = 1000;
    LENSAV = 2013;
    LENWRK = 2 * N;

    real RR, RI;

    static real WSAVE[2013];
	static real WORK[2000];

	srand(time(NULL));

	static integer INC = 1;
	complex C[1000];

	complex NEED = { RR, RI };
	real diff = 0.0;

	// Back-forward
	/**
	--- IDENTIFY TEST AND INITIALIZE FFT
	*/
	printf("Program cfft1i and related messages:\n");
	cfft1i_(&N, WSAVE, &LENSAV, &IER);

    fprintf(file, "Back-forward\n");
    fprintf(file, "IDENTIFY TEST AND INITIALIZE FFT\n");
    fprintf(file, "WSAVE = [");
	for(int i = 0; i < LENSAV; i++) {
		fprintf(file, "%f, ", WSAVE[i]);
	}
	fprintf(file, "]\n");

	/**
	--- GENERATE TEST VECTOR FOR BACKWARD-FORWARD TEST
	*/
    RR = (float)rand() / RAND_MAX;
    RI = (float)rand() / RAND_MAX;
	NEED.r = RR;
	NEED.i = RI;
	fprintf(file, "RR - %f, RI - %f\n", RR, RI);
	for(int i = 0; i < N; i++) { C[i].r = RR; C[i].i = RI; }

	/**
	--- PERFORM BACKWARD TRANSFORM
	*/
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "BACKWARD WORK = [");
	for(int i = 0; i < LENWRK; i++) {
		fprintf(file, "%f, ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PERFORM FORWARD TRANSFORM
	*/
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "FORWARD WORK = [");
	for(int i = 0; i < LENWRK; i++) {
		fprintf(file, "%f, ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PRINT TEST RESULTS
	*/
	diff = 0.0;
	for(int i = 0; i < N; i++) {
		diff = max(diff, abs(C[i].r - NEED.r) + abs(C[i].i - NEED.i));
	}
	printf("CFFT1 BACKWARD-FORWARD MAX ERROR = %f\n", diff);

	fprintf(file, "\n");

	// Forward-back
	/**
	--- IDENTIFY TEST AND INITIALIZE FFT
	*/
	printf("Program cfft1i and related messages:\n");
	cfft1i_(&N, WSAVE, &LENSAV, &IER);

    fprintf(file, "Forward-back\n");
    fprintf(file, "IDENTIFY TEST AND INITIALIZE FFT\n");
    fprintf(file, "WSAVE = [");
	for(int i = 0; i < LENSAV; i++) {
		fprintf(file, "%f, ", WSAVE[i]);
	}
	fprintf(file, "]\n");

	/**
	GENERATE TEST VECTOR FOR FORWARD-BACKWARD TEST
	*/
    RR = (float)rand() / RAND_MAX;
    RI = (float)rand() / RAND_MAX;
	NEED.r = RR;
	NEED.i = RI;
	fprintf(file, "RR - %f, RI - %f\n", RR, RI);
	for(int i = 0; i < N; i++) { C[i].r = RR; C[i].i = RI; }

	/**
	--- PERFORM FORWARD TRANSFORM
	*/
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "FORWARD WORK = [");
	for(int i = 0; i < LENWRK; i++) {
		fprintf(file, "%f, ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PERFORM BACKWARD TRANSFORM
	*/
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "BACKWARD WORK = [");
	for(int i = 0; i < LENWRK; i++) {
		fprintf(file, "%f ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PRINT TEST RESULTS
	*/
	diff = 0.0;
	for(int i = 0; i < N; i++) { diff = max(diff, abs(C[i].r - NEED.r) + abs(C[i].i - NEED.i)); }
	printf("CFFT1 FORWARD-BACKWARD MAX ERROR = %f\n", diff);


    fclose(file);
    return 0;
}
