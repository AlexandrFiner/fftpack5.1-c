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
	complex NEED[1000];
	real diff = 0.0;

	// Back-forward
	/**
	--- IDENTIFY TEST AND INITIALIZE FFT
	*/
	printf("Program cfft1i and related messages:\n");
	cfft1i_(&N, WSAVE, &LENSAV, &IER);

    fprintf(file, "Back-forward\n");
    fprintf(file, "IDENTIFY TEST AND INITIALIZE FFT\n");

	/**
	--- GENERATE TEST VECTOR FOR BACKWARD-FORWARD TEST
	*/
	for(int i = 1; i < N; i++) {
		RR = (float)rand() / RAND_MAX;
		RI = (float)rand() / RAND_MAX;
		C[i].r = RR;
		C[i].i = RI;
		NEED[i].r = RR;
		NEED[i].i = RI;
	}
    fprintf(file, "C = [");
	for(int i = 1; i < N; i++) {
		fprintf(file, "%f+%fj ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");

	/**
	--- PERFORM BACKWARD TRANSFORM
	*/
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "BACKWARD WORK = [");
	for(int i = 1; i < LENWRK; i++) {
		fprintf(file, "%f, ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PERFORM FORWARD TRANSFORM
	*/
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "FORWARD WORK = [");
	for(int i = 1; i < LENWRK; i++) {
		fprintf(file, "%f, ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PRINT TEST RESULTS
	*/
	diff = 0.0;
	for(int i = 1; i < N; i++) { diff = max(diff, abs(C[i].r - NEED[i].r) + abs(C[i].i - NEED[i].i)); }
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

	/**
	GENERATE TEST VECTOR FOR FORWARD-BACKWARD TEST
	*/
	for(int i = 1; i < N; i++) {
		RR = (float)rand() / RAND_MAX;
		RI = (float)rand() / RAND_MAX;
		C[i].r = RR;
		C[i].i = RI;
		NEED[i].r = RR;
		NEED[i].i = RI;
	}
    fprintf(file, "C = [");
	for(int i = 1; i < N; i++) {
		fprintf(file, "%f+%fj ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");

	/**
	--- PERFORM FORWARD TRANSFORM
	*/
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "FORWARD WORK = [");
	for(int i = 1; i < LENWRK; i++) {
		fprintf(file, "%f, ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PERFORM BACKWARD TRANSFORM
	*/
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "BACKWARD WORK = [");
	for(int i = 1; i < LENWRK; i++) {
		fprintf(file, "%f ", WORK[i]);
	}
	fprintf(file, "]\n");

	/**
	--- PRINT TEST RESULTS
	*/
	diff = 0.0;
	for(int i = 1; i < N; i++) { diff = max(diff, abs(C[i].r - NEED[i].r) + abs(C[i].i - NEED[i].i)); }
	printf("CFFT1 FORWARD-BACKWARD MAX ERROR = %f\n", diff);


    fclose(file);
    return 0;
}
