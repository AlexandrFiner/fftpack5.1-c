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
    FILE *file_complex = fopen("complex_data_combined.bin", "wb");
	FILE *file = fopen("result.txt", "w");

    // Проверить, был ли файл успешно открыт
    if (file == NULL) {
        fprintf(stderr, "Ошибка открытия файла.\n");
        return 1;  // Вернуть код ошибки
    }

    integer I, N, LENSAV, IER, LENWRK;
    N = 10;
    LENSAV = 2 * N + (int)(log2((double)N)) + 4;
    LENWRK = 2 * N;

    real RR, RI;


    real *WSAVE = (real *)malloc(LENSAV * sizeof(real));
    real *WORK = (real *)malloc(LENWRK * sizeof(real));

	srand(time(NULL));

	static integer INC = 1;
    complex *C = (complex *)malloc(N * sizeof(complex));
    complex *NEED = (complex *)malloc(N * sizeof(complex));
	real diff = 0.0;

	// Back-forward
	/**
	--- IDENTIFY TEST AND INITIALIZE FFT
	*/
	printf("Program cfft1i and related messages:\n");
	cfft1i_(&N, WSAVE, &LENSAV, &IER);
    fprintf(file, "WSAVE = [");
	for(int i = 0; i < LENSAV; i++) {
		fprintf(file, "%f ", WSAVE[i]);
	}
	// WSAVE[5] = 1.0;
	fprintf(file, "]\n");


    fprintf(file, "Back-forward\n");
    fprintf(file, "IDENTIFY TEST AND INITIALIZE FFT\n");

	/**
	--- GENERATE TEST VECTOR FOR BACKWARD-FORWARD TEST
	*/
	for(int i = 0; i < N; i++) {
		RR = (float)rand() / RAND_MAX;
		RI = 0.0;
		C[i].r = RR;
		C[i].i = RI;
		NEED[i].r = RR;
		NEED[i].i = RI;
	}
    fprintf(file, "Generated C = [");
	for(int i = 0; i < N; i++) {
		fprintf(file, "%f+%fj, ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");
    fwrite(C, sizeof(complex), N, file_complex);

	/**
	--- PERFORM BACKWARD TRANSFORM
	*/
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "Backward C = [");
	for(int i = 0; i < N; i++) {
		fprintf(file, "%f+%fj, ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");
	fwrite(C, sizeof(complex), N, file_complex);

	/**
	--- PERFORM FORWARD TRANSFORM
	*/
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "Forward C = [");
	for(int i = 0; i < N; i++) {
		fprintf(file, "%f+%fj, ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");

	/**
	--- PRINT TEST RESULTS
	*/
	diff = 0.0;
	for(int i = 0; i < N; i++) { diff = max(diff, abs(C[i].r - NEED[i].r) + abs(C[i].i - NEED[i].i)); }
	printf("CFFT1 BACKWARD-FORWARD MAX ERROR = %.20f\n", diff);

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
	for(int i = 0; i < N; i++) {
		RR = (float)rand() / RAND_MAX;
		RI = 0.0;
		C[i].r = RR;
		C[i].i = RI;
		NEED[i].r = RR;
		NEED[i].i = RI;
	}
    fprintf(file, "Generated C = [");
	for(int i = 0; i < N; i++) {
		fprintf(file, "%f+%fj, ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");
    fwrite(C, sizeof(complex), N, file_complex);

	/**
	--- PERFORM FORWARD TRANSFORM
	*/
	cfft1f_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "Forward C = [");
	for(int i = 0; i < N; i++) {
		fprintf(file, "%f+%fj, ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");
    fwrite(C, sizeof(complex), N, file_complex);

	/**
	--- PERFORM BACKWARD TRANSFORM
	*/
	cfft1b_(&N, &INC, C, &N, WSAVE, &LENSAV, WORK, &LENWRK, &IER);
    fprintf(file, "Backward C = [");
	for(int i = 0; i < N; i++) {
		fprintf(file, "%f+%fj, ", C[i].r, C[i].i);
	}
	fprintf(file, "]\n");

	/**
	--- PRINT TEST RESULTS
	*/
	diff = 0.0;
	for(int i = 0; i < N; i++) { diff = max(diff, abs(C[i].r - NEED[i].r) + abs(C[i].i - NEED[i].i)); }
	printf("CFFT1 FORWARD-BACKWARD MAX ERROR = %.20f\n", diff);


    fclose(file);
    return 0;
}
