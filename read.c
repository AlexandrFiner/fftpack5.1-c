#include <stdio.h>

// Define the complex type in C
typedef struct {
    double real;  // real part
    double imag;  // imaginary part
} Complex;

int main() {
    // Define and initialize C[100] array
    Complex C[100];
    for (int i = 0; i < 100; ++i) {
        C[i].real = i * 0.1;  // Replace with your desired values
        C[i].imag = i * 0.2;
    }

    // Define and initialize C2[10] array
    Complex C2[10];
    for (int i = 0; i < 10; ++i) {
        C2[i].real = 1.0;  // Replace with your desired values
        C2[i].imag = 0.0;
    }

    // Write both arrays to a single binary file
    FILE *file = fopen("complex_data_combined.bin", "wb");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        return 1;
    }

    // Write C[100] array
    fwrite(C, sizeof(Complex), 100, file);

    // Write C2[10] array
    fwrite(C2, sizeof(Complex), 10, file);

    fclose(file);

    return 0;
}
