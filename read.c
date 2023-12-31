#include <stdio.h>
#include <stdlib.h>

// Function to read float array from a binary file
float* readFloatArray(FILE *file, size_t *size) {
    fread(size, sizeof(size_t), 1, file);  // Read the size of the array
    float *floatArray = (float *)malloc(sizeof(float) * (*size));
    fread(floatArray, sizeof(float), *size, file);  // Read the array data
    return floatArray;
}

int main() {
    // Specify the file path
    const char *file_path = "float_data.bin";

    // Read binary data from file
    FILE *file = fopen(file_path, "rb");
    if (file != NULL) {
        size_t size1, size2;
        float *readArray1 = readFloatArray(file, &size1);
        float *readArray2 = readFloatArray(file, &size2);
        fclose(file);

        // Print the resulting arrays of floats
        printf("Array 1:\n");
        for (size_t i = 0; i < size1; i++) {
            printf("%f ", readArray1[i]);
        }
        printf("\n");

        printf("Array 2:\n");
        for (size_t i = 0; i < size2; i++) {
            printf("%f ", readArray2[i]);
        }
        printf("\n");

        // Free allocated memory
        free(readArray1);
        free(readArray2);

    } else {
        printf("Error opening file for reading.\n");
        return 1;
    }

    return 0;
}
