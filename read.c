#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main() {
    // Specify the file path
    const char *file_path = "float_data.bin";

    // Open the file in binary read mode
    FILE *file = fopen(file_path, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Determine the size of the file
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    // Allocate memory to store the binary data
    uint8_t *binary_data = (uint8_t *)malloc(file_size);
    if (binary_data == NULL) {
        perror("Memory allocation error");
        fclose(file);
        return 1;
    }

    // Read the binary data from the file
    fread(binary_data, 1, file_size, file);
    fclose(file);

    // Determine the number of floats in the binary data
    size_t num_floats = file_size / sizeof(float);

    // Unpack binary data into floats
    float *float_array = (float *)binary_data;

    // Print the resulting array of floats
    for (size_t i = 0; i < num_floats; ++i) {
        printf("%.2f\n", float_array[i]);
    }

    // Free allocated memory
    free(binary_data);

    return 0;
}
