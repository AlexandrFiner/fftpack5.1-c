import struct

# Specify the file path
file_path = 'float_data.bin'

# Read binary data from file
with open(file_path, 'rb') as file:
    binary_data = file.read()

# Determine the length of each array
len_float_array1 = len(binary_data) // (2 * struct.calcsize('f'))
len_float_array2 = len_float_array1

# Unpack binary data into arrays
float_array1 = struct.unpack('f' * len_float_array1, binary_data[:len_float_array1 * struct.calcsize('f')])
float_array2 = struct.unpack('f' * len_float_array2, binary_data[len_float_array1 * struct.calcsize('f'):])

# Print the resulting arrays of floats
print(float_array1)
print(float_array2)
