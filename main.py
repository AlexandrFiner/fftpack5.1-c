import struct

# Sample arrays of floats
float_array1 = [1.23, 4.56, 7.89, 10.11]
float_array2 = [2.34, 5.67, 8.90, 12.34]

# Pack floats from both arrays into binary data
binary_data = struct.pack('f' * len(float_array1), *float_array1) + struct.pack('f' * len(float_array2), *float_array2)

# Specify the file path
file_path = 'float_data.bin'

# Write binary data to file
with open(file_path, 'wb') as file:
    file.write(binary_data)
