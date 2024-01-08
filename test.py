import struct

# Define the complex type in Python
class ComplexStruct:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

# Read both arrays from the combined binary file
with open("complex_data_combined.bin", "rb") as file:
    data_combined = file.read()

# Calculate the number of complex values
num_values_total = len(data_combined) // struct.calcsize("ff")

# Split the data into four equal parts
part_size = num_values_total // 4
data_parts = [data_combined[i:i+part_size*8] for i in range(0, len(data_combined), part_size*8)]

# Process each part separately
for part_num, part_data in enumerate(data_parts):
    complex_array_part = [ComplexStruct(*struct.unpack("ff", part_data[i:i+8])) for i in range(0, len(part_data), 8)]

    # Print the results for each part
    print(f"Part {part_num + 1}:")
    for i, c in enumerate(complex_array_part):
        print(f"Element {i + 1}: Real = {c.real}, Imaginary = {c.imag}")
