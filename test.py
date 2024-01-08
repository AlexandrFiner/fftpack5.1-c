import struct

# Define the complex type in Python
class ComplexStruct:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

# Read both arrays from the combined binary file
with open("complex_data_combined.bin", "rb") as file:
    data_combined = file.read()

# Calculate the number of values in each array
num_values_C = 10
num_values_C2 = 10

# Read C[100] array from the combined binary file
complex_array_C = [ComplexStruct(*struct.unpack("dd", data_combined[i:i+16])) for i in range(0, num_values_C * 16, 16)]

# Read C2[10] array from the combined binary file
complex_array_C2 = [ComplexStruct(*struct.unpack("dd", data_combined[i:i+16])) for i in range(num_values_C * 16, len(data_combined), 16)]

# Print the results for C[100] array
print("C[10] array:")
for i, c in enumerate(complex_array_C):
    print(f"Element {i + 1}: Real = {c.real}, Imaginary = {c.imag}")

# Print the results for C2[10] array
print("\nC2[10] array:")
for i, c in enumerate(complex_array_C2):
    print(f"Element {i + 1}: Real = {c.real}, Imaginary = {c.imag}")
