import struct

# Sample arrays of floats
float_array1 = [1.23, 4.56, 7.89, 10.11]
float_array2 = [2.34, 5.67, 8.90, 12.34]
float_array3 = [8.90, 12.34,7.89, 10.11]



array = [float_array1,float_array2,float_array3]

for arr in range(len(array)):
    b_data = struct.pack('f'* len(array[arr]),*array[arr])
    file_path = f'float_data{arr}.bin'
    with open("tests/"+file_path, 'wb') as file:
        file.write(b_data)
    print(f"file {arr} was created successfully")






