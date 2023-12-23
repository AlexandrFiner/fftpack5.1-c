# import numpy as np
#
# def initialize_complex(size):
#     return np.random.randn(size) + 1j * np.random.randn(size)
#
# def backward_complex(arr):
#     return np.imag(arr)
#
# def forward_complex(arr):
#     return np.abs(arr)
#
# # Пример использования
# size = 5
#
# # Инициализация
# complex_array = initialize_complex(size)
# print("Initialized Complex Array:")
# print(complex_array)
#
# # Прямой проход
# forward_result = forward_complex(complex_array)
# print("\nForward Result:")
# print(forward_result)
#
# # Обратный проход
# backward_result = backward_complex(complex_array)
# print("\nBackward Result:")
# print(backward_result)
#
#
import numpy as np

def initialize_complex(size):
    # Заменяем случайную инициализацию статическими данными
    real_part = np.arange(1, size + 1)
    imaginary_part = np.arange(size, 0, -1)
    return real_part + 1j * imaginary_part

def backward_complex(arr):
    return np.imag(arr)

def forward_complex(arr):
    return np.abs(arr)

# Пример использования
size = 5

# Инициализация
complex_array = initialize_complex(size)
print("Initialized Complex Array:")
print(complex_array)

# Прямой проход
forward_result = forward_complex(complex_array)
print("\nForward Result:")
print(forward_result)

# Обратный проход
backward_result = backward_complex(complex_array)
print("\nBackward Result:")
print(backward_result)
