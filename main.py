import numpy as np
from numpy.fft import fft


def generate_test_data(size):
    # Генерируем случайные данные для теста
    data = np.random.random(size)
    return data


def perform_fft(data):
    result = fft(data)
    return result


def main():
    # Генерируем тестовые данные и выполняем fft
    data_size = 1024
    test_data = generate_test_data(data_size)

    fft_result = perform_fft(test_data)

    print("Исходные данные:")
    print(test_data)
    print("\nРезультат FFT:")
    print(fft_result)


if __name__ == "__main__":
    main()
