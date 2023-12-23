import unittest
import numpy as np
from fftpack5_1 import CFFT1I, CFFT1B, CFFT1F

class TestCFFT1(unittest.TestCase):

    def test_1d_complex_initialization(self):
        # Тестирование 1D комплексной инициализации
        N = 8  # Пример размера массива
        WSAVE = np.zeros(2 * N + 4, dtype=np.float64)  # Пример массива WSAVE
        IER = CFFT1I(N, WSAVE)
        self.assertEqual(IER, 0)  # Проверка успешной инициализации

    def test_1d_complex_backward_forward(self):
        # Тестирование 1D комплексного обратного и прямого преобразований
        N = 8  # Пример размера массива
        WSAVE = np.zeros(2 * N + 4, dtype=np.float64)  # Пример массива WSAVE
        C = np.random.rand(N) + 1j * np.random.rand(N)  # Пример комплексного входного массива

        # Сохранение оригинального массива для сравнения
        original_C = np.copy(C)

        # Обратное преобразование
        IER = CFFT1B(N, 1, C, N, WSAVE)
        self.assertEqual(IER, 0)  # Проверка успешного обратного преобразования

        # Прямое преобразование
        IER = CFFT1F(N, 1, C, N, WSAVE)
        self.assertEqual(IER, 0)  # Проверка успешного прямого преобразования

        # Проверка, что результат прямого-обратного преобразования близок к оригиналу
        np.testing.assert_allclose(C, original_C, atol=1e-8)

if __name__ == '__main__':
    unittest.main()
