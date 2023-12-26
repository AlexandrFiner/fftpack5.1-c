# test_fftpack5.py
import unittest
import numpy as np
from ctypes import CDLL, c_int, c_void_p, c_float, c_double

# Загрузка сконвертированной библиотеки C
fftpack5 = CDLL('./fftpack5.so')  # Замените './fftpack5.so' на путь к вашему .so файлу

# Объявление сигнатур функций
cfft1i = fftpack5.CFFT1I
cfft1i.argtypes = [c_int, c_void_p, c_void_p]
cfft1i.restype = None

cfft1b = fftpack5.CFFT1B
cfft1b.argtypes = [c_int, c_void_p, c_void_p]
cfft1b.restype = None

cfft1f = fftpack5.CFFT1F
cfft1f.argtypes = [c_int, c_void_p, c_void_p]
cfft1f.restype = None

class TestFFTPack5(unittest.TestCase):
    def test_cfft1i(self):
        n = 8
        x = np.random.rand(n) + 1j * np.random.rand(n)
        w = np.zeros_like(x)

        cfft1i(c_int(n), x.ctypes.data_as(c_void_p), w.ctypes.data_as(c_void_p))

        # Напишите здесь тестовые случаи для CFFT1I
        self.assertTrue(np.allclose(np.fft.ifft(x), w))

    def test_cfft1b(self):
        n = 8
        x = np.random.rand(n) + 1j * np.random.rand(n)
        w = np.zeros_like(x)

        cfft1b(c_int(n), x.ctypes.data_as(c_void_p), w.ctypes.data_as(c_void_p))

        # Напишите здесь тестовые случаи для CFFT1B
        self.assertTrue(np.allclose(np.fft.fft(x), w))

    def test_cfft1f(self):
        n = 8
        x = np.random.rand(n) + 1j * np.random.rand(n)
        w = np.zeros_like(x)

        cfft1f(c_int(n), x.ctypes.data_as(c_void_p), w.ctypes.data_as(c_void_p))

        # Напишите здесь тестовые случаи для CFFT1F
        self.assertTrue(np.allclose(np.fft.fft(x), w))

if __name__ == '__main__':
    unittest.main()
