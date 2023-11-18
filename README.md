# https://netlib.org/clapack/readme.maintain

Задание выполняли Мазяров Александр КС-20-1Б и Морохин Александр КС-20-1Б

Конвертирование библиотеки fftpack5 в С с помощью f2c и написание генератора тестов на IntelPython.
https://box.icmm.ru/index.php/s/4vhnShGHH9WkLkg



gcc main.c -o test -lfftpack -L libs -I include -v

ranlib libs/libfftpack.a

gcc main.c -o test -lm -lfftpack -L libs -I include/