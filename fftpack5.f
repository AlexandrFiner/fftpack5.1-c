! fftpack5.f90

SUBROUTINE CFFT1I(n, x, w)
  INTEGER, INTENT(IN) :: n
  COMPLEX, INTENT(IN) :: x(n)
  COMPLEX, INTENT(OUT) :: w(n)
  INTEGER :: i, j, k, m
  COMPLEX :: t, u
  INTEGER :: i1, j1, i2, j2
  REAL :: theta, sin_theta, cos_theta

  m = INT(LOG(FLOAT(n))/LOG(2.0))
  j = 1

  DO i = 1, n-1
    IF (i .LT. j) THEN
      t = x(j)
      x(j) = x(i)
      x(i) = t
    ENDIF

    k = n / 2

    DO WHILE (k .GE. 1 .AND. k .LT. j)
      j = j - k
      k = k / 2
    END DO

    j = j + k
  END DO

  DO i = 1, m
    theta = -2.0 * REAL(i) * 3.141592653589793238 / FLOAT(n)
    sin_theta = SIN(theta)
    cos_theta = COS(theta)

    DO k = 1, n, 2 * 2**i
      DO j = k, k + 2**i - 1
        i1 = j
        i2 = j + 2**i
        j1 = j + 2**(i-1)
        j2 = j1 + 2**i

        t = x(i2) * CMPLX(cos_theta, sin_theta)
        u = x(i1)
        x(i1) = u + t
        x(i2) = u - t
      END DO
    END DO
  END DO

  w = x
END SUBROUTINE

SUBROUTINE CFFT1B(n, x, w)
  INTEGER, INTENT(IN) :: n
  COMPLEX, INTENT(IN) :: x(n)
  COMPLEX, INTENT(OUT) :: w(n)
  INTEGER :: i, j, k, m
  COMPLEX :: t, u
  INTEGER :: i1, j1, i2, j2
  REAL :: theta, sin_theta, cos_theta

  m = INT(LOG(FLOAT(n))/LOG(2.0))
  j = 1

  DO i = 1, n-1
    IF (i .LT. j) THEN
      t = x(j)
      x(j) = x(i)
      x(i) = t
    ENDIF

    k = n / 2

    DO WHILE (k .GE. 1 .AND. k .LT. j)
      j = j - k
      k = k / 2
    END DO

    j = j + k
  END DO

  DO i = 1, m
    theta = 2.0 * REAL(i) * 3.141592653589793238 / FLOAT(n)
    sin_theta = SIN(theta)
    cos_theta = COS(theta)

    DO k = 1, n, 2 * 2**i
      DO j = k, k + 2**i - 1
        i1 = j
        i2 = j + 2**i
        j1 = j + 2**(i-1)
        j2 = j1 + 2**i

        t = x(i2) * CMPLX(cos_theta, sin_theta)
        u = x(i1)
        x(i1) = u + t
        x(i2) = u - t
      END DO
    END DO
  END DO

  w = x / FLOAT(n)
END SUBROUTINE

SUBROUTINE CFFT1F(n, x, w)
  INTEGER, INTENT(IN) :: n
  COMPLEX, INTENT(IN) :: x(n)
  COMPLEX, INTENT(OUT) :: w(n)
  INTEGER :: i, j, k, m
  COMPLEX :: t, u
  INTEGER :: i1, j1, i2, j2
  REAL :: theta, sin_theta, cos_theta

  m = INT(LOG(FLOAT(n))/LOG(2.0))
  j = 1

  DO i = 1, n-1
    IF (i .LT. j) THEN
      t = x(j)
      x(j) = x(i)
      x(i) = t
    ENDIF

    k = n / 2

    DO WHILE (k .GE. 1 .AND. k .LT. j)
      j = j - k
      k = k / 2
    END DO

    j = j + k
  END DO

  DO i = 1, m
    theta = -2.0 * REAL(i) * 3.141592653589793238 / FLOAT(n)
    sin_theta = SIN(theta)
    cos_theta = COS(theta)

    DO k = 1, n, 2 * 2**i
      DO j = k, k + 2**i - 1
        i1 = j
        i2 = j + 2**i
        j1 = j + 2**(i-1)
        j2 = j1 + 2**i

        t = x(i2) * CMPLX(cos_theta, sin_theta)
        u = x(i1)
        x(i1) = u + t
        x(i2) = u - t
      END DO
    END DO
  END DO

  w = x
END SUBROUTINE
