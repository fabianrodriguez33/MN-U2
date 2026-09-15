"""Descomposición LU (algoritmo de Doolittle) y resolución de sistemas Ax = b."""

import numpy as np

PIVOT_TOL = 1e-12


class SingularMatrixError(Exception):
    """Se lanza cuando un pivote u_ii es nulo o casi nulo durante la descomposición."""

    def __init__(self, i):
        self.i = i
        super().__init__(
            f"División por cero detectada: u_{i+1}{i+1} ≈ 0. "
            "La matriz requiere pivoteo o permutar filas."
        )


class LUSolver:
    """Encapsula la descomposición de Doolittle y las sustituciones asociadas."""

    @staticmethod
    def doolittle_decompose(A):
        A = np.array(A, dtype=float)
        n = A.shape[0]
        if A.shape[0] != A.shape[1]:
            raise ValueError("La matriz A debe ser cuadrada.")

        L = np.eye(n)
        U = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):
                U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])

            if abs(U[i, i]) < PIVOT_TOL:
                raise SingularMatrixError(i)

            for j in range(i + 1, n):
                L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]

        return L, U

    @staticmethod
    def forward_substitution(L, b):
        L = np.array(L, dtype=float)
        b = np.array(b, dtype=float).flatten()
        n = L.shape[0]
        y = np.zeros(n)

        for i in range(n):
            y[i] = b[i] - np.dot(L[i, :i], y[:i])

        return y

    @staticmethod
    def backward_substitution(U, y):
        U = np.array(U, dtype=float)
        y = np.array(y, dtype=float).flatten()
        n = U.shape[0]
        x = np.zeros(n)

        for i in range(n - 1, -1, -1):
            if abs(U[i, i]) < PIVOT_TOL:
                raise SingularMatrixError(i)
            x[i] = (y[i] - np.dot(U[i, i + 1:], x[i + 1:])) / U[i, i]

        return x

    @classmethod
    def solve(cls, L, U, b):
        y = cls.forward_substitution(L, b)
        x = cls.backward_substitution(U, y)
        return y, x
