"""Sesión 8: Interpolación polinómica de Lagrange."""

import numpy as np
import sympy as sp


class DuplicateXError(Exception):
    """Se lanza cuando existen valores x_i repetidos."""

    def __init__(self, value):
        self.value = value
        super().__init__(
            f"Error: el valor x = {value:g} está repetido. Los nodos x_i deben ser "
            "distintos para evitar división por cero en L_k(x)."
        )


class LagrangeSolver:
    @staticmethod
    def _validate(x_points, y_points=None):
        x = np.asarray(x_points, dtype=float)
        if x.ndim != 1 or len(x) < 2:
            raise ValueError("Se requieren al menos 2 puntos.")
        if y_points is not None and len(np.asarray(y_points)) != len(x):
            raise ValueError("x e y deben tener la misma longitud.")
        vals, counts = np.unique(x, return_counts=True)
        if np.any(counts > 1):
            raise DuplicateXError(vals[counts > 1][0])
        return x

    @staticmethod
    def compute_basis(x_points, x_eval):
        """Devuelve [L_0(x_eval), ..., L_n(x_eval)]."""
        x = LagrangeSolver._validate(x_points)
        basis = np.ones(len(x))
        for k in range(len(x)):
            for j in range(len(x)):
                if j != k:
                    basis[k] *= (x_eval - x[j]) / (x[k] - x[j])
        return basis

    @staticmethod
    def evaluate(x_points, y_points, x_eval):
        """Devuelve P_n(x_eval)."""
        LagrangeSolver._validate(x_points, y_points)
        basis = LagrangeSolver.compute_basis(x_points, x_eval)
        return float(np.dot(np.asarray(y_points, dtype=float), basis))

    @staticmethod
    def get_polynomial_expression(x_points, y_points):
        """Devuelve (P_n(x) expandido como expresión sympy, símbolo x)."""
        LagrangeSolver._validate(x_points, y_points)
        xs = [sp.Rational(str(float(v))) for v in x_points]
        ys = [sp.Rational(str(float(v))) for v in y_points]
        x = sp.Symbol("x")
        poly = 0
        for k in range(len(xs)):
            term = ys[k]
            for j in range(len(xs)):
                if j != k:
                    term *= (x - xs[j]) / (xs[k] - xs[j])
            poly += term
        return sp.expand(poly), x

    @staticmethod
    def polynomial_latex(x_points, y_points, digits=6):
        """P_n(x) simplificado en LaTeX, con coeficientes redondeados."""
        poly, x = LagrangeSolver.get_polynomial_expression(x_points, y_points)
        rounded = sum(
            sp.Float(float(c), digits) * x**int(m[0])
            for m, c in sp.Poly(poly, x).terms()
        )
        return "P(x) = " + sp.latex(rounded)
