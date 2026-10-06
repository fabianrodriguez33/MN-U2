"""Sesión 9: Interpolación polinómica de Newton mediante diferencias divididas."""

import numpy as np
import sympy as sp

from lagrange_solver import DuplicateXError


class NewtonSolver:
    @staticmethod
    def _validate(x, y):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.ndim != 1 or len(x) < 2:
            raise ValueError("Se requieren al menos 2 puntos.")
        if len(y) != len(x):
            raise ValueError("x e y deben tener la misma longitud.")
        vals, counts = np.unique(x, return_counts=True)
        if np.any(counts > 1):
            raise DuplicateXError(vals[counts > 1][0])
        return x, y

    @staticmethod
    def compute_divided_differences(x, y):
        """Matriz triangular n×n: columna k = diferencias divididas de orden k."""
        x, y = NewtonSolver._validate(x, y)
        n = len(x)
        table = np.full((n, n), np.nan)
        table[:, 0] = y
        for k in range(1, n):
            for i in range(n - k):
                table[i, k] = (table[i + 1, k - 1] - table[i, k - 1]) / (x[i + k] - x[i])
        return table

    @staticmethod
    def coefficients(table):
        """Primera fila de la tabla: a_k = f[x_0, ..., x_k]."""
        return np.asarray(table)[0, :].astype(float)

    @staticmethod
    def evaluate(x_eval, x_nodes, coefs):
        """Evalúa P_n(x_eval) con el esquema anidado (Horner) de Newton."""
        x_eval = np.asarray(x_eval, dtype=float)
        result = np.full_like(x_eval, coefs[-1], dtype=float)
        for k in range(len(coefs) - 2, -1, -1):
            result = coefs[k] + (x_eval - x_nodes[k]) * result
        return float(result) if result.ndim == 0 else result

    @staticmethod
    def get_symbolic_expressions(x_nodes, coefs):
        """Devuelve (P_n en forma de Newton factorizada, P_n canónico, símbolo x)."""
        x = sp.Symbol("x")
        xs = [sp.nsimplify(float(v), rational=True) for v in x_nodes]
        cs = [sp.nsimplify(float(c), rational=True, tolerance=1e-9) for c in coefs]
        newton = 0
        prod = 1
        for k, c in enumerate(cs):
            newton += c * prod
            if k < len(xs):
                prod = prod * (x - xs[k])
        return newton, sp.expand(newton), x

    @staticmethod
    def newton_latex(x_nodes, coefs, digits=6):
        """P_n(x) en forma de Newton, con términos factorizados (x - x_j)."""
        out = ""
        for k, c in enumerate(coefs):
            if k > 0 and abs(c) < 1e-12:
                continue
            factors = "".join(
                f"(x - {sp.latex(sp.nsimplify(float(x_nodes[j]), rational=True))})"
                for j in range(k)
            )
            body = sp.latex(sp.Float(abs(float(c)), digits)) + factors
            if not out:
                out = ("-" if c < 0 else "") + body
            else:
                out += (" - " if c < 0 else " + ") + body
        return "P(x) = " + out

    @staticmethod
    def canonical_latex(poly, x, digits=6):
        """P_n(x) simplificado en LaTeX, con coeficientes redondeados."""
        rounded = sum(
            sp.Float(float(c), digits) * x ** int(m[0])
            for m, c in sp.Poly(poly, x).terms()
        )
        return "P(x) = " + sp.latex(rounded)
