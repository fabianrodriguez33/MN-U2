"""Sesión 10: Interpolación de Newton-Gregory (diferencias finitas, nodos equidistantes)."""

import math

import numpy as np
import sympy as sp

from lagrange_solver import DuplicateXError


class UnequalSpacingError(Exception):
    """Se lanza cuando los nodos no están igualmente espaciados."""

    def __init__(self):
        super().__init__(
            "Error: Newton-Gregory requiere nodos igualmente espaciados "
            "(x_{i+1} - x_i = h constante). Use diferencias divididas (Sesión 9) en otro caso."
        )


class GregorySolver:
    @staticmethod
    def validate(x, y):
        """Valida los datos y devuelve (x, y, h); exige nodos crecientes y equidistantes."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.ndim != 1 or len(x) < 2:
            raise ValueError("Se requieren al menos 2 puntos.")
        if len(y) != len(x):
            raise ValueError("x e y deben tener la misma longitud.")
        vals, counts = np.unique(x, return_counts=True)
        if np.any(counts > 1):
            raise DuplicateXError(vals[counts > 1][0])
        steps = np.diff(x)
        h = steps[0]
        if h <= 0 or not np.allclose(steps, h, rtol=1e-9, atol=1e-12):
            raise UnequalSpacingError()
        return x, y, float(h)

    @staticmethod
    def difference_table(y):
        """Matriz n×n: columna k = diferencias finitas de orden k (Δ^k y_i = ∇^k y_{i+k})."""
        y = np.asarray(y, dtype=float)
        n = len(y)
        table = np.full((n, n), np.nan)
        table[:, 0] = y
        for k in range(1, n):
            for i in range(n - k):
                table[i, k] = table[i + 1, k - 1] - table[i, k - 1]
        return table

    @staticmethod
    def forward_coefficients(table):
        """Δ^k y_0 (primera fila de la tabla)."""
        return np.asarray(table)[0, :].astype(float)

    @staticmethod
    def backward_coefficients(table):
        """∇^k y_n (última diagonal: elemento [n-k, k])."""
        t = np.asarray(table)
        n = t.shape[0]
        return np.array([t[n - 1 - k, k] for k in range(n)], dtype=float)

    @staticmethod
    def _s_terms(s, k, mode):
        """Coeficiente binomial generalizado: Π(s - j)/k! (progresivo) o Π(s + j)/k! (regresivo)."""
        sign = -1 if mode == "forward" else 1
        out = np.ones_like(s, dtype=float)
        for j in range(k):
            out = out * (s + sign * j)
        return out / math.factorial(k)

    @staticmethod
    def evaluate(x_eval, x, coefs, h, mode="forward"):
        """P_n(x_eval) con Newton-Gregory progresivo (anclado en x_0) o regresivo (en x_n)."""
        x_eval = np.asarray(x_eval, dtype=float)
        anchor = x[0] if mode == "forward" else x[-1]
        s = (x_eval - anchor) / h
        total = np.zeros_like(s, dtype=float)
        for k, c in enumerate(coefs):
            total = total + c * GregorySolver._s_terms(s, k, mode)
        return float(total) if total.ndim == 0 else total

    @staticmethod
    def get_symbolic_expressions(x, coefs, h, mode="forward"):
        """Devuelve (P_n en forma de Newton-Gregory, P_n canónico, símbolo x)."""
        sym = sp.Symbol("x")
        anchor = sp.nsimplify(float(x[0] if mode == "forward" else x[-1]), rational=True)
        hh = sp.nsimplify(float(h), rational=True)
        s = (sym - anchor) / hh
        sign = -1 if mode == "forward" else 1
        poly = 0
        for k, c in enumerate(coefs):
            term = sp.nsimplify(float(c), rational=True, tolerance=1e-9)
            for j in range(k):
                term *= (s + sign * j)
            poly += term / sp.factorial(k)
        return poly, sp.expand(poly), sym

    @staticmethod
    def newton_gregory_latex(coefs, mode="forward", digits=6):
        """P_n(x) = Σ (Δ^k y_0 o ∇^k y_n) · C(s, k) en LaTeX, con s = (x - ancla)/h."""
        sign = "-" if mode == "forward" else "+"
        out = ""
        for k, c in enumerate(coefs):
            if k > 0 and abs(c) < 1e-12:
                continue
            factors = "".join(f"(s {sign} {j})" if j else "s" for j in range(k))
            body = sp.latex(sp.Float(abs(float(c)), digits))
            if k:
                body += f"\\frac{{{factors}}}{{{k}!}}" if k > 1 else factors
            out += ("-" if c < 0 else "") + body if not out else (" - " if c < 0 else " + ") + body
        return f"P(x) = {out}"

    @staticmethod
    def canonical_latex(poly, x, digits=6):
        rounded = sum(
            sp.Float(float(c), digits) * x ** int(m[0])
            for m, c in sp.Poly(poly, x).terms()
        )
        return "P(x) = " + sp.latex(rounded)
