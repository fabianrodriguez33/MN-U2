"""Métodos iterativos (Jacobi y Gauss-Seidel) para resolver sistemas Ax = b."""

import numpy as np
import pandas as pd


def check_diagonal_dominance(A):
    """Retorna (dominante_por_fila: list[bool], es_edd_global: bool)."""
    A = np.array(A, dtype=float)
    n = A.shape[0]
    per_row = []
    for i in range(n):
        diag = abs(A[i, i])
        off_sum = np.sum(np.abs(A[i, :])) - diag
        per_row.append(bool(diag > off_sum))
    return per_row, bool(all(per_row))


def check_sassenfeld(A):
    """Retorna (betas: np.ndarray, converge: bool) según el criterio de Sassenfeld."""
    A = np.array(A, dtype=float)
    n = A.shape[0]
    beta = np.zeros(n)
    for i in range(n):
        s_prev = sum(abs(A[i, j]) * beta[j] for j in range(i))
        s_next = sum(abs(A[i, j]) for j in range(i + 1, n))
        beta[i] = (s_prev + s_next) / abs(A[i, i])
    return beta, bool(np.max(beta) < 1)


def _relative_error(x_new, x_old):
    diff = np.abs(x_new - x_old)
    denom = np.abs(x_new)
    denom[denom < 1e-15] = 1e-15
    return np.max(diff / denom) * 100.0


def jacobi_solver(A, b, x0, tol=1e-4, max_iter=100):
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).flatten()
    x = np.array(x0, dtype=float).flatten().copy()
    n = A.shape[0]

    history = [{"k": 0, **{f"x{i+1}": x[i] for i in range(n)}, "Error (%)": np.nan}]
    converged = False
    it = 0

    for it in range(1, max_iter + 1):
        x_new = np.zeros(n)
        for i in range(n):
            s = np.dot(A[i, :], x) - A[i, i] * x[i]
            x_new[i] = (b[i] - s) / A[i, i]

        error = _relative_error(x_new, x)
        history.append({"k": it, **{f"x{i+1}": x_new[i] for i in range(n)}, "Error (%)": error})
        x = x_new

        if error < tol * 100.0:
            converged = True
            break

    r = b - A @ x
    return {
        "history": pd.DataFrame(history),
        "x": x,
        "r": r,
        "iterations": it,
        "converged": converged,
    }


def gauss_seidel_solver(A, b, x0, tol=1e-4, max_iter=100):
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).flatten()
    x = np.array(x0, dtype=float).flatten().copy()
    n = A.shape[0]

    history = [{"k": 0, **{f"x{i+1}": x[i] for i in range(n)}, "Error (%)": np.nan}]
    converged = False
    it = 0

    for it in range(1, max_iter + 1):
        x_new = x.copy()
        for i in range(n):
            s_prev = np.dot(A[i, :i], x_new[:i])
            s_next = np.dot(A[i, i + 1:], x[i + 1:])
            x_new[i] = (b[i] - s_prev - s_next) / A[i, i]

        error = _relative_error(x_new, x)
        history.append({"k": it, **{f"x{i+1}": x_new[i] for i in range(n)}, "Error (%)": error})
        x = x_new

        if error < tol * 100.0:
            converged = True
            break

    r = b - A @ x
    return {
        "history": pd.DataFrame(history),
        "x": x,
        "r": r,
        "iterations": it,
        "converged": converged,
    }
