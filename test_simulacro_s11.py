"""Validación del Simulacro de Examen 2 (S11 GT) contra los módulos del proyecto."""

import sys

import numpy as np

from iterative_solvers import check_diagonal_dominance, check_sassenfeld, gauss_seidel_solver, jacobi_solver
from gregory_solver import GregorySolver
from lagrange_solver import LagrangeSolver
from newton_solver import NewtonSolver


def ej1_jacobi():
    A = [[10, -2, 1], [-3, 9, 2], [2, -1, 7]]
    b = [27, -8, 15]
    per_row, edd = check_diagonal_dominance(A)
    assert edd and all(per_row)
    res = jacobi_solver(A, b, [0, 0, 0], tol=1e-12, max_iter=200)
    h = res["history"]
    esperado = {
        1: [2.700000, -0.888889, 2.142857],
        2: [2.307937, -0.465079, 1.244444],
        3: [2.482540, -0.396120, 1.417007],
    }
    for k, row in esperado.items():
        got = h.loc[h["k"] == k, ["x1", "x2", "x3"]].to_numpy()[0]
        assert np.allclose(got, row, atol=5e-7), (k, got, row)
    assert np.allclose(res["x"], np.linalg.solve(A, b), atol=1e-9)
    assert np.allclose(res["x"], [2.488889, -0.365812, 1.379487], atol=5e-7), res["x"]


def ej2_gauss_seidel():
    A = [[10, 1, 2], [4, 6, -1], [-2, 3, 8]]
    b = [3, 9, 51]
    beta, ok = check_sassenfeld(A)
    assert np.allclose(beta, [0.3, 0.36667, 0.2125], atol=1e-4) and ok
    res = gauss_seidel_solver(A, b, [0, 0, 0], tol=0.01)  # tol 1 %
    h = res["history"]
    esperado = {
        1: [0.3, 1.3, 5.9625],
        2: [-1.0225, 3.175417, 4.928594],
        3: [-1.003260, 2.990273, 5.002833],
        4: [-0.999594, 3.000201, 5.000026],
    }
    for k, row in esperado.items():
        got = h.loc[h["k"] == k, ["x1", "x2", "x3"]].to_numpy()[0]
        assert np.allclose(got, row, atol=5e-6), (k, got, row)
    assert res["converged"] and res["iterations"] == 4, res["iterations"]
    assert np.allclose(res["x"], [-1, 3, 5], atol=1e-3)


def ej3_lagrange():
    xs = [1.6, 2.0, 2.5, 3.2, 4.0, 4.5]
    ys = [2, 8, 14, 15, 8, 2]
    p2 = LagrangeSolver.evaluate(xs[1:4], ys[1:4], 2.8)
    p3 = LagrangeSolver.evaluate(xs[1:5], ys[1:5], 2.8)
    assert abs(p2 - 15.485714) < 5e-7, p2
    assert abs(p3 - 15.388571) < 5e-7, p3


def ej4_newton_dd():
    x = [0, 20, 40, 60, 80]
    y = [0, 15, 22, 35, 60]
    table = NewtonSolver.compute_divided_differences(x, y)
    c = NewtonSolver.coefficients(table)
    assert np.allclose(c, [0, 0.75, -0.01, 7 / 24000, -1 / 480000], atol=1e-9), c
    assert abs(NewtonSolver.evaluate(50, x, c) - 27.1875) < 1e-9
    assert abs(NewtonSolver.evaluate(70, x, c) - 45.9375) < 1e-9


def ej5_newton_gregory():
    t = np.array([0, 2, 4, 6], dtype=float)
    y = np.array([0.0, 1.6, 3.5, 6.3])
    t, y, h = GregorySolver.validate(t, y)
    assert h == 2.0
    table = GregorySolver.difference_table(y)
    fwd = GregorySolver.forward_coefficients(table)
    assert np.allclose(fwd, [0.0, 1.6, 0.3, 0.6]), fwd
    bwd = GregorySolver.backward_coefficients(table)
    assert np.allclose(bwd, [6.3, 2.8, 0.9, 0.6]), bwd
    assert abs(GregorySolver.evaluate(5.0, t, fwd, h, "forward") - 4.75) < 1e-12
    assert abs(GregorySolver.evaluate(5.0, t, bwd, h, "backward") - 4.75) < 1e-12
    # Nota: ins11.md daba 0.875t en la forma expandida; la expansión correcta es 0.825t.
    for tv in (0.0, 1.0, 3.7, 5.0):
        ref = 0.0125 * tv**3 - 0.0375 * tv**2 + 0.825 * tv
        assert abs(GregorySolver.evaluate(tv, t, fwd, h, "forward") - ref) < 1e-12
    # Coincide con Newton por diferencias divididas
    c = NewtonSolver.coefficients(NewtonSolver.compute_divided_differences(t, y))
    assert abs(NewtonSolver.evaluate(5.0, t, c) - 4.75) < 1e-12


def ej6_cota_error():
    f = lambda v: 1 / (1 + v)
    nodes = [0, 1, 2]
    prod = abs(np.prod([1.3 - xi for xi in nodes]))
    assert abs(prod - 0.273) < 1e-12
    ts = np.linspace(0, 2, 20001)
    m3 = np.max(np.abs(-6 / (1 + ts) ** 4))
    assert abs(m3 - 6.0) < 1e-12
    cota = prod / 6 * m3
    assert abs(cota - 0.2730) < 1e-12
    # El error real debe respetar la cota
    ys = [f(v) for v in nodes]
    real = abs(f(1.3) - LagrangeSolver.evaluate(nodes, ys, 1.3))
    assert real <= cota, (real, cota)


TESTS = [
    ("1. Jacobi (distribución de agua)", ej1_jacobi),
    ("2. Gauss-Seidel + Sassenfeld", ej2_gauss_seidel),
    ("3. Lagrange (selección de puntos)", ej3_lagrange),
    ("4. Newton diferencias divididas", ej4_newton_dd),
    ("5. Newton-Gregory progresivo", ej5_newton_gregory),
    ("6. Cota del error de truncamiento", ej6_cota_error),
]

if __name__ == "__main__":
    fallos = 0
    for nombre, fn in TESTS:
        try:
            fn()
            print(f"PASSED  {nombre}")
        except AssertionError as e:
            fallos += 1
            print(f"FAILED  {nombre}: {e}")
    sys.exit(1 if fallos else 0)
