"""Aplicación Streamlit: Sesión 6 (LU Doolittle) y Sesión 7 (Jacobi / Gauss-Seidel)
para balanceo de carga en clústeres Cloud."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from iterative_solvers import (
    check_diagonal_dominance,
    check_sassenfeld,
    gauss_seidel_solver,
    jacobi_solver,
)
from lu_solver import LUSolver, SingularMatrixError

st.set_page_config(page_title="Métodos Numéricos - MN-U2", layout="wide")

# --- Caso de prueba predeterminado (Sesión 6) ---
DEFAULT_A = np.array([
    [4, 2, 1],
    [12, 10, 5],
    [-8, 8, 7],
], dtype=float)
DEFAULT_B1 = np.array([14, 46, 26], dtype=float)
DEFAULT_B2 = np.array([20, 62, 30], dtype=float)

# --- Caso de prueba predeterminado (Sesión 7) ---
DEFAULT_A7 = np.array([
    [10, -2, -1, 0],
    [-1, 8, 0, -2],
    [-2, 0, 12, -3],
    [0, -1, -2, 9],
], dtype=float)
DEFAULT_B7 = np.array([15, 18, 25, 20], dtype=float)


def matrix_to_latex(M, name):
    rows = " \\\\ ".join(" & ".join(f"{v:.4g}" for v in row) for row in M)
    return f"{name} = \\begin{{pmatrix}} {rows} \\end{{pmatrix}}"


def vector_to_latex(v, name):
    rows = " \\\\ ".join(f"{val:.4g}" for val in v)
    return f"{name} = \\begin{{pmatrix}} {rows} \\end{{pmatrix}}"


# ============================================================
# SESIÓN 6 — Factorización LU (Doolittle)
# ============================================================

def init_state_s6(n):
    st.session_state.A_df = pd.DataFrame(
        np.zeros((n, n)),
        columns=[f"a{j+1}" for j in range(n)],
        index=[f"f{i+1}" for i in range(n)],
    )
    st.session_state.b_df = pd.DataFrame(
        np.zeros((n, 1)), columns=["b"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.n = n


def load_preset_s6():
    n = 3
    st.session_state.n = n
    st.session_state.A_df = pd.DataFrame(
        DEFAULT_A,
        columns=[f"a{j+1}" for j in range(n)],
        index=[f"f{i+1}" for i in range(n)],
    )
    st.session_state.b_df = pd.DataFrame(
        DEFAULT_B1, columns=["b"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.pop("lu_result", None)


def render_sesion6():
    st.title("🖥️ Balanceo de Carga en Clústeres Cloud — Factorización LU (Doolittle)")
    st.markdown(
        "Esta aplicación resuelve sistemas de ecuaciones lineales $A\\mathbf{x} = \\mathbf{b}$ "
        "mediante la **factorización LU (algoritmo de Doolittle)**, útil para modelar el "
        "balanceo de carga entre nodos de un clúster Cloud: cada ecuación representa la "
        "distribución de tráfico/recursos entre microservicios, y resolver el sistema permite "
        "calcular la asignación óptima ($\\mathbf{x}$) para un vector de demanda ($\\mathbf{b}$) dado."
    )

    with st.sidebar:
        st.header("⚙️ Configuración (Sesión 6)")

        if "n" not in st.session_state:
            init_state_s6(3)

        n = st.number_input(
            "Dimensión del sistema (N)", min_value=2, max_value=10,
            value=st.session_state.n, step=1, key="n_input_s6"
        )
        if n != st.session_state.n:
            init_state_s6(n)

        if st.button("📌 Cargar Ejemplo Predeterminado (Sesión 6)", use_container_width=True):
            load_preset_s6()

        st.divider()
        st.caption(
            "Validación de pivotes: si $|u_{ii}| < 10^{-12}$ durante la descomposición, "
            "la aplicación detiene el cálculo y muestra un error amigable."
        )

    n = st.session_state.n

    st.subheader("1️⃣ Matriz de coeficientes A y vector de carga b")
    col_a, col_b = st.columns([3, 1])

    with col_a:
        st.markdown("**Matriz A**")
        A_edited = st.data_editor(
            st.session_state.A_df, key="A_editor", num_rows="fixed", use_container_width=True
        )
        st.session_state.A_df = A_edited

    with col_b:
        st.markdown("**Vector b (inicial)**")
        b_edited = st.data_editor(
            st.session_state.b_df, key="b_editor", num_rows="fixed", use_container_width=True
        )
        st.session_state.b_df = b_edited

    A = st.session_state.A_df.to_numpy(dtype=float)
    b1 = st.session_state.b_df.to_numpy(dtype=float).flatten()

    st.divider()

    if st.button("🚀 Calcular Descomposición LU y Resolver", type="primary"):
        try:
            L, U = LUSolver.doolittle_decompose(A)
            y, x = LUSolver.solve(L, U, b1)
            st.session_state.lu_result = {"L": L, "U": U, "y": y, "x": x, "b": b1}
        except SingularMatrixError as e:
            st.session_state.pop("lu_result", None)
            st.error(str(e))
        except ValueError as e:
            st.session_state.pop("lu_result", None)
            st.error(f"Error de validación: {e}")

    if "lu_result" in st.session_state:
        res = st.session_state.lu_result
        L, U, y, x, b_used = res["L"], res["U"], res["y"], res["x"], res["b"]

        st.subheader("2️⃣ Resultados de la Descomposición")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Matriz $L$ (Triangular Inferior)**")
            st.latex(matrix_to_latex(L, "L"))
            st.dataframe(pd.DataFrame(L).round(6), use_container_width=True)
        with c2:
            st.markdown("**Matriz $U$ (Triangular Superior)**")
            st.latex(matrix_to_latex(U, "U"))
            st.dataframe(pd.DataFrame(U).round(6), use_container_width=True)

        st.subheader("3️⃣ Proceso de Sustitución")
        c3, c4 = st.columns(2)
        with c3:
            st.markdown("**Sustitución hacia adelante** ($L\\mathbf{y} = \\mathbf{b}$)")
            st.latex(vector_to_latex(y, "\\mathbf{y}"))
        with c4:
            st.markdown("**Sustitución hacia atrás** ($U\\mathbf{x} = \\mathbf{y}$)")
            st.latex(vector_to_latex(x, "\\mathbf{x}"))

        st.subheader("4️⃣ Vector Solución Final")
        st.latex(vector_to_latex(x, "\\mathbf{x}"))
        st.dataframe(
            pd.DataFrame({"Variable": [f"x{i+1}" for i in range(len(x))], "Valor": np.round(x, 6)}),
            use_container_width=True, hide_index=True,
        )

        st.divider()
        st.subheader("5️⃣ Reevaluación rápida con un nuevo vector b (reutilizando L y U)")
        st.caption("No se recalcula la matriz LU: solo se repiten las sustituciones hacia adelante y hacia atrás.")

        default_b2 = DEFAULT_B2 if len(x) == 3 else np.zeros(len(x))
        b2_df = st.data_editor(
            pd.DataFrame(default_b2, columns=["b2"], index=[f"f{i+1}" for i in range(len(x))]),
            key="b2_editor", num_rows="fixed", use_container_width=True,
        )

        if st.button("🔁 Reevaluar con b2"):
            try:
                b2 = b2_df.to_numpy(dtype=float).flatten()
                y2, x2 = LUSolver.solve(L, U, b2)
                st.session_state.lu_result2 = {"y": y2, "x": x2}
            except SingularMatrixError as e:
                st.error(str(e))

        if "lu_result2" in st.session_state:
            y2, x2 = st.session_state.lu_result2["y"], st.session_state.lu_result2["x"]
            c5, c6 = st.columns(2)
            with c5:
                st.markdown("**Vector $Y_2$**")
                st.latex(vector_to_latex(y2, "\\mathbf{y}_2"))
            with c6:
                st.markdown("**Vector $X_2$ (Nueva solución)**")
                st.latex(vector_to_latex(x2, "\\mathbf{x}_2"))


# ============================================================
# SESIÓN 7 — Métodos Iterativos (Jacobi y Gauss-Seidel)
# ============================================================

def init_state_s7(n):
    st.session_state.A7_df = pd.DataFrame(
        np.zeros((n, n)),
        columns=[f"a{j+1}" for j in range(n)],
        index=[f"f{i+1}" for i in range(n)],
    )
    st.session_state.b7_df = pd.DataFrame(
        np.zeros((n, 1)), columns=["b"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.x0_df = pd.DataFrame(
        np.zeros((n, 1)), columns=["x0"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.n7 = n


def load_preset_s7():
    n = 4
    st.session_state.n7 = n
    st.session_state.A7_df = pd.DataFrame(
        DEFAULT_A7,
        columns=[f"a{j+1}" for j in range(n)],
        index=[f"f{i+1}" for i in range(n)],
    )
    st.session_state.b7_df = pd.DataFrame(
        DEFAULT_B7, columns=["b"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.x0_df = pd.DataFrame(
        np.zeros(n), columns=["x0"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.pop("iter_result", None)


def render_sesion7():
    st.title("🖥️ Balanceo de Carga en Clústeres Cloud — Métodos Iterativos (Jacobi / Gauss-Seidel)")
    st.markdown(
        "Esta aplicación resuelve sistemas de ecuaciones lineales $A\\mathbf{x} = \\mathbf{b}$ "
        "mediante los **métodos iterativos de Jacobi y Gauss-Seidel**, comparando su velocidad "
        "de convergencia para modelar el balanceo de carga entre nodos de un clúster Cloud."
    )

    with st.sidebar:
        st.header("⚙️ Configuración (Sesión 7)")

        if "n7" not in st.session_state:
            init_state_s7(4)

        n = st.number_input(
            "Dimensión del sistema (N)", min_value=2, max_value=10,
            value=st.session_state.n7, step=1, key="n_input_s7"
        )
        if n != st.session_state.n7:
            init_state_s7(n)

        if st.button("📌 Cargar Ejemplo Predeterminado (Sesión 7 - Clúster de 4 Servidores)", use_container_width=True):
            load_preset_s7()

        st.divider()
        tol = st.number_input("Tolerancia de error ε", min_value=1e-8, max_value=1.0, value=1e-4, format="%.8f")
        max_iter = st.number_input("Iteraciones máximas $K_{max}$", min_value=1, max_value=1000, value=100, step=1)

    n = st.session_state.n7

    st.subheader("1️⃣ Matriz A, vector b y vector inicial x⁽⁰⁾")
    col_a, col_b, col_x0 = st.columns([3, 1, 1])

    with col_a:
        st.markdown("**Matriz A**")
        A7_edited = st.data_editor(
            st.session_state.A7_df, key="A7_editor", num_rows="fixed", use_container_width=True
        )
        st.session_state.A7_df = A7_edited

    with col_b:
        st.markdown("**Vector b**")
        b7_edited = st.data_editor(
            st.session_state.b7_df, key="b7_editor", num_rows="fixed", use_container_width=True
        )
        st.session_state.b7_df = b7_edited

    with col_x0:
        st.markdown("**Vector $x^{(0)}$**")
        x0_edited = st.data_editor(
            st.session_state.x0_df, key="x0_editor", num_rows="fixed", use_container_width=True
        )
        st.session_state.x0_df = x0_edited

    A = st.session_state.A7_df.to_numpy(dtype=float)
    b = st.session_state.b7_df.to_numpy(dtype=float).flatten()
    x0 = st.session_state.x0_df.to_numpy(dtype=float).flatten()

    st.divider()

    st.subheader("2️⃣ Verificación de Convergencia")
    per_row, is_edd = check_diagonal_dominance(A)
    beta, sassenfeld_ok = check_sassenfeld(A)

    edd_table = pd.DataFrame({
        "Fila (i)": [i + 1 for i in range(n)],
        "|a_ii|": [abs(A[i, i]) for i in range(n)],
        "Σ|a_ij| (j≠i)": [np.sum(np.abs(A[i, :])) - abs(A[i, i]) for i in range(n)],
        "Cumple EDD": ["✅" if v else "❌" for v in per_row],
    })
    c_edd, c_sass = st.columns(2)
    with c_edd:
        st.markdown("**Dominancia Diagonal Estricta (EDD) por filas**")
        st.dataframe(edd_table, use_container_width=True, hide_index=True)
        if is_edd:
            st.success("La matriz A cumple EDD: convergencia garantizada para Jacobi y Gauss-Seidel.")
        else:
            st.warning(
                "La matriz A NO cumple EDD. La dominancia diagonal es condición suficiente pero "
                "no necesaria: se verifica el criterio de Sassenfeld y la convergencia experimental."
            )
    with c_sass:
        st.markdown("**Criterio de Sassenfeld** (específico para Gauss-Seidel)")
        st.dataframe(
            pd.DataFrame({"i": [i + 1 for i in range(n)], "β_i": np.round(beta, 6)}),
            use_container_width=True, hide_index=True,
        )
        if sassenfeld_ok:
            st.success(f"max(β_i) = {np.max(beta):.6f} < 1: convergencia de Gauss-Seidel garantizada.")
        else:
            st.warning(f"max(β_i) = {np.max(beta):.6f} ≥ 1: Sassenfeld no garantiza convergencia.")

    st.divider()

    if st.button("🚀 Ejecutar Jacobi y Gauss-Seidel", type="primary"):
        jac = jacobi_solver(A, b, x0, tol=tol, max_iter=int(max_iter))
        gs = gauss_seidel_solver(A, b, x0, tol=tol, max_iter=int(max_iter))
        st.session_state.iter_result = {"jacobi": jac, "gs": gs}

    if "iter_result" in st.session_state:
        jac = st.session_state.iter_result["jacobi"]
        gs = st.session_state.iter_result["gs"]

        st.subheader("3️⃣ Tablas Iterativas")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"**Jacobi** — {'convergió' if jac['converged'] else 'NO convergió'} en {jac['iterations']} iteraciones")
            st.dataframe(jac["history"].round(6), use_container_width=True, hide_index=True)
        with c2:
            st.markdown(f"**Gauss-Seidel** — {'convergió' if gs['converged'] else 'NO convergió'} en {gs['iterations']} iteraciones")
            st.dataframe(gs["history"].round(6), use_container_width=True, hide_index=True)

        st.subheader("4️⃣ Comparación de Métodos")
        comp_table = pd.DataFrame({
            "Característica": ["Uso de Datos", "Iteraciones Requeridas", "Paralelización"],
            "Método de Jacobi": [
                "Emplea valores de la iteración anterior (k)",
                f"{jac['iterations']} iteraciones ({'convergió' if jac['converged'] else 'no convergió'})",
                "Alta (componentes independientes)",
            ],
            "Método de Gauss-Seidel": [
                "Emplea valores actualizados inmediatamente (k+1)",
                f"{gs['iterations']} iteraciones ({'convergió' if gs['converged'] else 'no convergió'})",
                "Secuencial (dependencia interna)",
            ],
        })
        st.dataframe(comp_table, use_container_width=True, hide_index=True)

        st.subheader("5️⃣ Gráfico de Convergencia (Error vs. Iteración)")
        fig, ax = plt.subplots()
        jac_err = jac["history"]["Error (%)"].to_numpy()[1:]
        gs_err = gs["history"]["Error (%)"].to_numpy()[1:]
        ax.plot(range(1, len(jac_err) + 1), jac_err, marker="o", label="Jacobi")
        ax.plot(range(1, len(gs_err) + 1), gs_err, marker="s", label="Gauss-Seidel")
        ax.set_xlabel("Iteración (k)")
        ax.set_ylabel("Error relativo (%)")
        ax.set_yscale("log")
        ax.set_title("Convergencia: Jacobi vs. Gauss-Seidel")
        ax.legend()
        ax.grid(True, which="both", linestyle="--", alpha=0.5)
        st.pyplot(fig)

        st.subheader("6️⃣ Vector Solución Final y Residuo")
        c3, c4 = st.columns(2)
        with c3:
            st.markdown("**Jacobi**")
            st.latex(vector_to_latex(jac["x"], "\\mathbf{x}"))
            st.latex(vector_to_latex(jac["r"], "\\mathbf{r}"))
        with c4:
            st.markdown("**Gauss-Seidel**")
            st.latex(vector_to_latex(gs["x"], "\\mathbf{x}"))
            st.latex(vector_to_latex(gs["r"], "\\mathbf{r}"))


# ============================================================
# NAVEGACIÓN PRINCIPAL
# ============================================================

st.sidebar.title("📚 Navegación")
modulo = st.sidebar.radio(
    "Selecciona el módulo",
    ["Sesión 6: Factorización LU (Doolittle)", "Sesión 7: Métodos Iterativos (Jacobi / Gauss-Seidel)"],
)
st.sidebar.divider()

if modulo.startswith("Sesión 6"):
    render_sesion6()
else:
    render_sesion7()
