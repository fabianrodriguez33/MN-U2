"""Aplicación Streamlit: Descomposición LU (Doolittle) para balanceo de carga en clústeres Cloud."""

import numpy as np
import pandas as pd
import streamlit as st

from lu_solver import LUSolver, SingularMatrixError

st.set_page_config(page_title="Descomposición LU - Doolittle", layout="wide")

# --- Caso de prueba predeterminado (Sesión 6) ---
DEFAULT_A = np.array([
    [4, 2, 1],
    [12, 10, 5],
    [-8, 8, 7],
], dtype=float)
DEFAULT_B1 = np.array([14, 46, 26], dtype=float)
DEFAULT_B2 = np.array([20, 62, 30], dtype=float)


def matrix_to_latex(M, name):
    rows = " \\\\ ".join(" & ".join(f"{v:.4g}" for v in row) for row in M)
    return f"{name} = \\begin{{pmatrix}} {rows} \\end{{pmatrix}}"


def vector_to_latex(v, name):
    rows = " \\\\ ".join(f"{val:.4g}" for val in v)
    return f"{name} = \\begin{{pmatrix}} {rows} \\end{{pmatrix}}"


def init_state(n):
    st.session_state.A_df = pd.DataFrame(
        np.zeros((n, n)),
        columns=[f"a{j+1}" for j in range(n)],
        index=[f"f{i+1}" for i in range(n)],
    )
    st.session_state.b_df = pd.DataFrame(
        np.zeros((n, 1)), columns=["b"], index=[f"f{i+1}" for i in range(n)]
    )
    st.session_state.n = n


def load_preset():
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


st.title("🖥️ Balanceo de Carga en Clústeres Cloud — Factorización LU (Doolittle)")
st.markdown(
    "Esta aplicación resuelve sistemas de ecuaciones lineales $A\\mathbf{x} = \\mathbf{b}$ "
    "mediante la **factorización LU (algoritmo de Doolittle)**, útil para modelar el "
    "balanceo de carga entre nodos de un clúster Cloud: cada ecuación representa la "
    "distribución de tráfico/recursos entre microservicios, y resolver el sistema permite "
    "calcular la asignación óptima ($\\mathbf{x}$) para un vector de demanda ($\\mathbf{b}$) dado."
)

# --- Sidebar ---
with st.sidebar:
    st.header("⚙️ Configuración")

    if "n" not in st.session_state:
        init_state(3)

    n = st.number_input(
        "Dimensión del sistema (N)", min_value=2, max_value=10,
        value=st.session_state.n, step=1, key="n_input"
    )
    if n != st.session_state.n:
        init_state(n)

    if st.button("📌 Cargar Ejemplo Predeterminado (Sesión 6)", use_container_width=True):
        load_preset()

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
