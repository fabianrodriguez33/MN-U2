# Prompt / Instrucciones para Claude Code: Validación y Prueba de Simulacro de Examen 2 (S11 GT - Métodos Numéricos)

## 🎯 Objetivo General
Ejecutar una suite de **pruebas de validación y verificación** sobre la aplicación Streamlit y los módulos numéricos del proyecto `D:\MN-U2` (`lu_solver.py`, `iterative_solvers.py`, `lagrange_solver.py`, `newton_solver.py`), utilizando como benchmark los **6 ejercicios del Simulacro de Examen 2 (`u2_s11_gt_interpolación`)**. 

---

## 🛠️ Estructura del Proyecto a Validar
```text
D:\MN-U2/
├── app.py                  # UI Streamlit (Sesiones 6, 7, 8, 9)
├── lu_solver.py            # Factorización LU Doolittle
├── iterative_solvers.py    # Jacobi, Gauss-Seidel, EDD, Sassenfeld
├── lagrange_solver.py      # Interpolación de Lagrange
├── newton_solver.py        # Diferencias Divididas y Newton-Gregory
├── test_simulacro_s11.py   # NUEVO: Script de prueba automatizada para los 6 ejercicios
├── README.md
└── requirements.txt
```

---

## 📐 Benchmark del Simulacro (Ejercicios, Procedimiento y Resultados Esperados)

### 1. Sistema de Distribución de Agua (Método de Jacobi)
* **Contexto:** Caudales $Q_1, Q_2, Q_3$ en red de tuberías.
* **Sistema $A\mathbf{Q} = \mathbf{B}$:**
  $$\begin{pmatrix} 10 & -2 & 1 \ -3 & 9 & 2 \ 2 & -1 & 7 \end{pmatrix} \begin{pmatrix} Q_1 \\ Q_2 \\ Q_3 \end{pmatrix} = \begin{pmatrix} 27 \\ -8 \\ 15 \end{pmatrix}$$
* **Verificación EDD (Dominancia Diagonal Estricta):**
  * Fila 1: $|10| > |-2| + |1| = 3$ (Cumple)
  * Fila 2: $|9| > |-3| + |2| = 5$ (Cumple)
  * Fila 3: $|7| > |2| + |-1| = 3$ (Cumple)
* **Iteraciones de Jacobi ($\mathbf{Q}^{(0)} = [0, 0, 0]^T$):**

| Iteración $k$ | $Q_1^{(k)}$ | $Q_2^{(k)}$ | $Q_3^{(k)}$ |
| :---: | :---: | :---: | :---: |
| **0** | 0.000000 | 0.000000 | 0.000000 |
| **1** | 2.700000 | -0.888889 | 2.142857 |
| **2** | 2.307937 | -0.465079 | 1.244444 |
| **3** | 2.482540 | -0.396120 | 1.417007 |
| **Exacto** | **2.488889** | **-0.365812** | **1.379487** |

---

### 2. Sistema Lineal (Gauss-Seidel y Criterio de Sassenfeld)
* **Sistema $A\mathbf{x} = \mathbf{b}$:**
  $$\begin{pmatrix} 10 & 1 & 2 \ 4 & 6 & -1 \ -2 & 3 & 8 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \\ x_3 \end{pmatrix} = \begin{pmatrix} 3 \\ 9 \\ 51 \end{pmatrix}$$
* **Criterio de Sassenfeld:**
  $$\beta_1 = \frac{|1| + |2|}{|10|} = 0.3000, \quad \beta_2 = \frac{|4|(0.3000) + |-1|}{|6|} = 0.3667, \quad \beta_3 = \frac{|-2|(0.3000) + |3|(0.3667)}{|8|} = 0.2125$$
  * $\max(\beta_i) = 0.3667 < 1 \implies$ **Garantiza convergencia.**
* **Iteraciones Gauss-Seidel ($\mathbf{x}^{(0)} = [0, 0, 0]^T$, Tol $< 1\%$):**

| Iteración $k$ | $x_1^{(k)}$ | $x_2^{(k)}$ | $x_3^{(k)}$ | Error Relativo Máx |
| :---: | :---: | :---: | :---: | :---: |
| **1** | 0.300000 | 1.300000 | 5.962500 | 100.0% |
| **2** | -1.022500 | 3.175417 | 4.928594 | 129.3% |
| **3** | -1.003260 | 2.990273 | 5.002833 | 6.19% |
| **4** | **-0.999594** | **3.000201** | **5.000026** | **0.37% (< 1%)** |
| **Exacto** | **-1.000000** | **3.000000** | **5.000000** | **0.00%** |

---

### 3. Interpolación de Lagrange (Selección Óptima de Puntos)
* **Datos Tabulados:** $x \in [1.6, 2.0, 2.5, 3.2, 4.0, 4.5]$, $f(x) \in [2, 8, 14, 15, 8, 2]$.
* **Evaluación en $x = 2.8$:**
  * **Grado 2 (3 puntos más cercanos: $x = [2.0, 2.5, 3.2]$):** $P_2(2.8) = \mathbf{15.485714}$
  * **Grado 3 (4 puntos más cercanos: $x = [2.0, 2.5, 3.2, 4.0]$):** $P_3(2.8) = \mathbf{15.388571}$

---

### 4. Deflexión de Cable de Puente Colgante (Newton Diferencias Divididas)
* **Datos Experimentales:** Posición $x$ (m): $[0, 20, 40, 60, 80]$, Deflexión $f(x)$ (cm): $[0, 15, 22, 35, 60]$.
* **Tabla de Diferencias Divididas:**
  * $f[x_0] = 0.0$
  * $f[x_0, x_1] = 0.75$
  * $f[x_0, x_1, x_2] = -0.01$
  * $f[x_0, x_1, x_2, x_3] = 0.000291667$ ($rac{7}{24000}$)
  * $f[x_0, x_1, x_2, x_3, x_4] = -0.000002083$ ($-rac{1}{480000}$)
* **Estimaciones:**
  * $f(50) = \mathbf{27.1875	ext{ cm}}$
  * $f(70) = \mathbf{45.9375	ext{ cm}}$

---

### 5. Asentamiento de Zapata Superficial (Newton-Gregory Progresivo)
* **Datos (Intervalos iguales $h = 2$):** Tiempo $t \in [0, 2, 4, 6]$, Asentamiento $y \in [0.0, 1.6, 3.5, 6.3]$.
* **Diferencias Finitas Progresivas:**
  * $\Delta^0 y_0 = 0.0, \quad \Delta^1 y_0 = 1.6, \quad \Delta^2 y_0 = 0.3, \quad \Delta^3 y_0 = 0.6$
* **Polinomio Interpolante:**
  $$P_3(t) = 0.8t + 0.0375t(t-2) + 0.0125t(t-2)(t-4) = 0.0125t^3 - 0.0375t^2 + 0.825t$$
* **Estimación a las 5 horas ($t = 5$):**
  $$P_3(5) = \mathbf{4.7500	ext{ mm}}$$

---

### 6. Cota Superior del Error de Truncamiento
* **Función:** $f(x) = rac{1}{1+x}$ en $x_0 = 0, x_1 = 1, x_2 = 2$. Evaluar en $x = 1.3$.
* **Cota Teórica del Error:**
  $$|E_2(x)| \le rac{|(x-x_0)(x-x_1)(x-x_2)|}{3!} \max_{t \in [0,2]} |f'''(t)|$$
  * $f'''(t) = -rac{6}{(1+t)^4} \implies \max_{[0,2]} |f'''(t)| = |f'''(0)| = 6.0$
  * Para $x = 1.3$: $|(1.3 - 0)(1.3 - 1)(1.3 - 2)| = |1.3 	imes 0.3 	imes (-0.7)| = 0.273$
  * **Cota de Error:** $|E_2(1.3)| \le rac{0.273}{6} 	imes 6.0 = \mathbf{0.2730}$

---

## 🚀 Instrucciones de Ejecución para Claude Code

1. **Crear Script de Test Unitario (`test_simulacro_s11.py`):**
   * Crear un script independiente que importe los módulos `iterative_solvers.py`, `lagrange_solver.py` y `newton_solver.py`.
   * Ejecutar y validar programáticamente que las salidas coincidan exactamente con los valores de benchmark especificados en este documento.
2. **Verificar Interfaz Streamlit (`app.py`):**
   * Iniciar la app con `streamlit run app.py` y probar la carga manual de los 6 casos de estudio.
   * Asegurar que las tablas de iteración, expresiones LaTeX y gráficos generados reflejen los mismos resultados teóricos.
3. **Reporte:**
   * Mostrar en consola la confirmación de aprobación (`PASSED`) de los 6 ejercicios del simulacro.
