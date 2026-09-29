# Prompt / Instrucciones para Claude Code: Aplicación Streamlit de Interpolación Polinómica de Lagrange (Sesión 8)

## 🎯 Objetivo General
Integrar el módulo de la **Sesión 8 (Interpolación Polinómica de Lagrange)** en la aplicación web interactiva en **Streamlit** (`MN-U2`). La aplicación debe permitir ingresar conjuntos dinámicos de puntos discretos $(x_i, y_i)$, calcular los polinomios base $L_k(x)$, evaluar cualquier valor $x_{	ext{eval}}$ deseado, construir la representación simbólica del polinomio interpolante $P_n(x)$, y graficar la curva interpolada junto con los puntos experimentales.

---

## 🏗️ Arquitectura del Proyecto (`D:\MN-U2`)

```text
MN-U2/
├── app.py                      # UI Streamlit con menú de navegación (Sesión 6, 7 y 8)
├── lu_solver.py                # Módulo Sesión 6: Factorización LU (Doolittle)
├── iterative_solvers.py        # Módulo Sesión 7: Jacobi y Gauss-Seidel
├── lagrange_solver.py          # NUEVO MÓDULO Sesión 8: Interpolación de Lagrange
├── requisitos.txt              # streamlit, numpy, pandas, matplotlib, sympy
└── .claude/launch.json         # Servidor dev (Streamlit en puerto 8501)
```

---

## 📐 Fundamento Teórico y Procedimiento de Cólculo

### 1. Polinomio Interpolante de Lagrange
Dado un conjunto de $n+1$ puntos distintos $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$, existe un único polinomio $P_n(x)$ de grado menor o igual a $n$ tal que $P_n(x_i) = y_i$:

$$P_n(x) = \sum_{k=0}^{n} y_k \cdot L_k(x)$$

### 2. Polinomios Base de Lagrange
Para cada $k \in \{0, 1, \dots, n\}$, el polinomio base $L_k(x)$ se define como:

$$L_k(x) = \prod_{\substack{j=0 \ j 
eq k}}^{n} rac{x - x_j}{x_k - x_j} = rac{(x - x_0)(x - x_1)\cdots(x - x_{k-1})(x - x_{k+1})\cdots(x - x_n)}{(x_k - x_0)(x_k - x_1)\cdots(x_k - x_{k-1})(x_k - x_{k+1})\cdots(x_k - x_n)}$$

### 3. Propiedad de Ortogonalidad Implícita
$$L_k(x_i) = \delta_{ki} = egin{cases} 1 & 	ext{si } k = i \ 0 & 	ext{si } k 
eq i \end{cases}$$

---

## 📊 Caso Benchmark / Test de Verificación (Guía Autónoma Sesión 8)

### Contexto del Problema: Optimización de Latencia en Microservicios Cloud
Un equipo DevOps analiza la relación entre la **Memoria RAM asignada $x$ (GB)** a un contenedor y la **Latencia Media $y = f(x)$ (milisegundos)** de las peticiones HTTP.

### Tabla de Datos Experimentales ($n=3$, Grado $\le 3$):

| Punto ($k$) | Memoria $x_k$ (GB) | Latencia $y_k$ (ms) |
| :---: | :---: | :---: |
| $x_0$ | 2.0 | 150.0 |
| $x_1$ | 4.0 | 85.0 |
| $x_2$ | 8.0 | 50.0 |
| $x_3$ | 12.0 | 70.0 |

**Punto a evaluar ($x_{	ext{eval}}$):** $6.0 	ext{ GB}$

### Paso a Paso del Procedimiento Numérico ($x = 6.0$):

1. **Polinomios Base $L_k(6.0)$:**
   * $L_0(6) = rac{(6-4)(6-8)(6-12)}{(2-4)(2-8)(2-12)} = rac{24}{-120} = -0.200000$
   * $L_1(6) = rac{(6-2)(6-8)(6-12)}{(4-2)(4-8)(4-12)} = rac{48}{64} = 0.750000$
   * $L_2(6) = rac{(6-2)(6-4)(6-12)}{(8-2)(8-4)(8-12)} = rac{-48}{-96} = 0.500000$
   * $L_3(6) = rac{(6-2)(6-4)(6-8)}{(12-2)(12-4)(12-8)} = rac{-16}{320} = -0.050000$

2. **Evaluación de la Solución $P_3(6.0)$:**

| $k$ | $x_k$ | $y_k$ | $L_k(6.0)$ | Término $y_k \cdot L_k(6.0)$ |
| :---: | :---: | :---: | :---: | :---: |
| 0 | 2.0 | 150.0 | -0.200000 | -30.000000 |
| 1 | 4.0 | 85.0 | 0.750000 | +63.750000 |
| 2 | 8.0 | 50.0 | 0.500000 | +25.000000 |
| 3 | 12.0 | 70.0 | -0.050000 | -3.500000 |
| **Suma** | -- | -- | **1.000000** | **$P_3(6.0) = 55.250000 	ext{ ms}$** |

---

## 🛠️ Requisitos Técnicos e Implementación en Python

### 1. Módulo Numérico (`lagrange_solver.py`)
Implementar la clase `LagrangeSolver`:
* `compute_basis(x_points, x_eval)`: Devuelve la lista/array de valores $[L_0(x_{	ext{eval}}), \dots, L_n(x_{	ext{eval}})]$.
* `evaluate(x_points, y_points, x_eval)`: Devuelve el valor interpolado $P_n(x_{	ext{eval}})$.
* `get_polynomial_expression(x_points, y_points)`: Usa `sympy` o formateo explícito para construir la expresión simplificada del polinomio $P_n(x)$ (ej. $P(x) = a_n x^n + \dots + a_0$).
* **Validación de Errores:** Lanzar un error personalizado `DuplicateXError` si existen puntos con valores de $x_i$ repetidos (evita división por cero en $L_k$).

### 2. Interfaz en Streamlit (`app.py`)
* **Menú Lateral (`st.sidebar`):**
  * Radio button para seleccionar la sesión:
    1. **Sesión 6:** Factorización LU (Doolittle)
    2. **Sesión 7:** Métodos Iterativos (Jacobi y Gauss-Seidel)
    3. **Sesión 8:** Interpolación Polinómica de Lagrange (NUEVO)
* **Contenido de la Sesión 8:**
  * **Configuración de Datos:**
    * Selector dinámico del número de puntos ($n+1 \ge 2$).
    * Tabla editable (`st.data_editor`) para ingresar o modificar libremente las coordenadas $(x_i, y_i)$.
    * Input numérico para ingresar el valor $x_{	ext{eval}}$.
    * Botón **"Cargar Caso Predeterminado (Sesión 8 - Microservicios)"** para autocompletar la tabla con el benchmark ($2, 4, 8, 12$ GB y $x_{	ext{eval}}=6$).
  * **Resultados e Interfaz:**
    * Muestra la tabla de coeficientes $L_k(x_{	ext{eval}})$ y términos ponderados en una tabla interactiva.
    * Muestra el polinomio interpolante simplificado $P_n(x)$ renderizado en **LaTeX** con `st.latex`.
    * Muestra la estimación final $P_n(x_{	ext{eval}})$ en una tarjeta de métrica (`st.metric`).
  * **Visualización Gráfica (`matplotlib`):**
    * Graficar los puntos discretos experimentales $(x_i, y_i)$ como puntos destacados (markers).
    * Graficar la curva continua $P_n(x)$ generada con alta resolución (ej. 200 puntos en el intervalo $[x_{\min}-1, x_{\max}+1]$).
    * Marcar de forma diferenciada el punto interpolado $(x_{	ext{eval}}, P_n(x_{	ext{eval}}))$ con líneas punteadas hacia los ejes.

---

## 🚀 Instrucción de Ejecución para Claude Code
> "Integra el módulo de la Sesión 8 (Interpolación Polinómica de Lagrange) en la aplicación Streamlit `D:\MN-U2`. Crea el archivo `lagrange_solver.py` con la lógica numérica pura y amplía `app.py` añadiendo un selector de módulo en la barra lateral para alternar entre las Sesiones 6, 7 y 8. La UI debe permitir el ingreso dinámico de $n+1$ puntos $(x_i, y_i)$, la carga del caso predeterminado de la guía autónoma (RAM vs Latencia), la visualización de la tabla de polinomios de base $L_k(x)$, la expresión simbólica en LaTeX del polinomio, la métrica del valor interpolado $P_n(x_{	ext{eval}})$ y un gráfico comparativo con Matplotlib."
