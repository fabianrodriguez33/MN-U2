# Prompt / Instrucciones para Claude Code: Aplicación Streamlit de Interpolación Polinómica de Newton mediante Diferencias Divididas (Sesión 9)

## 🎯 Objetivo General
Integrar el módulo de la **Sesión 9: Interpolación Polinómica de Newton mediante Diferencias Divididas** dentro del proyecto existente `D:\MN-U2`. La aplicación web desarrollada en **Python** con **Streamlit** debe permitir conmutar fluidamente entre cuatro sesiones de Métodos Numéricos (Sesión 6, Sesión 7, Sesión 8 y la nueva Sesión 9) a través de la barra lateral, preservando intacto el código de las sesiones previas.

---

## 📂 Estructura del Proyecto Actualizada (`D:\MN-U2`)

```text
D:\MN-U2/
├── app.py                  # Streamlit UI con navegación multi-módulo (Sesiones 6, 7, 8 y 9)
├── lu_solver.py            # Sesión 6: Factorización LU (Doolittle)
├── iterative_solvers.py    # Sesión 7: Métodos Iterativos (Jacobi, Gauss-Seidel, EDD, Sassenfeld)
├── lagrange_solver.py      # Sesión 8: Interpolación de Lagrange
├── newton_solver.py        # NUEVO (Sesión 9): Lógica numérica de Diferencias Divididas
├── README.md               # Actualizado profesionalmente con la documentación de las 4 sesiones
├── .gitignore              # NUEVO: Ignorar __pycache__/, *.pyc, etc.
├── requirements.txt        # Dependencias: streamlit, numpy, pandas, matplotlib, sympy
└── .claude/launch.json     # Configuración del servidor de desarrollo (puerto 8501)
```

---

## 📐 Fundamento Teórico y Algorítmico (Sesión 9)

### 1. Definición de Diferencias Divididas
Dado un conjunto de $n+1$ puntos discretos $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$ con $x_i$ distintos:
* **Orden 0:**
  $$f[x_i] = y_i$$
* **Orden 1:**
  $$f[x_i, x_{i+1}] = \frac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i}$$
* **Orden $k$ general:**
  $$f[x_i, x_{i+1}, \dots, x_{i+k}] = \frac{f[x_{i+1}, \dots, x_{i+k}] - f[x_i, \dots, x_{i+k-1}]}{x_{i+k} - x_i}$$

### 2. Construcción del Polinomio de Newton $P_n(x)$
El polinomio interpolante en la forma de Newton se construye utilizando la primera fila de la Tabla de Diferencias Divididas (coeficientes $a_k = f[x_0, x_1, \dots, x_k]$):
$$P_n(x) = f[x_0] + f[x_0, x_1](x - x_0) + f[x_0, x_1, x_2](x - x_0)(x - x_1) + \dots + f[x_0, x_1, \dots, x_n] \prod_{j=0}^{n-1}(x - x_j)$$

### 3. Excepción de Validación
Si se ingresan puntos con valores $x_i$ repetidos ($x_i = x_j$ para $i \neq j$), el algoritmo debe lanzar una excepción personalizada `DuplicateXError` para evitar divisiones por cero.

---

## 🧪 Benchmark de Validación / Caso Predeterminado (Guía Autónoma Sesión 9)

### Contexto de Ingeniería SRE (Site Reliability Engineering)
Un equipo SRE evalúa la latencia media de respuesta $y$ (en ms) de una API REST crítica desplegada en Cloud según la carga concurrente de peticiones $x$ (en unidades de $100\text{ req/s}$).

### Datos Experimentales ($n=3$, 4 puntos)

| Punto ($i$) | Carga Concurrente $x_i$ ($100\text{ req/s}$) | Latencia Medida $y_i = f(x_i)$ (ms) |
| :---: | :---: | :---: |
| 0 | 1 | 45 |
| 1 | 2 | 65 |
| 2 | 4 | 110 |
| 3 | 7 | 220 |

---

### Procedimiento Numérico y Tabla Completa de Diferencias Divididas

1. **Tabla Triangular de Diferencias Divididas:**

| $i$ | $x_i$ | $y_i = f[x_i]$ | 1ª Dif. $f[x_i, x_{i+1}]$ | 2ª Dif. $f[x_i, x_{i+1}, x_{i+2}]$ | 3ª Dif. $f[x_i, x_{i+1}, x_{i+2}, x_{i+3}]$ |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **0** | **1** | **45.0000** | **20.0000** | **0.8333** ($rac{5}{6}$) | **0.3333** ($rac{1}{3}$) |
| 1 | 2 | 65.0000 | 22.5000 | 2.8333 ($rac{17}{6}$) | -- |
| 2 | 4 | 110.0000 | 36.6667 | -- | -- |
| 3 | 7 | 220.0000 | -- | -- | -- |

2. **Coeficientes Principales (Primera Fila):**
   * $a_0 = f[x_0] = 45$
   * $a_1 = f[x_0, x_1] = 20$
   * $a_2 = f[x_0, x_1, x_2] = \frac{5}{6} \approx 0.833333$
   * $a_3 = f[x_0, x_1, x_2, x_3] = \frac{1}{3} \approx 0.333333$

3. **Polinomio de Newton Explicito:**
   $$P_3(x) = 45 + 20(x-1) + \frac{5}{6}(x-1)(x-2) + \frac{1}{3}(x-1)(x-2)(x-4)$$

4. **Polinomio Simplificado en Forma Canónica Standard:**
   $$P_3(x) = \frac{1}{3}x^3 - \frac{3}{2}x^2 + \frac{133}{6}x + 24 \approx 0.333333x^3 - 1.5x^2 + 22.166667x + 24$$

5. **Interpolación en $x_{\text{eval}} = 5$ ($500\text{ req/s}$):**
   $$P_3(5) = 45 + 20(4) + \frac{5}{6}(4)(3) + \frac{1}{3}(4)(3)(1) = 45 + 80 + 10 + 4 = 139.0\text{ ms}$$

6. **Tabla de Verificación de Puntos Experimentales:**

| Nivel de Carga $x_i$ | Latencia Medida $y_i$ (ms) | $P_3(x_i)$ Evaluado (ms) | Error $|y_i - P_3(x_i)|$ |
| :---: | :---: | :---: | :---: |
| 1.0 | 45.0 | 45.0 | 0.00 |
| 2.0 | 65.0 | 65.0 | 0.00 |
| 4.0 | 110.0 | 110.0 | 0.00 |
| 7.0 | 220.0 | 220.0 | 0.00 |

---

## 🛠️ Requisitos Técnicos e Implementación en Python

### 1. Archivo `newton_solver.py`
Crear el módulo numérico con la clase `NewtonSolver`:
* `compute_divided_differences(x, y)`: Genera la matriz de $n \times n$ de diferencias divididas.
* `evaluate(x_eval, x_nodes, coefs)`: Implementa la regla de Horner / evaluación polinomial iterativa de Newton.
* `get_symbolic_expressions(x_nodes, coefs)`: Utiliza `sympy` para retornar:
  * Expresión de Newton con términos factorizados $(x - x_j)$.
  * Expresión simplificada en forma canónica $a_n x^n + \dots + a_0$.

### 2. Interfaz de Usuario en `app.py`
* **Navegación:** `st.sidebar.radio` para seleccionar la **Sesión 9: Interpolación Polinómica de Newton**.
* **Carga de Datos Flexible:**
  * Selector numérico para la cantidad de puntos ($n+1 \ge 2$).
  * Componente `st.data_editor` que permite ingresar y editar puntos $(x_i, y_i)$ de forma sencilla y directa.
  * Botón **"Cargar Ejemplo Predeterminado (Sesión 9 - SRE REST API)"** para autocompletar la tabla benchmark.
  * Selector / Input numérico para evaluar $x_{\text{eval}}$.
* **Despliegue de Resultados:**
  * Visualización formateada de la **Tabla Completa de Diferencias Divididas** usando `pandas.DataFrame`.
  * Ecuaciones en **LaTeX** del polinomio en su forma de Newton y en forma simplificada.
  * Indicador destacado (`st.metric`) con el valor de latencia interpolado $P_n(x_{\text{eval}})$.
  * **Gráfico Interactivo Matplotlib:**
    * Curva suave $P_n(x)$ trazada en un rango continuo.
    * Puntos experimentales marcados en rojo.
    * Punto interpolado $(x_{\text{eval}}, P_n(x_{	ext{eval}}))$ destacado en verde con líneas punteadas proyectadas a los ejes.

---

## 🔧 Actualización del Repositorio Git y Documentación

1. **Gestión de Archivos e Ignorados:**
   * Crear el archivo `.gitignore` incluyendo:
     ```text
     __pycache__/
     *.pyc
     .pytest_cache/
     .vscode/
     .DS_Store
     ```
   * Remover del seguimiento de Git los archivos en cache compilados: `git rm -r --cached __pycache__`.

2. **Actualización del `README.md`:**
   * Documentar la arquitectura general de la aplicación con los 4 módulos activos (Factorización LU Doolittle, Métodos Iterativos Jacobi/Gauss-Seidel, Interpolación de Lagrange e Interpolación de Newton).
   * Detallar las instrucciones de instalación y ejecución (`streamlit run app.py`).

3. **Normas de Commit y Push:**
   * Guardar cambios y realizar el commit con un mensaje estructurado.
   * **REGLA ESTRICTA DE GIT:** No agregar firmas ni trailers de coautoría (`Co-authored-by:`). Realizar el commit limpio y hacer `git push` a la rama `main`.

---

## 🚀 Prompt Directo para Claude Code

> *"Integra el módulo de la Sesión 9 (Interpolación Polinómica de Newton mediante Diferencias Divididas) en la aplicación Streamlit `D:\MN-U2`. Crea el archivo `newton_solver.py` con la lógica numérica pura y extiende `app.py` agregando la Sesión 9 en el menú de navegación de la barra lateral. Asegúrate de incluir la carga dinámica de puntos $(x_i, y_i)$, la tabla completa de diferencias divididas, las expresiones LaTeX del polinomio en forma de Newton y simplificada, la interpolación en $x_{\text{eval}}$, la gráfica continua con Matplotlib y la carga predeterminada del caso SRE REST API. Finalmente, crea el archivo `.gitignore` para desvincular `__pycache__`, actualiza profesionalmente el `README.md` documentando los 4 módulos y realiza commit y push sin ningún tráiler de coautoría."*
