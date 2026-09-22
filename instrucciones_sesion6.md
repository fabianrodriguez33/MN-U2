# Prompt / Instrucciones para Claude Code: Aplicación Streamlit de Métodos Iterativos (Jacobi y Gauss-Seidel - Sesión 7)

## 🎯 Objetivo General
Extender la aplicación Streamlit existente (`MN-U2`) conservando el módulo de la **Sesión 6 (Factorización LU - Doolittle)** y agregando el módulo de la **Sesión 7: Métodos Iterativos (Jacobi y Gauss-Seidel)** para la resolución de Sistemas de Ecuaciones Lineales ($A \mathbf{x} = \mathbf{b}$). La aplicación debe permitir configurar dimensiones parametrizables ($N \ge 2$), ingresar matrices personalizadas y vectores iniciales $\mathbf{x}^{(0)}$, verificar automáticamente la dominancia diagonal y el criterio de Sassenfeld, desplegar tablas iterativas completas paso a paso, comparar la velocidad de convergencia entre Jacobi y Gauss-Seidel, y graficar la evolución del error.

---

## 🏗️ Estructura del Proyecto Actualizada

```text
MN-U2/
├── app.py                  # Interfaz Streamlit con navegador lateral (Sesión 6 vs Sesión 7)
├── lu_solver.py            # Módulo numérico para Factorización LU Doolittle (Sesión 6)
├── iterative_solvers.py    # NUEVO: Módulo numérico para Jacobi y Gauss-Seidel (Sesión 7)
├── instrucciones.md        # Especificaciones de la Sesión 6
├── instrucciones_sesion7.md# Especificaciones de la Sesión 7
├── requirements.txt        # Dependencias: streamlit, numpy, pandas, matplotlib
└── .claude/
    └── launch.json         # Configuración del servidor de desarrollo Streamlit
```

---

## 📐 Fundamento Teórico y Algoritmos (Sesión 7)

### 1. Descomposición Matricial Base
Dada una matriz $A \in \mathbb{R}^{n \times n}$, se descompone como $A = D - L - U$ (o $A = L + D + R$), donde $D$ es la diagonal principal, $L$ la triangular inferior estricta y $U$ la triangular superior estricta.

### 2. Algoritmo de Jacobi (Desplazamientos Simultáneos)
Utiliza únicamente los valores de la iteración anterior $k$ para calcular todos los componentes de la nueva iteración $k+1$:
$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j \neq i} a_{ij} x_j^{(k)} \right)$$

### 3. Algoritmo de Gauss-Seidel (Actualización Inmediata)
Emplea de forma inmediata los valores recién actualizados en la iteración corriente $k+1$ para los componentes anteriores ($j < i$):
$$x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)$$

### 4. Criterios de Convergencia
* **Estricta Dominancia Diagonal (EDD) por Filas:** Condición suficiente para garantizar convergencia en ambos métodos.
  $$|a_{ii}| > \sum_{j \neq i} |a_{ij}| \quad \forall i = 1, \dots, n$$
* **Criterio de Sassenfeld (Específico para Gauss-Seidel):**
  $$\beta_1 = \frac{1}{|a_{11}|} \sum_{j=2}^{n} |a_{1j}|, \quad \beta_i = \frac{\sum_{j=1}^{i-1} |a_{ij}| \beta_j + \sum_{j=i+1}^{n} |a_{ij}|}{|a_{ii}|}$$
  Si $\max_i \beta_i < 1$, se garantiza la convergencia de Gauss-Seidel.
* **Nota sobre no dominancia:** Si no existe EDD, la app debe advertir que la dominancia diagonal es condición suficiente pero no necesaria, procediendo a verificar Sassenfeld y la convergencia experimental.

### 5. Cálculo del Error y Residuo
* **Error Relativo Aproximado:**
  $$E = \max_i \left| \frac{x_i^{(k+1)} - x_i^{(k)}}{x_i^{(k+1)}} \right| \times 100\%$$
* **Vector de Residuo:** $\mathbf{r} = \mathbf{b} - A \mathbf{x}^{(k)}$

---

## 📊 Caso Benchmark / Preset de Validación (Guía Autónoma Sesión 7)

### Contexto: Balanceo de Carga en Clúster de 4 Servidores
Modelado en estado estacionario para $N=4$ nodos ($x_1, x_2, x_3, x_4$ en miles de peticiones/minuto).

### Sistema de Ecuaciones:
$$\begin{cases}
10x_1 - 2x_2 - x_3 = 15 \\
-x_1 + 8x_2 - 2x_4 = 18 \\
-2x_1 + 12x_3 - 3x_4 = 25 \\
-x_2 - 2x_3 + 9x_4 = 20
\end{cases}$$

$$\text{Matriz } A = \begin{pmatrix} 10 & -2 & -1 & 0 \\ -1 & 8 & 0 & -2 \\ -2 & 0 & 12 & -3 \\ 0 & -1 & -2 & 9 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 15 \\ 18 \\ 25 \\ 20 \end{pmatrix}, \quad \mathbf{x}^{(0)} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}, \quad \epsilon = 10^{-4} \text{ (0.01\%)}$$

### Verificación Teórica de Dominancia Diagonal (Tabla Integrada):
| Fila ($i$) | Elemento Diagonal $|a_{ii}|$ | Suma Fuera de Diagonal $\sum_{j \neq i} |a_{ij}|$ | Condición $|a_{ii}| > \sum |a_{ij}|$ | Resultado EDD |
| :---: | :---: | :---: | :---: | :---: |
| 1 | $|10| = 10$ | $|-2| + |-1| + |0| = 3$ | $10 > 3$ | **Cumple** |
| 2 | $|8| = 8$ | $|-1| + |0| + |-2| = 3$ | $8 > 3$ | **Cumple** |
| 3 | $|12| = 12$ | $|-2| + |0| + |-3| = 5$ | $12 > 5$ | **Cumple** |
| 4 | $|9| = 9$ | $|0| + |-1| + |-2| = 3$ | $9 > 3$ | **Cumple** |

---

## 🛠️ Requisitos Técnicos e Implementación en Python/Streamlit

### 1. Nuevo Módulo Numérico (`iterative_solvers.py`)
Crear la clase o funciones puras para métodos iterativos:
* `check_diagonal_dominance(A)`: Retorna booleano por fila y global.
* `check_sassenfeld(A)`: Retorna vector de coeficientes $\beta_i$ y booleano $\max \beta_i < 1$.
* `jacobi_solver(A, b, x0, tol, max_iter)`: Ejecuta iteraciones de Jacobi. Retorna historial en DataFrame, solución final $\mathbf{x}$, residuo $\mathbf{r}$, número de iteraciones y bandera de convergencia.
* `gauss_seidel_solver(A, b, x0, tol, max_iter)`: Ejecuta iteraciones de Gauss-Seidel. Retorna la misma estructura de resultados.

### 2. Actualización de Interfaz (`app.py`)
* **Navegación Múltiple (Sidebar):** Implementar un menú `st.sidebar.radio` para cambiar de módulo:
  1. **Sesión 6:** Factorización LU (Doolittle) — mantener intacto.
  2. **Sesión 7:** Métodos Iterativos (Jacobi vs. Gauss-Seidel).
* **Parámetros Personalizables para Sesión 7:**
  - Selector de dimensión $N \ge 2$.
  - Editor interactivo para la matriz $A$ y el vector $\mathbf{b}$ (`st.data_editor`).
  - Vector inicial $\mathbf{x}^{(0)}$ personalizable (por defecto ceros).
  - Tolerancia de error $\epsilon$ (ej. $10^{-4}$) e Iteraciones Máximas $K_{\max}$ (ej. 100).
  - Botón **"Cargar Ejemplo Predeterminado (Sesión 7 - Clúster de 4 Servidores)"**.
  - Opción para ingresar sistemas/ecuaciones personalizadas manualmente.
* **Despliegue de Resultados:**
  - **Sección de Verificación de Convergencia:** Muestra la tabla de Dominancia Diagonal y los Coeficientes de Sassenfeld $\beta_i$ con tarjetas `st.metric` o `st.success`/`st.warning`.
  - **Tabla Historial Iterativo:** Despliega paso a paso $k, x_1^{(k)}, \dots, x_N^{(k)}, \text{Error Relative (\%)}$.
  - **Tabla Comparativa de Métodos:**
    | Característica | Método de Jacobi | Método de Gauss-Seidel |
    | :--- | :--- | :--- |
    | **Uso de Datos** | Emplea valores de la iteración anterior ($k$) | Emplea valores actualizados inmediatamente ($k+1$) |
    | **Iteraciones Requeridas** | Generalmente más iteraciones | Suele converger en aprox. la mitad de iteraciones |
    | **Paralelización** | Alta (componentes independientes) | Secuencial (dependencia interna) |
  - **Gráfico de Convergencia:** Gráfico de líneas (Error vs. Iteración $k$) comparando Jacobi vs Gauss-Seidel en el mismo eje.
  - **Vector de Residuo $\mathbf{r}$:** Cálculo de $\mathbf{r} = \mathbf{b} - A\mathbf{x}$.

---

## 🚀 Instrucción de Ejecución para Claude Code
> "Integra el módulo de la Sesión 7 (Métodos Iterativos: Jacobi y Gauss-Seidel) en la aplicación Streamlit existente `MN-U2`. Crea el archivo `iterative_solvers.py` con las funciones numéricas puras y actualiza `app.py` con un menú de navegación en la barra lateral para alternar entre la Sesión 6 y la Sesión 7. Asegúrate de incluir la verificación teórica de dominancia diagonal, el criterio de Sassenfeld, el caso de prueba del clúster de 4 servidores de la guía autónoma, tablas iterativas completas, comparación gráfica de convergencia y soporte para sistemas de $N \ge 2$ variables con vectores iniciales personalizados."
