# Prompt / Instrucciones para Claude Code: Aplicación Streamlit de Descomposición LU (Métodos Numéricos)

## 🎯 Objetivo General
Desarrollar una aplicación web interactiva en **Python** utilizando **Streamlit** que implemente y automatice la **Factorización LU mediante el algoritmo de Doolittle** para la resolución de Sistemas de Ecuaciones Lineales ($A \mathbf{x} = \mathbf{b}$). La aplicación debe permitir configurar dimensiones parametrizables ($n \ge 2$), resolver pasos intermedios (triangulación, sustitución hacia adelante y hacia atrás), y permitir reevaluar múltiples vectores de carga/tráfico $\mathbf{b}$ sin recalcular la matriz $LU$.

---

## 📐 Fundamento Teórico y Caso de Verificación

### 1. Algoritmo de Doolittle
Dada una matriz cuadrada $A$ de orden $n \times n$, se descompone en $A = LU$:
- **$L$ (Matriz Triangular Inferior):** Diagonal principal con unos ($l_{ii} = 1$).
  $$l_{ji} = \frac{a_{ji} - \sum_{k=1}^{i-1} l_{jk} u_{ki}}{u_{ii}} \quad (j > i)$$
- **$U$ (Matriz Triangular Superior):**
  $$u_{ij} = a_{ij} - \sum_{k=1}^{i-1} l_{ik} u_{kj} \quad (j \ge i)$$

### 2. Resolución en 2 Pasos
1. **Sustitución hacia adelante ($L\mathbf{y} = \mathbf{b}$):** Determina el vector intermedio $\mathbf{y}$.
2. **Sustitución hacia atrás ($U\mathbf{x} = \mathbf{y}$):** Determina el vector de solución $\mathbf{x}$.

### 3. Caso Test / Benchmark de Validación (Caso predeterminado de la guía)
* **Matriz de coeficientes $A$:**
  $$A = \begin{pmatrix} 4 & 2 & 1 \\ 12 & 10 & 5 \\ -8 & 8 & 7 \end{pmatrix}$$
* **Vector SLA inicial $\mathbf{b}_1$:** $[14, 46, 26]^T$
* **Resultados Esperados:**
  * **$L$:** $\begin{pmatrix} 1 & 0 & 0 \\ 3 & 1 & 0 \\ -2 & 3 & 1 \end{pmatrix}$
  * **$U$:** $\begin{pmatrix} 4 & 2 & 1 \\ 0 & 4 & 2 \\ 0 & 0 & 3 \end{pmatrix}$
  * **Vector Intermedio $Y$ ($L\mathbf{y} = \mathbf{b}_1$):** $a=14,\ b=4,\ c=42 \implies \mathbf{y} = [14, 4, 42]^T$
  * **Vector Solución $X$ ($U\mathbf{x} = \mathbf{y}$):** $x_1=3,\ x_2=-6,\ x_3=14 \implies \mathbf{x} = [3, -6, 14]^T$
* **Vector de tráfico alternativo $\mathbf{b}_2$ (Reevaluación en tiempo real sin recalcular $LU$):**
  * $\mathbf{b}_2 = [20, 62, 30]^T$

---

## 🛠️ Requisitos Técnicos y de Implementación

1. **Lenguaje y Librerías:**
   - Python 3.9+
   - `streamlit`
   - `numpy`
   - `pandas`

2. **Arquitectura del Código:**
   - Todo el cálculo numérico debe estar **100% automatizado en Python** (prohibido hardcodear resultados o solicitar cálculos manuales).
   - Crear una clase o módulo independiente `LUSolver`:
     - `doolittle_decompose(A)`: retorna $L$ y $U$.
     - `forward_substitution(L, b)`: retorna $\mathbf{y}$.
     - `backward_substitution(U, y)`: retorna $\mathbf{x}$.

3. **Flexibilidad e Interactividad:**
   - Selector numérico para la dimensión del sistema $N$ (mínimo $N=2$, configurable a $N=3, 4, \dots$).
   - Formulario/Edición matricial dinámica (vía `st.data_editor` o inputs numéricos individuales).
   - Botón de **"Cargar Ejemplo Predeterminado (Sesión 6)"** que llene automáticamente la matriz $A$ y el vector $\mathbf{b}_1$ del problema del clúster de microservicios.
   - Opción para ingresar un nuevo vector $\mathbf{b}_2$ y resolverlo reutilizando las matrices $L$ y $U$ ya calculadas.

4. **Validación y Robustez Numérica (Manejo de Pivote Nulo):**
   - El algoritmo de Doolittle puro requiere división por los pivotes $u_{ii}$. Si durante la descomposición se detecta un pivote nulo o casi nulo ($|u_{ii}| < 10^{-12}$):
     - Capturar la excepción numéricamente y mostrar una advertencia/error amigable en la interfaz de Streamlit (`st.error("División por cero detectada: $u_{ii} \approx 0$. La matriz requiere pivoteo o permutar filas.")`).
     - *(Opcional)* Implementar comprobación de determinante o Pivoteo Parcial ($PA = LU$) si se desea admitir cualquier matriz arbitraria sin colapsar la aplicación.

5. **Interfaz de Usuario (Streamlit UI):**
   - **Encabezado y Contexto:** Título claro, explicación breve del balanceo de carga en clústeres Cloud.
   - **Panel Lateral (Sidebar):** Configuración de dimensiones ($N$), presets y parámetros.
   - **Panel Principal:**
     - Tablas interactivas para $A$ y $\mathbf{b}$.
     - Despliegue de resultados organizados en columnas/tarjetas:
       - Matrices $L$ y $U$.
       - Proceso de triangulación.
       - Vector $Y$ (Sustitución hacia adelante).
       - Vector $X$ (Solución final).
     - Renderizado formal con LaTeX para fórmulas y expresiones matemáticas.

---

## 📋 Estructura de Archivos Sugerida
```text
project/
├── app.py
├── lu_solver.py
└── requirements.txt
```

---

## 🚀 Instrucción de Ejecución para Claude Code
> "Genera la aplicación en Streamlit implementando el algoritmo de Doolittle para la descomposición LU en Python de forma totalmente automatizada. Asegúrate de cumplir los resultados numéricos del caso predeterminado ($A, L, U, Y, X$) especificados en este documento, permitir sistemas dinámicos de $N \ge 2$ variables con reevaluación rápida de vectores $\mathbf{b}$, e incluir validación de pivotes nulos ($u_{ii} \approx 0$) para evitar errores en tiempo de ejecución."
