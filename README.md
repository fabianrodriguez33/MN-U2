# MN-U2 — Métodos Numéricos aplicados a Cloud / SRE

Aplicación web en **Streamlit** con cuatro módulos de Métodos Numéricos: resolución de sistemas lineales $A\mathbf{x} = \mathbf{b}$ (métodos directos e iterativos) e interpolación polinómica (Lagrange y Newton). Cada módulo se aplica a un escenario de infraestructura Cloud. La navegación entre sesiones se hace desde la barra lateral.

## Módulos

| Sesión | Módulo | Archivo | Descripción |
| :---: | --- | --- | --- |
| 6 | Factorización LU (Doolittle) | `lu_solver.py` | Descompone $A = LU$, resuelve por sustitución hacia adelante y hacia atrás, y reutiliza $L$ y $U$ con nuevos vectores $\mathbf{b}$. Valida pivotes nulos ($\lvert u_{ii}\rvert < 10^{-12}$). |
| 7 | Métodos iterativos (Jacobi y Gauss-Seidel) | `iterative_solvers.py` | Verifica dominancia diagonal estricta y criterio de Sassenfeld; muestra la tabla de iteraciones, el gráfico de convergencia del error y el residuo. |
| 8 | Interpolación de Lagrange | `lagrange_solver.py` | Ingreso dinámico de puntos $(x_i, y_i)$, tabla de $L_k(x_{eval})$, polinomio $P_n(x)$ simbólico (sympy), valor interpolado y gráfico. |
| 9 | Interpolación de Newton (diferencias divididas) | `newton_solver.py` | Tabla triangular de diferencias divididas, polinomio en forma de Newton y en forma canónica (LaTeX), evaluación con esquema de Horner, verificación en los nodos y gráfico. |

Las sesiones 8 y 9 lanzan `DuplicateXError` si existen valores $x_i$ repetidos.

## Estructura

```text
MN-U2/
├── app.py                  # Interfaz Streamlit (navegación Sesiones 6 / 7 / 8 / 9)
├── lu_solver.py            # Sesión 6: Factorización LU (Doolittle)
├── iterative_solvers.py    # Sesión 7: Jacobi, Gauss-Seidel, EDD, Sassenfeld
├── lagrange_solver.py      # Sesión 8: Interpolación de Lagrange
├── newton_solver.py        # Sesión 9: Diferencias divididas de Newton
├── instruccioness8.md      # Especificación Sesión 8
├── inss9.md                # Especificación Sesión 9
├── requirements.txt
├── .gitignore
└── .claude/launch.json     # Configuración del servidor de desarrollo (puerto 8501)
```

## Instalación

```bash
pip install -r requirements.txt
```

Dependencias: `streamlit`, `numpy`, `pandas`, `matplotlib`, `sympy`.

## Ejecución

```bash
streamlit run app.py
```

La aplicación se abre en `http://localhost:8501`.

## Casos de prueba predeterminados

Cada sesión incluye un botón en la barra lateral para cargar su ejemplo.

- **Sesión 6:** matriz $A$ de $3\times3$ con $\mathbf{b}_1 = [14, 46, 26]^T$ y $\mathbf{b}_2 = [20, 62, 30]^T$.
- **Sesión 7:** clúster de 4 servidores, $\mathbf{x}^{(0)} = \mathbf{0}$, $\epsilon = 10^{-4}$.
- **Sesión 8:** RAM $x = 2, 4, 8, 12$ GB, latencia $y = 150, 85, 50, 70$ ms, $x_{eval}=6$ → $P_3(6) = 55.25$ ms.
- **Sesión 9:** carga de una API REST $x = 1, 2, 4, 7$ (×100 req/s), latencia $y = 45, 65, 110, 220$ ms, $x_{eval}=5$ → $P_3(5) = 139$ ms, con
  $P_3(x) = 45 + 20(x-1) + \tfrac{5}{6}(x-1)(x-2) + \tfrac{1}{3}(x-1)(x-2)(x-4)$.
