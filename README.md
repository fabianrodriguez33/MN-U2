# MN-U2 — Métodos Numéricos: Sistemas de Ecuaciones Lineales

Aplicación Streamlit que resuelve sistemas de ecuaciones lineales $A\mathbf{x} = \mathbf{b}$ mediante métodos directos e iterativos, aplicados al balanceo de carga en clústeres Cloud.

## Módulos

- **Sesión 6 — Factorización LU (Doolittle):** descompone $A = LU$, resuelve por sustitución hacia adelante y hacia atrás, y permite reevaluar nuevos vectores $\mathbf{b}$ sin recalcular $L$ y $U$. Valida pivotes nulos ($|u_{ii}| < 10^{-12}$).
- **Sesión 7 — Métodos Iterativos (Jacobi y Gauss-Seidel):** verifica dominancia diagonal estricta y el criterio de Sassenfeld, ejecuta ambos métodos con tabla de iteraciones, gráfico de convergencia del error y cálculo del residuo.

## Estructura

```text
MN-U2/
├── app.py                     # Interfaz Streamlit (navegación Sesión 6 / Sesión 7)
├── lu_solver.py                # Factorización LU (Doolittle)
├── iterative_solvers.py        # Jacobi, Gauss-Seidel, EDD, Sassenfeld
├── instrucciones.md            # Especificación Sesión 6
├── instrucciones_sesion6.md    # Especificación Sesión 7
├── requirements.txt
└── .claude/launch.json         # Configuración del servidor de desarrollo
```

## Instalación

```bash
pip install -r requirements.txt
```

## Ejecución

```bash
streamlit run app.py
```

La aplicación se abre en `http://localhost:8501`.

## Casos de prueba predeterminados

- **Sesión 6:** matriz $A$ de $3\times3$ con $\mathbf{b}_1 = [14, 46, 26]^T$ y $\mathbf{b}_2 = [20, 62, 30]^T$ (botón "Cargar Ejemplo Predeterminado").
- **Sesión 7:** clúster de 4 servidores, $\mathbf{x}^{(0)} = \mathbf{0}$, $\epsilon = 10^{-4}$ (botón "Cargar Ejemplo Predeterminado - Clúster de 4 Servidores").
