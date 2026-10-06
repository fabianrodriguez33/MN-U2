const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell,
  HeadingLevel, AlignmentType, WidthType, ShadingType, BorderStyle, LevelFormat,
  Footer, PageNumber,
} = require("docx");

const FONT = "Calibri";
const EV = "D:/MN-U2/evidencias_s9/";

const p = (text, opts = {}) =>
  new Paragraph({ spacing: { after: 120 }, ...opts, children: [new TextRun({ text, font: FONT, size: 22, ...(opts.run || {}) })] });
const rich = (runs, opts = {}) =>
  new Paragraph({ spacing: { after: 120 }, ...opts, children: runs.map(r => typeof r === "string" ? new TextRun({ text: r, font: FONT, size: 22 }) : new TextRun({ font: FONT, size: 22, ...r })) });
const eq = (text) =>
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 60, after: 160 },
    children: [new TextRun({ text, font: "Cambria Math", size: 23, italics: true })] });
const h1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text: t, font: FONT })] });
const h2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun({ text: t, font: FONT })] });
const bullet = (t) => new Paragraph({ numbering: { reference: "bul", level: 0 }, spacing: { after: 60 },
  children: [new TextRun({ text: t, font: FONT, size: 22 })] });
const code = (lines) => lines.map(l => new Paragraph({
  spacing: { after: 0 }, shading: { type: ShadingType.CLEAR, fill: "F2F2F2" },
  children: [new TextRun({ text: l, font: "Consolas", size: 19 })] }));

const border = { style: BorderStyle.SINGLE, size: 4, color: "BFBFBF" };
const borders = { top: border, bottom: border, left: border, right: border };
function table(widths, header, rows, boldRow = -1) {
  const total = widths.reduce((a, b) => a + b, 0);
  const cell = (t, w, head, bold) => new TableCell({
    width: { size: w, type: WidthType.DXA }, borders,
    shading: head ? { type: ShadingType.CLEAR, fill: "1F3864" } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({
      text: String(t), font: FONT, size: 20, bold: head || bold, color: head ? "FFFFFF" : "000000" })] })],
  });
  return new Table({
    width: { size: total, type: WidthType.DXA }, columnWidths: widths,
    rows: [
      new TableRow({ tableHeader: true, children: header.map((t, i) => cell(t, widths[i], true, true)) }),
      ...rows.map((r, ri) => new TableRow({ children: r.map((t, i) => cell(t, widths[i], false, ri === boldRow)) })),
    ],
  });
}

function figure(file, caption, w = 380) {
  const h = Math.round(w * 770 / 800);
  return [
    new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 }, keepNext: true,
      children: [new ImageRun({ type: "jpg", data: fs.readFileSync(EV + file),
        transformation: { width: w, height: h }, altText: { title: caption, description: caption, name: file } })] }),
    new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 },
      children: [new TextRun({ text: caption, font: FONT, size: 19, italics: true, color: "595959" })] }),
  ];
}

const children = [
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 600, after: 120 },
    children: [new TextRun({ text: "Informe de Guía Autónoma — Sesión 9", font: FONT, size: 40, bold: true, color: "1F3864" })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 80 },
    children: [new TextRun({ text: "Interpolación Polinómica de Newton mediante Diferencias Divididas", font: FONT, size: 28, bold: true })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 80 },
    children: [new TextRun({ text: "Latencia de una API REST en Cloud (Ingeniería SRE)", font: FONT, size: 24 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 60 },
    children: [new TextRun({ text: "Proyecto MN-U2 — Métodos Numéricos", font: FONT, size: 22, color: "595959" })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 300 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 8, color: "1F3864", space: 6 } },
    children: [new TextRun({ text: "Autor: Julio Fabián Rodríguez Bazán  |  6 de octubre de 2026", font: FONT, size: 20, color: "595959" })] }),

  h1("1. Objetivo"),
  p("Integrar en la aplicación Streamlit MN-U2 el módulo de la Sesión 9: interpolación polinómica de Newton mediante diferencias divididas. El módulo permite ingresar puntos (xᵢ, yᵢ) de forma dinámica, construir la tabla triangular de diferencias divididas, obtener el polinomio interpolante en forma de Newton y en forma canónica, estimar la latencia en un nivel de carga no medido y graficar el resultado. Las Sesiones 6, 7 y 8 se conservan sin cambios y se accede a las cuatro desde la barra lateral."),

  h1("2. Fundamento Teórico"),
  h2("2.1 Diferencias divididas"),
  p("Dados n+1 puntos con abscisas distintas, las diferencias divididas se definen recursivamente:"),
  eq("f[xᵢ] = yᵢ"),
  eq("f[xᵢ, xᵢ₊₁, …, xᵢ₊ₖ] = ( f[xᵢ₊₁, …, xᵢ₊ₖ] − f[xᵢ, …, xᵢ₊ₖ₋₁] ) / ( xᵢ₊ₖ − xᵢ )"),
  h2("2.2 Polinomio de Newton"),
  p("Con los coeficientes aₖ = f[x₀, …, xₖ] (primera fila de la tabla), el polinomio interpolante es:"),
  eq("Pₙ(x) = a₀ + a₁(x − x₀) + a₂(x − x₀)(x − x₁) + … + aₙ ∏ⱼ₌₀ⁿ⁻¹ (x − xⱼ)"),
  p("Su evaluación se realiza con el esquema anidado de Horner, que requiere solo n multiplicaciones:"),
  eq("Pₙ(x) = a₀ + (x − x₀)( a₁ + (x − x₁)( a₂ + … + (x − xₙ₋₁) aₙ ) )"),
  h2("2.3 Validación"),
  p("Si existen xᵢ = xⱼ con i ≠ j, el denominador (xᵢ₊ₖ − xᵢ) puede anularse. El módulo lanza la excepción personalizada DuplicateXError antes de calcular, evitando la división por cero."),

  h1("3. Caso Benchmark: API REST en Cloud"),
  p("Un equipo SRE mide la latencia media de respuesta y (ms) de una API REST crítica según la carga concurrente x, expresada en unidades de 100 req/s. Se dispone de 4 mediciones (n = 3):"),
  table([1500, 3600, 3600], ["Punto i", "Carga xᵢ (×100 req/s)", "Latencia yᵢ (ms)"],
    [["0", "1", "45"], ["1", "2", "65"], ["2", "4", "110"], ["3", "7", "220"]]),
  p("", { spacing: { after: 60 } }),
  p("Se desea estimar la latencia para una carga de 500 req/s (x_eval = 5), valor no medido."),

  h1("4. Procedimiento Numérico"),
  h2("4.1 Tabla de diferencias divididas (cálculo manual)"),
  bullet("1.º orden: (65 − 45)/(2 − 1) = 20;  (110 − 65)/(4 − 2) = 22.5;  (220 − 110)/(7 − 4) = 36.6667"),
  bullet("2.º orden: (22.5 − 20)/(4 − 1) = 0.8333;  (36.6667 − 22.5)/(7 − 2) = 2.8333"),
  bullet("3.er orden: (2.8333 − 0.8333)/(7 − 1) = 0.3333"),
  p("", { spacing: { after: 60 } }),
  table([700, 800, 1500, 1800, 1800, 1800],
    ["i", "xᵢ", "f[xᵢ]", "1.ª dif.", "2.ª dif.", "3.ª dif."],
    [["0", "1", "45.0000", "20.0000", "0.8333", "0.3333"],
     ["1", "2", "65.0000", "22.5000", "2.8333", "--"],
     ["2", "4", "110.0000", "36.6667", "--", "--"],
     ["3", "7", "220.0000", "--", "--", "--"]], 0),
  p("", { spacing: { after: 60 } }),
  p("Coeficientes de Newton (fila 0): a₀ = 45, a₁ = 20, a₂ = 5/6 ≈ 0.833333, a₃ = 1/3 ≈ 0.333333."),
  h2("4.2 Polinomio interpolante"),
  p("Forma de Newton:"),
  eq("P₃(x) = 45 + 20(x − 1) + (5/6)(x − 1)(x − 2) + (1/3)(x − 1)(x − 2)(x − 4)"),
  p("Forma canónica, al expandir y agrupar términos (x³: 1/3; x²: 5/6 − 7/3 = −3/2; x: 20 − 5/2 + 14/3 = 133/6; constante: 45 − 20 + 5/3 − 8/3 = 24):"),
  eq("P₃(x) = (1/3)x³ − (3/2)x² + (133/6)x + 24 ≈ 0.333333x³ − 1.5x² + 22.1667x + 24"),
  h2("4.3 Interpolación en x_eval = 5"),
  eq("P₃(5) = 45 + 20(4) + (5/6)(4)(3) + (1/3)(4)(3)(1) = 45 + 80 + 10 + 4 = 139 ms"),
  p("La latencia estimada para 500 req/s es de 139 ms."),

  h1("5. Implementación"),
  h2("5.1 Módulo newton_solver.py"),
  p("Contiene la lógica numérica pura en la clase NewtonSolver:"),
  bullet("compute_divided_differences(x, y): devuelve la matriz n × n de diferencias divididas (triangular superior; el resto, NaN). Valida duplicados y lanza DuplicateXError."),
  bullet("coefficients(table): extrae la primera fila, es decir, los coeficientes aₖ."),
  bullet("evaluate(x_eval, x_nodes, coefs): evaluación por Horner; acepta un escalar o un arreglo NumPy (usado para trazar la curva)."),
  bullet("get_symbolic_expressions(x_nodes, coefs): con sympy devuelve la expresión de Newton factorizada y la forma canónica expandida, con coeficientes racionales exactos."),
  bullet("newton_latex y canonical_latex: generan las ecuaciones LaTeX mostradas en la interfaz."),
  h2("5.2 Interfaz en app.py"),
  bullet("Barra lateral (st.sidebar.radio) con las cuatro sesiones: 6, 7, 8 y 9."),
  bullet("Selector numérico de cantidad de puntos (n+1 ≥ 2) y st.data_editor para editar los pares (xᵢ, yᵢ)."),
  bullet("Botón «Cargar Ejemplo Predeterminado (Sesión 9 - SRE REST API)» que autocompleta el benchmark."),
  bullet("Resultados: tabla de diferencias divididas (pandas), polinomio en LaTeX (Newton y canónico), st.metric con Pₙ(x_eval), tabla de verificación en los nodos y gráfico Matplotlib."),
  bullet("Aviso de extrapolación si x_eval cae fuera del intervalo [mín xᵢ, máx xᵢ]."),
  h2("5.3 Repositorio"),
  bullet("Se añadió .gitignore (__pycache__/, *.pyc, .pytest_cache/, .vscode/, .DS_Store) y se retiró __pycache__ del seguimiento de Git."),
  bullet("Se actualizó el README.md con la documentación de los cuatro módulos, la instalación y la ejecución (streamlit run app.py)."),

  h1("6. Evidencias de Ejecución"),
  p("Las siguientes capturas corresponden a la aplicación en ejecución local (streamlit run app.py, http://localhost:8501) con el caso predeterminado cargado."),
  h2("6.1 Navegación multi-módulo y carga de datos"),
  p("La barra lateral permite conmutar entre las Sesiones 6, 7, 8 y 9. Se muestran el botón del ejemplo predeterminado, el número de puntos (4) y x_eval = 5."),
  ...figure("01_navegacion_y_datos.jpg", "Figura 1. Navegación lateral, configuración de la Sesión 9 y tabla de puntos editable."),
  h2("6.2 Puntos experimentales"),
  ...figure("02_encabezado_puntos.jpg", "Figura 2. Descripción del método y puntos (xᵢ, yᵢ) del caso SRE REST API en st.data_editor."),
  h2("6.3 Tabla de diferencias divididas y polinomios"),
  p("La tabla reproduce los valores calculados manualmente (20, 22.5, 36.6667; 0.8333, 2.8333; 0.3333). Los coeficientes de Newton y la forma canónica P₃(x) = 0.333333x³ − 1.5x² + 22.1667x + 24 coinciden con el cálculo de la sección 4. La línea de la forma de Newton aparece recortada por el ancho de la ventana de captura; la expresión completa es la de la sección 4.2."),
  ...figure("03_tabla_dd_polinomios.jpg", "Figura 3. Tabla de diferencias divididas, coeficientes de Newton y polinomio interpolante."),
  h2("6.4 Valor interpolado y verificación en los nodos"),
  p("El indicador st.metric muestra P₃(5) = 139.000000 ms, igual al resultado manual. La verificación en los cuatro nodos arroja error 0.0000 en todos los casos, como corresponde a un polinomio interpolante."),
  ...figure("04_valor_interpolado_verificacion.jpg", "Figura 4. Latencia interpolada P₃(5) = 139 ms y tabla de verificación |yᵢ − P₃(xᵢ)|."),
  h2("6.5 Gráfico de la interpolación"),
  p("Curva continua de P₃(x), puntos experimentales en rojo y punto interpolado (5, 139) en verde con líneas punteadas proyectadas a los ejes."),
  ...figure("05_grafico.jpg", "Figura 5. Gráfico Matplotlib de la interpolación polinómica de Newton."),
  h2("6.6 Validación de nodos repetidos"),
  p("Al ingresar x = [1, 2, 2, 7], el módulo detiene el cálculo y lanza DuplicateXError, evitando la división por cero:"),
  ...code([
    ">>> NewtonSolver.compute_divided_differences([1, 2, 2, 7], [45, 65, 70, 220])",
    "DuplicateXError: Error: el valor x = 2 está repetido. Los nodos x_i deben ser",
    "distintos para evitar división por cero en L_k(x).",
  ]),
  p("", { spacing: { after: 120 } }),
  p("En la interfaz, esta excepción se captura y se muestra como un mensaje de error (st.error) sin detener la aplicación."),

  h1("7. Comparación de Resultados"),
  table([3200, 2000, 2000, 2000],
    ["Magnitud", "Cálculo manual", "Aplicación", "Diferencia"],
    [["a₀", "45", "45.000000", "0"],
     ["a₁", "20", "20.000000", "0"],
     ["a₂", "0.833333", "0.833333", "0"],
     ["a₃", "0.333333", "0.333333", "0"],
     ["P₃(5) (ms)", "139", "139.000000", "0"],
     ["Error máx. en nodos", "0", "0.0000", "0"]]),
  p("", { spacing: { after: 120 } }),

  h1("8. Conclusiones"),
  bullet("El módulo reproduce exactamente el benchmark: a₀…a₃ = 45, 20, 5/6, 1/3 y P₃(5) = 139 ms."),
  bullet("La forma de Newton permite agregar un nuevo punto sin recalcular los coeficientes anteriores, ventaja frente a Lagrange, que obliga a reconstruir todos los polinomios base."),
  bullet("El esquema de Horner evalúa Pₙ(x) con solo n multiplicaciones y se usa tanto para x_eval como para trazar la curva."),
  bullet("Con nodos distintos el polinomio pasa por todos los puntos medidos (error cero en la verificación). Fuera del intervalo [1, 7] se trata de extrapolación y la interfaz lo advierte."),
  bullet("Un polinomio de grado alto puede oscilar entre nodos (fenómeno de Runge); para muchos puntos conviene interpolación por tramos."),
  bullet("La aplicación integra ahora cuatro módulos (LU, Jacobi/Gauss-Seidel, Lagrange y Newton) con navegación lateral y sin alterar las sesiones previas."),
];

const doc = new Document({
  creator: "Julio Fabián Rodríguez Bazán",
  title: "Informe Sesión 9 — Interpolación de Newton",
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 30, bold: true, font: FONT, color: "1F3864" },
        paragraph: { spacing: { before: 320, after: 140 }, outlineLevel: 0, keepNext: true } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 25, bold: true, font: FONT, color: "2E5597" },
        paragraph: { spacing: { before: 220, after: 100 }, outlineLevel: 1, keepNext: true } },
    ],
  },
  numbering: { config: [{ reference: "bul", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•",
    alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1300, right: 1300, bottom: 1300, left: 1300 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [new TextRun({ text: "MN-U2 — Sesión 9 — Página ", font: FONT, size: 18, color: "7F7F7F" }),
        new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18, color: "7F7F7F" })] })] }) },
    children,
  }],
});

Packer.toBuffer(doc).then(b => { fs.writeFileSync("D:/MN-U2/Informe_Sesion9_Newton_Diferencias_Divididas.docx", b); console.log("ok"); });

