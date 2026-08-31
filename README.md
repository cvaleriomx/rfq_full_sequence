# RFQ pipeline: de TRANSOPTR a Warp

Repositorio reproducible para organizar la creación, validación y simulación de
una RFQ. La fuente de verdad del diseño ISAC2 es
[`data/designs/isac2/table1.txt`](data/designs/isac2/table1.txt), una tabla
exportada por TRANSOPTR. A partir de ella se derivan los coeficientes del modelo
de dos términos, los perfiles de las vanes, un mapa de campo externo y la
configuración de tracking para Warp.

Este proyecto separa deliberadamente tres cosas que antes estaban mezcladas:

1. **Entradas canónicas:** la tabla de celdas y los archivos TOML.
2. **Código reproducible:** transformaciones, validaciones y simulación.
3. **Resultados derivados:** HDF5, CSV, figuras y partículas, siempre dentro de
   `outputs/` y fuera de Git.

La tabla `table2_dans.txt` se conserva únicamente en
`data/experiments/quick_tests/`. No participa en el flujo ISAC2 salvo que se
cree de forma explícita otra configuración que la seleccione.

## Diagrama correcto del flujo

En el diagrama, `table1.txt` debe aparecer **después de TRANSOPTR**, porque es
un producto de ese programa y el contrato de entrada del resto del repositorio.
Desde esa tabla nacen dos ramas independientes: el campo analítico externo y la
geometría/solución electrostática de Warp.

```mermaid
flowchart LR
    A[Requisitos del haz y RF] --> B[Diseño y tracking en TRANSOPTR]
    B --> C[Exportar tabla canónica<br/>table1.txt]
    C --> D{Validar esquema,<br/>unidades y longitud}
    D -->|válida| E[Derivar A01, A10, k y fase]
    D -->|válida| F[Crear perfiles de las 4 vanes]
    E --> G[Mapa 3D externo<br/>modelo de dos términos]
    F --> H[Geometría 3D de vanes en Warp]
    H --> I[Aplicar +V/2 y -V/2<br/>y resolver campo DC]
    G --> J[Comparar Ez sobre el eje]
    I --> J
    J --> K{¿Error aceptable?}
    K -->|no| B
    K -->|sí| L[Aplicar dependencia RF<br/>cos(2πft + φ)]
    L --> M[Tracking de partículas en Warp]
    M --> N[Transmisión, energía,<br/>emitancia y pérdidas]
    N --> O{¿Cumple objetivos?}
    O -->|no| B
    O -->|sí| P[Diseño validado]
```

### Estado de las etapas

| Etapa | Estado en este repositorio |
|---|---|
| Ejecutar el diseño dentro de TRANSOPTR | Externo; su salida se versiona como `table1.txt` |
| Validar la tabla y convertir a SI | Implementado |
| Derivar coeficientes de dos términos | Implementado |
| Generar perfiles de vanes | Implementado |
| Construir el mapa externo 3D | Implementado |
| Validar simetrías y graficar el campo | Implementado |
| Importar un eje calculado por Warp DC y comparar | Implementado |
| Construcción/solución DC completa de las vanes en Warp | Interfaz pendiente de consolidar; se conserva como una frontera explícita |
| Tracking en el mapa RF externo con Warp | Implementado; debe validarse primero con `--dry-run` |
| Ejemplo de modelo subrogado de Coulomb | Incluido en `examples/surrogate_coulomb/` |

La etapa DC no se marca artificialmente como terminada: primero hay que
consolidar y verificar la construcción de conductores de Warp, especialmente
la orientación de cada vane, la resolución de malla y las condiciones de
frontera. El comando `compare-fields` ya fija el formato que deberá exportar
esa etapa.

## Resultado verificado para `table1.txt`

La validación incluida comprueba los siguientes valores:

| Propiedad | Valor |
|---|---:|
| Filas totales | 88 |
| Celdas físicas | 87 |
| Última etiqueta | `87T` |
| Longitud final | 2.0582 m |
| Voltaje de la tabla | 41.53 kV |
| Energía síncrona inicial | 0.035 MeV |
| Energía síncrona final | 1.0167 MeV |
| SHA-256 de la tabla canónica | `024c2e5212decf899900d7ae50a29a12aabc223a795a9b18d1fea413d5a96adf` |

El hash permite detectar que alguien sustituyó la tabla aunque conserve el
mismo nombre.

## Estructura

```text
rfq_pipeline_github/
├── configs/
│   ├── isac2_quick.toml          # malla pequeña para desarrollo
│   └── isac2_production.toml     # malla fina y 100 000 partículas
├── data/
│   ├── designs/isac2/table1.txt  # entrada canónica
│   └── experiments/quick_tests/  # entradas no canónicas
├── docs/                         # contratos y notas de validación
├── examples/surrogate_coulomb/   # prueba ML independiente
├── outputs/                      # resultados regenerables, ignorados por Git
├── src/rfq_pipeline/
│   ├── design.py                 # lectura, validación y procedencia
│   ├── coefficients.py           # modelo KT de dos términos
│   ├── vanes.py                  # perfiles geométricos
│   ├── fieldmap.py               # HDF5 3D normalizado
│   ├── validation.py             # pruebas físicas básicas
│   ├── diagnostics.py            # gráficas y comparación de campos
│   ├── tracking.py               # adaptador para Warp
│   └── cli.py                    # interfaz de línea de comandos
├── tests/
├── pyproject.toml
└── environment.yml
```

## Instalación en el ambiente `warp_2`

Warp ya está instalado en ese ambiente, por lo que no se intenta descargarlo
ni reemplazarlo. Desde esta carpeta:

```bash
conda activate warp_2
python -m pip install -e .
```

El modo editable hace que los cambios en `src/` se usen inmediatamente. Para
comprobar únicamente el núcleo sin Warp:

```bash
python -m unittest discover -s tests -v
```

También se puede crear un ambiente separado para generación de tablas, mapas y
ML:

```bash
conda env create -f environment.yml
conda activate rfq_pipeline_core
```

Ese ambiente **no incluye Warp**; sirve para las etapas independientes del
simulador.

## Primera ejecución recomendada

La configuración rápida genera un mapa pequeño y permite comprobar el flujo en
segundos:

```bash
python -m rfq_pipeline --config configs/isac2_quick.toml all
```

Equivalente con `make`:

```bash
make quick
```

Los productos aparecen en `outputs/isac2_quick/`. Deben existir, entre otros:

| Archivo | Contenido |
|---|---|
| `design_validation.json` | dimensiones, energías y resultado de validación |
| `design_provenance.json` | ruta y hash de la tabla que produjo los resultados |
| `design_si.csv` | tabla normalizada a unidades SI |
| `coefficients.csv` | A01, A10, k y fase por celda |
| `fort.75` | formato legado derivado; nunca es la fuente de verdad |
| `vane_profiles.csv` | coordenadas de las cuatro puntas |
| `vane_profiles.png` | inspección geométrica |
| `fieldmap.h5` | potencial y campo 3D comprimidos |
| `field_validation.json` | errores relativos de simetría |
| `field_diagnostics.png` | campo sobre el eje y plano `yz` |
| `pipeline_summary.json` | resumen de toda la ejecución |

Después de validar la ejecución rápida se puede construir el mapa fino:

```bash
python -m rfq_pipeline --config configs/isac2_production.toml all
```

La malla de producción solicitada es aproximadamente `101 × 101 × 4117`. El
HDF5 almacena los ejes una sola vez y comprime los campos por bloques; no crea
un `DataFrame` con una fila por punto ni duplica un `pickle` de varios GB.

## Comandos por etapa

```bash
# 1. Tabla y procedencia
python -m rfq_pipeline --config configs/isac2_quick.toml validate-design

# 2. Coeficientes
python -m rfq_pipeline --config configs/isac2_quick.toml coefficients

# 3. Perfiles geométricos
python -m rfq_pipeline --config configs/isac2_quick.toml vanes

# 4. Campo externo 3D
python -m rfq_pipeline --config configs/isac2_quick.toml fieldmap

# 5. Validaciones y gráficas
python -m rfq_pipeline --config configs/isac2_quick.toml validate-field
python -m rfq_pipeline --config configs/isac2_quick.toml plot-field
```

Ejecutar por etapas facilita identificar si un cambio provino de TRANSOPTR, de
los coeficientes, de la geometría o de Warp.

## Qué significa cada parte de la tabla

El lector conserva todas las columnas, pero el núcleo utiliza principalmente:

| Columna | Interpretación | Unidad de entrada |
|---|---|---|
| `Cell` | número/etiqueta de celda; `0` es condición inicial | — |
| `V` | voltaje de diseño | kV |
| `Wsyn` | energía síncrona | MeV |
| `Phi` | fase síncrona | grados |
| `a` | apertura mínima | cm |
| `m` | modulación | adimensional |
| `L` | longitud de celda | cm |
| `Z` | extremo acumulado de la celda | cm |
| `A10` | coeficiente incluido por TRANSOPTR | adimensional |

`design_si.csv` cambia longitudes a metros y mantiene nombres con la unidad en
el encabezado. La fila `Cell=0` se usa como condición inicial; no se cuenta
como celda física ni produce una oscilación adicional de las vanes.

## Modelo físico del campo externo

Para cada celda física se usa

```text
k = π/L
A10 = (m² - 1) / [m² I0(k a) + I0(k m a)]
A01 = [1 - A10 I0(k a)] / a²
```

El potencial normalizado por 1 V de diferencia intervane es

```text
Φ(x,y,z) = 1/2 [A01(z)(x²-y²) + A10(z) I0(k(z)r) cos(ψ(z))]
E = -∇Φ
```

La fase suma exactamente `π` radianes por celda. Los coeficientes se
interpolan sobre la malla longitudinal regular. El archivo conserva tanto el
`A10` derivado de la geometría como `a10_transoptr`, para que la diferencia sea
visible y pueda estudiarse; no se reemplaza silenciosamente uno por otro.

### Convención de voltaje

`V_intervane` es la diferencia entre una familia de vanes y la familia
opuesta. Para una solución DC coherente se debe aplicar:

```text
vanes x: +V_intervane/2
vanes y: -V_intervane/2
```

Así, para `V_intervane = 41.53 kV`, cada familia está a `±20.765 kV`. Aplicar
`+41.53 kV` y `-41.53 kV` produciría una diferencia de `83.06 kV` y duplicaría
el campo.

El mapa HDF5 guarda `E_normalized` en `V/m por V_intervane`. En tracking se
reconstruye el campo temporal con

```text
E(x,y,z,t) = V_intervane cos(2π f t + φ) E_normalized(x,y,z)
```

La configuración actual usa `f = 162 MHz` y `V_intervane = 41.53 kV`.

## Formato de `fieldmap.h5`

```text
/axes/x_m
/axes/y_m
/axes/z_m
/fields/phi_per_v
/fields/ex_vpm_per_v
/fields/ey_vpm_per_v
/fields/ez_vpm_per_v
```

Los campos tienen orden `(nx, ny, nz)`. Los atributos del archivo registran
las unidades, el orden de ejes y la convención de voltaje. La validación prueba
que no haya `NaN` o infinitos y verifica las simetrías esperadas:

- `Ex` impar respecto a `x`;
- `Ey` impar respecto a `y`;
- `Ez` par respecto a `x` y a `y`.

Estas pruebas detectan varios errores de forma, orden de ejes y signo, pero no
sustituyen la comparación contra la solución 3D de las vanes.

## Comparación con el campo DC de Warp

La solución DC debe exportar un CSV sobre el eje con dos columnas:

```csv
z_m,ez_vpm_per_v
0.0000,0.0
0.0005,12.3
```

El campo tiene que estar dividido por el **mismo** `V_intervane` usado en Warp.
Entonces se ejecuta:

```bash
python -m rfq_pipeline \
  --config configs/isac2_production.toml \
  compare-fields ruta/al/warp_dc_axis.csv
```

Se generan `field_comparison.png` y `field_comparison.json`. La figura compara
la forma normalizada y el JSON reporta RMSE, correlación y amplitudes antes de
normalizar. Conviene inspeccionar ambas cosas: una correlación alta puede
ocultar un error de factor dos en el voltaje.

## Tracking con Warp

Antes de reservar memoria o mover partículas, revisar el plan:

```bash
python -m rfq_pipeline \
  --config configs/isac2_production.toml \
  track --dry-run
```

El reporte incluye tamaño del campo, velocidad inicial, `dt`, número de pasos,
ciclos RF y memoria aproximada para los tres componentes. Si es coherente:

```bash
python -m rfq_pipeline \
  --config configs/isac2_production.toml \
  track
```

La configuración de producción usa 100 000 macropartículas. El número de
partículas no cambia el mapa externo: cambia el muestreo del haz, la memoria de
partículas y, si se activa, el ruido/costo de la carga espacial. Se recomienda
esta progresión antes de una corrida final:

```text
1 000 → 10 000 → 30 000 → 100 000 partículas
```

Comparar en cada nivel transmisión, energía final, emitancia y pérdidas. Usar
100 000 sólo cuando el cambio de 30 000 a 100 000 sea menor que la tolerancia
física definida para el estudio.

Advertencias del adaptador actual:

- inyecta un solo haz uniforme en disco, no una distribución casada completa;
- `space_charge = false` por defecto para aislar primero el campo externo;
- no instala aún la superficie modulada de las vanes como *scraper*;
- la configuración rápida tiene una malla longitudinal deliberadamente gruesa
  y no debe usarse para resultados físicos finales.

Por ello, una corrida exitosa es una prueba de integración, no una validación
final de la RFQ.

## Configuración reproducible

No edite números dentro de los módulos Python. Copie uno de los TOML y cambie
allí los parámetros. Por ejemplo:

```bash
cp configs/isac2_production.toml configs/mi_estudio.toml
```

Las secciones son:

- `[project]`: nombre y directorio de salida;
- `[design]`: tabla canónica y valores esperados;
- `[vanes]`: resolución longitudinal del perfil;
- `[field]`: límites, pasos, bloques HDF5 y tolerancia de simetría;
- `[tracking]`: partículas, semilla, haz, RF, paso y malla de Warp.

Cada estudio debe usar un `output_dir` diferente para evitar confundir
resultados de mallas o voltajes distintos.

## Modelo subrogado de Coulomb

El ejemplo previo de ML está preservado como un proyecto autocontenido en
[`examples/surrogate_coulomb/`](examples/surrogate_coulomb/README.md). Incluye:

- esfera uniforme, gaussiana 3D y proyección espacial tipo KV;
- estudio de convergencia con el número de macropartículas;
- generación del dataset;
- entrenamiento y gráfica de paridad;
- predicción y gráfica de `Ey` sobre el eje `y` usando sólo el subrogado.

En ese ejemplo, una fila del **dataset** es una configuración física completa,
no una partícula. Las siete entradas describen geometría, carga y tres tamaños;
las diez salidas son nueve valores de campo en sondas y el potencial central.
El número de macropartículas se usa internamente para calcular cada etiqueta.

El subrogado de Coulomb y el flujo RFQ están separados porque resuelven
problemas diferentes. Para integrarlos en una fase posterior, el contrato
correcto será: estado compacto del haz/elemento → campo o métricas predichas,
siempre acompañado por un conjunto de prueba independiente y comparación con
Warp.

## Pruebas y criterio mínimo antes de aceptar un cambio

```bash
python -m unittest discover -s tests -v
python -m rfq_pipeline --config configs/isac2_quick.toml all
python -m rfq_pipeline --config configs/isac2_quick.toml track --dry-run
```

Un cambio no debe aceptarse si ocurre cualquiera de estos casos:

- cambia la tabla canónica sin actualizar la procedencia y justificarlo;
- falla el número de celdas o la longitud esperada;
- aparecen valores no finitos;
- se rompen las simetrías del campo;
- cambia la convención de voltaje sin convertir los datos antiguos;
- se versiona un resultado pesado que puede regenerarse.

GitHub Actions ejecuta las pruebas y el flujo rápido en cada `push` y *pull
request*. Warp no se ejecuta en CI porque requiere el ambiente especializado.

## Preparación para GitHub

Los archivos grandes ya están excluidos mediante `.gitignore`. Antes del primer
commit:

```bash
git status
find . -type f -size +95M -print
git add .
git status
git commit -m "Initial reproducible RFQ pipeline"
```

Después de crear un repositorio vacío en GitHub:

```bash
git remote add origin git@github.com:USUARIO/NOMBRE.git
git push -u origin main
```

GitHub bloquea archivos mayores de 100 MiB en Git normal; los mapas HDF5,
partículas, modelos y `pickle` deben permanecer en `outputs/`, en almacenamiento
de artefactos o, si realmente deben versionarse, en Git LFS. Véase la
[documentación oficial de archivos grandes de GitHub](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).

Antes de publicar también falta elegir una licencia compatible con el código y
con Warp. No se añadió una licencia automáticamente porque esa decisión depende
del propietario y de las condiciones de las dependencias.

## Problemas frecuentes

### `No module named rfq_pipeline`

Ejecute `python -m pip install -e .` desde la raíz del repositorio o use
temporalmente `PYTHONPATH=src`.

### `Warp no está disponible`

Active `warp_2` y confirme:

```bash
python -c "import warp; print(warp.__file__)"
```

### El campo tiene el doble de amplitud

Revise la convención `±V_intervane/2`. No use `±V_intervane` en las vanes.

### La fase RF termina antes que el tracking

No construya una señal con un número fijo de muestras. El adaptador calcula
`steps + 2` tiempos usando exactamente el `dt` de la corrida.

### Cambié `table1.txt` y no sé qué resultados son válidos

Borre sólo la carpeta específica bajo `outputs/` y regenérela. Compare el hash
de `design_provenance.json` con el de la nueva tabla.

