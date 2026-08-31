# Prueba ML: distribución de carga y subrogado de Coulomb

Este prototipo reproduce la idea central de Edelen et al.: generar una muestra
dispersa con el modelo físico y entrenar una red neuronal que lo evalúe mucho
más rápido. No intenta reproducir OPAL ni un acelerador completo.

El código genera tres familias de carga:

- esfera uniforme;
- gaussiana 3D, esférica o elipsoidal;
- proyección espacial tipo KV, uniforme dentro de un elipsoide.

La última opción **no es un KV completo**. La distribución KV rigurosa vive en
el espacio de fases transversal 4D. Aquí sólo se modela su densidad espacial,
que es la parte relevante para esta primera prueba electrostática.

## Idea general del programa

El programa separa el trabajo en dos partes:

1. **Simulador físico:** coloca macropartículas, reparte la carga total entre
   ellas y suma sus contribuciones mediante la ley de Coulomb.
2. **Modelo subrogado:** aprende la relación entre los parámetros iniciales y
   los resultados del simulador. Después puede aproximar esos resultados sin
   volver a sumar las contribuciones de todas las partículas.

El flujo completo es:

```text
Parámetros físicos (geometría, Q, a, b, c)
                  │
                  ▼
       Generación de macropartículas
                  │
                  ▼
       Suma de Coulomb en las sondas
                  │
                  ▼
            Dataset: X → y
                  │
                  ▼
      Entrenamiento de la red neuronal
                  │
                  ▼
  Predicciones rápidas para parámetros nuevos
```

El subrogado no mueve partículas ni resuelve nuevamente la ley de Coulomb. Lo
que hace es interpolar la función física que observó durante el entrenamiento.
Por ello solamente debe usarse dentro de los intervalos de carga y tamaños que
aparecen en su dataset.

## ¿Qué es el dataset?

Un **dataset** es una colección organizada de ejemplos usados para enseñar y
evaluar el modelo. En este proyecto, una fila del dataset representa **una
configuración física completa**, no una partícula individual.

Por ejemplo, la entrada

```text
gaussian, Q=3 nC, sigma_x=2.5 mm, sigma_y=3.5 mm, sigma_z=5.0 mm
```

se almacena en el vector:

```text
X = [0, 1, 0, 3.0, 2.5, 3.5, 5.0]
```

Los tres primeros valores codifican la geometría:

| Geometría | `is_sphere` | `is_gaussian` | `is_kv` |
|---|---:|---:|---:|
| Esfera | 1 | 0 | 0 |
| Gaussiana | 0 | 1 | 0 |
| KV espacial | 0 | 0 | 1 |

Los cuatro valores restantes son la magnitud de la carga total en nC y los
tres tamaños en mm. Para una esfera se cumple `a = b = c = R`. Para una
gaussiana son `sigma_x`, `sigma_y` y `sigma_z`. Para KV espacial son los
semiejes del elipsoide.

Para esa entrada, el simulador produce el vector objetivo `y`:

| Columna de `y` | Significado |
|---|---|
| `E_x_0p5_V_m` | Componente E_x en `(0.5a, 0, 0)` |
| `E_x_1p0_V_m` | Componente E_x en `(a, 0, 0)` |
| `E_x_2p0_V_m` | Componente E_x en `(2a, 0, 0)` |
| `E_y_0p5_V_m` | Componente E_y en `(0, 0.5b, 0)` |
| `E_y_1p0_V_m` | Componente E_y en `(0, b, 0)` |
| `E_y_2p0_V_m` | Componente E_y en `(0, 2b, 0)` |
| `E_z_0p5_V_m` | Componente E_z en `(0, 0, 0.5c)` |
| `E_z_1p0_V_m` | Componente E_z en `(0, 0, c)` |
| `E_z_2p0_V_m` | Componente E_z en `(0, 0, 2c)` |
| `V_center_V` | Potencial eléctrico en el origen |

Por tanto, un dataset de 500 configuraciones tiene normalmente:

```text
X.shape = (500, 7)    # 500 casos, 7 entradas
y.shape = (500, 10)   # 500 casos, 10 resultados físicos
```

Esto no significa que se usaron sólo 500 partículas. Si se ejecuta con
`--particles 10000`, cada una de las 500 filas se calcula internamente con
10 000 macropartículas.

El archivo `dataset_coulomb.npz` es un contenedor comprimido de NumPy y guarda:

- `X`: entradas de la red;
- `y`: respuestas calculadas con Coulomb;
- `metadata_json`: semilla, número de partículas, intervalos de parámetros,
  suavizado, nombres de columnas y demás información necesaria para saber
  cómo se produjo el dataset.

El número de partículas `N` no aparece como columna de `X`. El subrogado
representa el valor de `N` usado al generar sus etiquetas. Para estudiar otro
valor, por ejemplo 100 000 partículas, se debe producir un nuevo dataset o
primero demostrar mediante convergencia que ambos valores son equivalentes
dentro de la precisión requerida.

## Cómo se calcula una fila del dataset

Para cada configuración, el código hace lo siguiente:

1. Selecciona `Q`, `a`, `b` y `c` mediante un muestreo Latin hypercube. Las
   tres geometrías aparecen de forma balanceada.
2. Genera una nube unidad con una secuencia Halton de baja discrepancia.
3. Refleja los puntos en los ocho octantes (`quiet_start`) para que el centro
   de carga sea exactamente cero y reducir ruido estadístico.
4. Escala la nube unidad para obtener la esfera, gaussiana o elipsoide físico.
5. Divide la carga total entre las partículas:

   ```text
   q_macro = Q / N
   ```

6. Coloca las nueve sondas y suma el campo suavizado de Coulomb:

   ```text
   E(r) = k Σ q_macro (r-r_i) / (|r-r_i|² + epsilon²)^(3/2)
   ```

7. Calcula también el potencial eléctrico en el centro.
8. Guarda los parámetros en una fila de `X` y los resultados en la fila
   correspondiente de `y`.

Se reutiliza la misma nube unidad para todas las configuraciones de una
geometría. De este modo, al modificar `Q`, `a`, `b` o `c`, cambia la física
pero no aparece una nueva fluctuación aleatoria que confunda a la red.

La energía electrostática de todos los pares puede añadirse con
`--include-energy`, pero esa operación requiere O(N²) interacciones. Los diez
observables normales sólo evalúan un número fijo de sondas y cuestan O(N).

## Estructura del código

| Archivo | Responsabilidad |
|---|---|
| `prueba_ml/distributions.py` | Genera y escala esfera, gaussiana y KV espacial |
| `prueba_ml/physics.py` | Campo, potencial, energía y soluciones analíticas |
| `prueba_ml/dataset.py` | Muestrea parámetros y construye `X`, `y` y metadatos |
| `prueba_ml/model.py` | Construye, entrena, evalúa y guarda la red neuronal |
| `estudio_convergencia.py` | Examina cómo cambia la solución al aumentar N |
| `generar_dataset.py` | Interfaz para producir un archivo `.npz` |
| `entrenar_modelo.py` | Entrena el subrogado y crea métricas y paridad |
| `predecir.py` | Carga un modelo y predice una configuración nueva |
| `graficar_campo_eje_y.py` | Grafica los puntos de E_y sobre el eje y usando sólo el modelo |
| `ejecutar_demo.py` | Ejecuta generación, entrenamiento y predicción juntos |

## Ambiente

El ambiente encontrado en esta máquina es:

```bash
conda activate warp_2
cd examples/surrogate_coulomb
```

También se puede llamar directamente a su Python:

```bash
python ejecutar_demo.py
```

## Primera ejecución

La prueba rápida usa 300 configuraciones y 2 000 macropartículas:

```bash
python ejecutar_demo.py
```

La configuración recomendada usa 500 simulaciones y 10 000 partículas:

```bash
python ejecutar_demo.py --full
```

Los resultados aparecen en `resultados/demo/`: dataset, modelo serializado,
métricas JSON y gráficas de paridad.

## Flujo por etapas

### 1. Estudio de convergencia

```bash
python estudio_convergencia.py
```

Para esfera y gaussiana esférica se compara contra el campo continuo
analítico. Para el elipsoide tipo KV se usa el mayor N como referencia. Los
resultados quedan en `resultados/convergencia/`.

Se puede cambiar la secuencia de partículas:

```bash
python estudio_convergencia.py --particles 1024,2048,5000,10000,20000
```

### 2. Dataset de producción

```bash
python generar_dataset.py \
  --samples 500 \
  --particles 10000 \
  --output datos/dataset_coulomb.npz
```

Cada ejemplo contiene siete entradas, explicadas con detalle en la sección
“¿Qué es el dataset?”:

```text
tipo de geometría, Q, a, b, c
```

El tipo se codifica con tres variables *one-hot*. Las salidas son el campo
axial en sondas situadas a 0.5, 1 y 2 tamaños característicos sobre cada eje,
más el potencial en el centro.

Las sondas se mueven con los semiejes de la distribución. Para una gaussiana,
por ejemplo, se evalúan en 0.5 sigma, 1 sigma y 2 sigma.

Opcionalmente se puede añadir la energía de todos los pares:

```bash
python generar_dataset.py --samples 100 --particles 2000 --include-energy
```

Esta opción cuesta O(N²). Se recomienda probarla primero con N pequeño. El
cálculo estándar del campo cuesta O(N) porque el número de sondas es fijo.

### 3. Entrenamiento

```bash
python entrenar_modelo.py \
  --dataset datos/dataset_coulomb.npz \
  --output resultados/modelo_subrogado.joblib
```

La red es una MLP con cuatro capas ocultas de 32 nodos y activación `tanh`.
Las entradas se estandarizan. Como los observables físicos son positivos, las
salidas se transforman con logaritmo y después se estandarizan; esto mejora el
error relativo y evita predicciones de campo negativas. Se reserva 20 % de los
ejemplos para una prueba que no participa en el entrenamiento.

En otras palabras, con 500 filas se usan 400 para ajustar la red y 100 para
medir si puede predecir casos que no vio. Dentro de las 400 filas de
entrenamiento, scikit-learn reserva además una pequeña fracción para detener el
entrenamiento cuando el error deja de mejorar.

La gráfica `paridad.png` muestra simulación en el eje horizontal y predicción
del subrogado en el vertical. Una predicción perfecta cae sobre la diagonal.
Un punto arriba de la diagonal es una sobreestimación y un punto abajo es una
subestimación.

Una buena paridad sólo demuestra que la red imita al simulador. El estudio de
convergencia es el que determina si el simulador de macropartículas se aproxima
bien a la solución física continua.

### 4. Predicción

```bash
python predecir.py \
  --model resultados/modelo_subrogado.joblib \
  --geometry gaussian \
  --charge-nc 3.0 \
  --a-mm 2.5 --b-mm 3.5 --c-mm 5.0
```

Para una esfera sólo se proporciona el radio:

```bash
python predecir.py \
  --model resultados/modelo_subrogado.joblib \
  --geometry sphere --charge-nc 2.0 --a-mm 5.0
```

### 5. Gráfica 2D del campo sobre el eje y

El siguiente comando carga exclusivamente el modelo entrenado; no genera ni
suma macropartículas:

```bash
python graficar_campo_eje_y.py \
  --model resultados/demo/modelo_subrogado.joblib \
  --geometry gaussian \
  --charge-nc 3.0 \
  --a-mm 2.5 --b-mm 3.5 --c-mm 5.0 \
  --output resultados/campo_eje_y.png
```

El modelo predice `E_y` en `y=0.5b`, `y=b` y `y=2b`. Como todas las
distribuciones de este prototipo están centradas y son simétricas, el script
usa

```text
E_y(-y) = -E_y(y),    E_y(0) = 0
```

para mostrar también la mitad negativa. Guarda la figura PNG y los siete
puntos en `resultados/campo_eje_y.csv`.

La línea discontinua sólo conecta los puntos como guía visual. No representa
una predicción continua entre ellos. Para obtener cien o más posiciones, o un
mapa espacial en el plano `xy`, se necesita entrenar otro subrogado que reciba
la coordenada de consulta:

```text
(geometría, Q, a, b, c, x, y, z) -> (E_x, E_y, E_z, V)
```

Interpolar los tres puntos actuales no añade información física nueva y puede
ser especialmente incorrecto cerca del centro o de la superficie de una
esfera.

## Decisiones numéricas importantes

- `quiet_start` refleja partículas en los ocho octantes. Esto elimina el
  desplazamiento del centro de carga y reduce mucho el ruido dipolar.
- Se conserva la carga total y cada macropartícula tiene carga Q/N.
- El suavizado de Plummer por defecto es 1.25 veces el espaciamiento medio
  característico. Se guarda en los metadatos del dataset.
- Se reutiliza una nube unidad por geometría para todo el dataset. Esta técnica
  de números aleatorios comunes evita entregar etiquetas ruidosas a la red.
- N no es una entrada del subrogado. Debe fijarse después del estudio de
  convergencia y mantenerse igual al generar el dataset.

## Pruebas

```bash
python -m unittest discover -s tests -v
```

Las pruebas verifican simetría, RMS de la esfera, energía de dos partículas,
campo de una esfera contra la solución analítica y generación del dataset.

## Siguiente extensión física

Una vez validado este ejemplo, el paso natural es sustituir los diez
observables por cantidades del RFQ: tamaños RMS, emitancias, transmisión,
halo, pérdidas y energía. Los parámetros geométricos y de operación del RFQ
serían entonces las entradas del subrogado.

Referencia: [Machine learning for orders of magnitude speedup in multiobjective
optimization of particle accelerator systems](https://doi.org/10.1103/PhysRevAccelBeams.23.044601).
