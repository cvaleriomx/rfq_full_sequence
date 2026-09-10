# Secuencia de trabajo: TRANSOPTR → vanes en Warp → beam dynamics

## Cambio actual: apertura de salida en tres celdas equivalentes

El YAML `configs/mirfq1_sections.yaml` ahora define **TC → OM → FF**:
TC termina la modulación; OM abre las vanes sin modulación hasta el radio de
entrada del radial matcher, actualmente **1.2 cm = 12 mm**; FF queda sin metal.

```yaml
exit:
  enabled: true
  transition_cells: 1
  transition_length_factor: 0.5
  opening_cells: 3
  final_aperture_cm: rm_entrance
  unmodulated_length_cm: 2.0  # Solo se usa cuando opening_cells es 0.
  fringe_length_cm: 1.0
```

OM sustituye el tramo UM constante. Su longitud se fija al terminar TC como
`3 * beta(W_salida_TC) * lambda / 2` y se divide en tres segmentos iguales.
La apertura sigue `a = ai + (af-ai)*(10*u^3-15*u^4+6*u^5)`, con `u` medido a
lo largo de todo OM. Se evalúa en cada nodo del campo; no es una interpolación
lineal entre tres aperturas. La pendiente y segunda derivada de esta ley se
anulan en ambos extremos. El enfoque disminuye al crecer la apertura.

`final_aperture_cm: rm_entrance` sigue automáticamente la apertura inicial del
RM; también acepta un radio numérico en cm. Cambiar solo `transition_cells`
no configura esta apertura. El JSON generado contará `OM: 3`, omitirá UM y
reportará radios y longitud en `exit_opening`.

**Comprobado el 9 de septiembre de 2026:** pasaron las 20 pruebas y la corrida
real de TRANSOPTR terminó con `target_reached`, a **1.00423 MeV**. La apertura
OM mide **12.83474 cm** y crece de **4.48860 mm a 12 mm** de radio. Se regeneraron
los JSON, tablas y figuras de `outputs/mirfq1_sections/`. Warp también regeneró las vanes y `fieldmap.h5`, con residual DC de
**0.000842 V**, inferior a la tolerancia de 0.02 V. El matching del haz
completo sigue pendiente de validación.

Ejecutar desde la raíz, en este orden:

```bash
conda activate idp
source /home/cvalerio/work1/transoptr/transoptr-master/activate_transoptr.sh
python -m unittest discover -s tests -v
python -m rfq_transop_design.design_rfq configs/mirfq1_sections.yaml
python -m rfq_pipeline --config configs/mirfq1_warp.toml warp-vanes --dry-run
python -m rfq_pipeline --config configs/mirfq1_warp.toml warp-vanes
```

Revisar primero `phase_validation.json` y `vane_geometry.png` antes de construir
en Warp. El TOML de Warp ya no fija el extremo longitudinal: se obtiene del
perfil metálico nuevo más 2 cm de vacío. Los intervalos `nz` siguen siendo
configurables; revisar el nuevo paso longitudinal en el dry-run. La construcción
de Warp incluye OM porque sus filas representan metal, y excluye FF.

La secuencia del nuevo diseño tiene dos etapas separadas:

1. **Generar y verificar la RFQ con TRANSOPTR.** Implementado, incluidos perfiles
   ideales de las puntas y una terminación de salida aproximada.
2. **Construir las vanes 3D y preparar el campo en Warp.** Implementado con
   superficies trianguladas cerradas y solución DC; se comprueba su diferencia
   respecto al diseño de referencia.

Después se realiza el **beam dynamics** sobre la geometría y el campo validados.
Se ha comprobado el sincronismo de referencia con TRANSOPTR y la construcción
y solución DC en Warp. El tracking ya puede cargar el campo y las superficies
de pérdida. La equivalencia física del campo y la dinámica del haz completo
todavía requieren validación.

## Paso 1. Generar la RFQ con TRANSOPTR

### 1.1. Preparar el ambiente

Desde la raíz del repositorio, no desde esta subcarpeta:

```bash
cd /home/cvalerio/work1/machinlearning/rfq_full_sequence
conda activate idp
source /home/cvalerio/work1/transoptr/transoptr-master/activate_transoptr.sh
```

La primera vez que uses el proyecto en ese ambiente:

```bash
python -m pip install -e .
```

El script de activación define `optr` como alias. El diseñador llama al
`runoptr.sh` original por su ruta absoluta para usarlo desde Python.

### 1.2. Configurar el diseño por secciones

Edita [`configs/mirfq1_sections.yaml`](../configs/mirfq1_sections.yaml).
Este es el estudio actual con radial matcher y terminación de salida.

| Parámetro | Valor del estudio |
|---|---:|
| Energía inicial | 0.030 MeV = 30 keV |
| Energía objetivo | 1.000 MeV |
| Tolerancia energética | 0.010 MeV = 1 % |
| Frecuencia | 162 MHz |
| Voltaje intervane | 55.2 kV |
| Radial matcher | 4 celdas, configurable |
| Ejecutor | `runner: optr` — TRANSOPTR real |
| Salida | `outputs/mirfq1_sections/` |

Las secciones configuradas son:

| Sección | Función y parámetros |
|---|---|
| RM | Apertura decreciente, `m=1`, fase de −90°; adaptación radial propuesta |
| MS | Preparación inicial; crecimiento suave de modulación y enfoque |
| MB | Agrupamiento principal; evolución de fase y modulación |
| MBA | Agrupamiento y aceleración; transición del enfoque |
| MA | Aceleración principal con fase, modulación y enfoque objetivos constantes |
| TC | Transición del perfil modulado hacia simetría cuadrupolar |
| OM | Apertura suave sin modulación en tres segmentos; sustituye UM constante |
| FF | Región de campo de borde aproximado después del metal |

`scheme: nfsp` selecciona una **parametrización inspirada en NFSP**, no una
optimización NFSP completa del haz. `scheme: fsp` permite programar RM, SH, GB
 y ACC. Los tamaños, aperturas y rampas son parámetros de diseño que deben
validarse con el haz real. El enfoque B modifica la apertura y los coeficientes
del campo; no es únicamente una curva de reporte. Los indicadores aproximados
`Sig0T` y `Sig0L` del diseñador no sustituyen el cálculo de estabilidad del haz.

Las 4 celdas de RM son una elección inicial, no un máximo físico universal:
la literatura describe habitualmente 4–6 y también diseños más largos. Las
referencias se encuentran al final de esta guía.

Para otro estudio, copia el YAML y cambia también `output_dir`. Repetir una
corrida en la misma carpeta sobrescribe los archivos que vuelve a generar.

### 1.3. Ejecutar el diseñador

```bash
python -m rfq_transop_design.design_rfq configs/mirfq1_sections.yaml
```

Para cada celda, el diseñador usa la energía calculada anteriormente, propone
la geometría y ajusta su longitud con la fase medida por TRANSOPTR. La RFQ
acumulada se vuelve a calcular desde la energía inicial. La fase RF de
alimentación se mantiene fija; no se reinicia el reloj en cada celda.

La energía objetivo se comprueba después de incorporar la salida completa,
incluida la región FF. La corrida final también se verifica con distintos
muestreos del campo y una tolerancia de integración más estricta.

**No uses `example_design.yaml` para una corrida real:** selecciona `mock`, que
suma una estimación analítica y no ejecuta TRANSOPTR.

### 1.4. Revisar los resultados antes de pasar a Warp

Abre los archivos de [`outputs/mirfq1_sections/`](../outputs/mirfq1_sections/):

| Archivo | Qué revisar |
|---|---|
| `phase_validation.json` | `target_reached` y `phase_and_refinement_passed`; alcance de la validación |
| `rfq_sections.png` | Fase, modulación, enfoque y energía por sección |
| `phase_validation.png` | Fase y energía de una corrida completa |
| `vane_geometry.png` | Perfiles ideales de las cuatro puntas y detalle de salida |
| `vane_geometry.csv` | Coordenadas en cm; `metal_profile` distingue metal y región FF |
| `sections.json` | Límites y cantidades de cada sección |
| `fort.75` | Campo denso que se utilizó en TRANSOPTR |
| `table1.txt` | Tabla en el formato interno del diseñador |
| `table_pipeline.txt` | Esquema de tabla para el pipeline, con energía de la corrida completa |
| `history.csv` | Historial de las iteraciones, no el perfil de una única corrida final |
| `verification/` | Entradas y salidas de las comprobaciones numéricas |

La corrida comprobada obtuvo **1.00423 MeV**, dentro de la tolerancia del 1 %.
Tiene **127 celdas de RM/MS/MB/MBA/MA**, una transición TC, tres segmentos OM y un tramo FF:
132 registros en total. FF representa vacío con campo de borde; no es una
celda periódica ni una prolongación metálica de las vanes.

La longitud total modelada es **2.400545 m**, incluido FF. El extremo metálico
está en **2.390545 m**. El error máximo de fase en las regiones controladas es
menor de **0.145°**; la diferencia energética entre las dos comprobaciones más
finas fue **20 eV**. La fase síncrona no se exige en los tramos de salida sin
aceleración. Pasaron las 20 pruebas del proyecto durante esta implementación.

**Criterio para cerrar este paso:** aceptar el diseño de referencia y conservar
juntos su YAML, tabla, `fort.75`, geometría y validaciones. Esto no demuestra
transmisión, adaptación transversal del haz ni viabilidad de fabricación.

## Paso 2. Generar las vanes 3D y preparar el campo en Warp

### 2.1. Revisar la configuración

Desde la raíz, en el mismo ambiente `idp`:

```bash
conda activate idp
python -m rfq_pipeline --config configs/mirfq1_warp.toml warp-vanes --dry-run
```

[`configs/mirfq1_warp.toml`](../configs/mirfq1_warp.toml) selecciona el
`vane_geometry.csv` generado en el paso 1. Este comando revisa el dominio y
reporta dimensiones y almacenamiento del campo, sin importar Warp.

La configuración inicial usa puntas de radio transversal **2 mm**, cuerpo hasta
**18 mm** del eje, malla de **48 × 48 × 1160 intervalos** y frontera exterior
rectangular a tierra a **±22 mm**. Son elecciones iniciales para estudiar
convergencia, no dimensiones de fabricación optimizadas.

### 2.2. Construir e instalar los conductores y resolver el campo

```bash
python -m rfq_pipeline --config configs/mirfq1_warp.toml warp-vanes
```

El comando:

1. Lee exclusivamente las filas `metal_profile=true`, convierte cm a m y
   conserva RM, cuerpo, TC y OM. No construye metal en FF.
2. Forma cuatro sólidos cerrados: una cara de punta semicircular aproximada
   por segmentos, un cuerpo rectangular y tapas planas en los extremos.
   Comprueba el cierre de las mallas y exporta NPZ y STL en metros.
3. Instala las cuatro superficies con `Warp.Triangles`: **+27.6 kV** en ±x y
   **−27.6 kV** en ±y, correspondientes a 55.2 kV intervane.
4. Resuelve el problema electrostático con `MultiGrid3D` y comprueba su residual.
5. Exporta potencial y campo como `fieldmap.h5`, normalizados por voltaje
   intervane, y compara el campo axial con el potencial del diseño denso.

Warp se importa en modo serial mediante su opción `warpoptions`, evitando la
inicialización MPI. El comando debe ejecutarse en un proceso nuevo, como en los
ejemplos. La instalación de TRANSOPTR no hace falta en esta segunda etapa.

### 2.3. Revisar los archivos producidos

Todos se guardan en [`outputs/mirfq1_warp/`](../outputs/mirfq1_warp/):

| Archivo | Contenido |
|---|---|
| `vanes_3d.png` | Vista de los cuatro sólidos; los ejes usan escalas diferentes |
| `x_plus.stl`, `x_minus.stl`, `y_plus.stl`, `y_minus.stl` | Superficies cerradas para inspección en un visor CAD; unidades m |
| Archivos `.npz` del mismo nombre | Vértices y caras con orientación para Warp |
| `vanes_manifest.json` | Unidades, hash del perfil, extremos, geometría y voltajes |
| `warp_preparation.json` | Dominio, malla y convergencia del solucionador |
| `fieldmap.h5` | Campo DC calculado por los conductores, listo para cargar en tracking |
| `warp_dc_axis.csv` | Campo axial normalizado |
| `warp_dc_field.png` | Gráficas del campo DC |
| `warp_design_comparison.json` y `.png` | Diferencia frente al diseño de dos términos |

En la comprobación con OM, el solucionador convergió con residual de
aproximadamente **0.000842 V**. Las pruebas previas de pertenencia confirmaron
que el eje está libre y cada sólido ocupa su lado correcto. El campo axial
actual presentó correlación de **0.99905** y diferencia RMS relativa de **14.15 %**
respecto al diseño de dos términos. Una buena correlación no garantiza la
amplitud correcta.

**La construcción y la carga del campo están implementadas; la equivalencia
física todavía no está validada.** Hay que refinar malla, dominio y geometría
transversal antes de esperar que Warp reproduzca la energía de TRANSOPTR. Las
puntas circulares finitas y las tapas planas no son la superficie equipotencial
ideal del modelo de dos términos. El campo de borde ahora procede del problema
DC de esos sólidos y de la frontera exterior, en lugar de copiar la envolvente
FF aproximada del diseñador.

El comando antiguo `vanes` sigue siendo solo un exportador de perfiles;
**el nuevo comando que construye los conductores es `warp-vanes`**. No ejecutar
`all` con este TOML: es una configuración específica para la etapa Warp.

## Después: beam dynamics

Una vez generado el campo, revisar el plan:

```bash
python -m rfq_pipeline --config configs/mirfq1_warp.toml track --dry-run
```

Para ejecutar una prueba del haz:

```bash
python -m rfq_pipeline --config configs/mirfq1_warp.toml track
```

Este TOML carga `fieldmap.h5` y el `geometry_manifest`. El tracking instala un
`ParticleScraper` con pruebas de pertenencia equivalentes a los sólidos
poligonales de las vanes, para eliminar partículas que entran en metal. Se usan
pruebas de pertenencia rápidas para los barridos; el solucionador DC usa las
superficies trianguladas. El resumen registra `vane_scraper_installed`.

La configuración de prueba inyecta 1000 protones de 30 keV en z=0, con corriente
cero y radio de haz de 1 mm. La dependencia temporal del mapa es coseno; se usa
`rf_phase_deg=-180` para representar el seno con fase −90° de TRANSOPTR. La
referencia temporal está en la entrada de las vanes, no en el borde de la malla.

Se comprobó la inicialización y avance de **10 partículas durante 2 pasos**, con
el campo DC y el scraper activos. No se presenta esa prueba corta como una
validación de transmisión o energía de salida del haz completo.

Los resultados quedan en `particles_final.npz` y `particles_final.json`. Las
partículas que salen del dominio son absorbidas: `surviving_particles` no mide
por sí solo la transmisión. Falta consolidar los diagnósticos que separen salida
útil, pérdidas en metal y pérdidas por otras fronteras. Para resultados físicos
hay que definir el haz casado, estudiar carga espacial y convergencia y resolver
la diferencia de campo observada en el paso 2.

## Otros ejemplos del repositorio

- `configs/mirfq1_transoptr.yaml`: estudio anterior de fase continua sin la
  separación de secciones; obtuvo 1.00707 MeV. Se conserva como referencia.
- `configs/transoptr_smoke.yaml`: prueba de tres celdas con TRANSOPTR real.
- `configs/isac2_quick.toml`: ejemplo independiente sobre la tabla ISAC2 histórica.
  Sus comandos están en el [README principal](../README.md).
- `configs/mirfq1.toml`: no se ha conectado automáticamente al nuevo diseño.

## Referencias y límites del modelo

- [Zhang et al., PRAB 23, 042003 (2020)](https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.23.042003):
  ejemplo de RM típico de 4–6 celdas y uso de una región no modulada de salida.
- [Beam Dynamics of the TRIUMF ISAC RFQ, LINAC 1996](https://proceedings.jacow.org/l96/PAPERS/TUP32.pdf):
  ejemplo de adaptación radial extendida de 8 a 10 celdas.
- [Beam Dynamics Studies on the ISAC RFQ, PAC 1997](https://proceedings.jacow.org/pac97/papers/pdf/2W028.PDF):
  transición de salida hacia vanes no moduladas y necesidad de estudiar el campo de extremo.
- [HSI RFQ Upgrade, IPAC 2016](https://proceedings.jacow.org/ipac2016/papers/mopoy016.pdf):
  estrategias FSP/NFSP y optimización del diseño.
- [TRI-BN-22-07](https://beamphys.triumf.ca/~oshelb/physnotes/RFQISAC/RFQ.pdf):
  entradas RFQ y estudios de fase/energía con TRANSOPTR.

Las rampas programadas y el perfil de puntas son un modelo de diseño de dos
términos. La transición TC es una aproximación suave; FF usa una envolvente de
campo aproximada. No se presentan como una solución electrostática 3D del extremo
ni como una optimización NFSP completa de aceptación y estabilidad del haz.

## Diagnósticos de óptica del haz con TRANSOPTR

Cada verificación final del diseño por secciones produce ahora:

| Archivo en `outputs/mirfq1_sections/` | Contenido |
|---|---|
| `beam_optics.csv` | Envolventes, divergencias, correlaciones, Twiss proyectados, emitancias normalizadas, apertura y ocupación, fase transversal acumulada y distribución longitudinal lineal |
| `beam_phase_advance.csv` | Avance transversal nativo por pares de celdas regulares, con posición y sección |
| `beam_optics_summary.json` | Máximos globales y por sección, posiciones críticas, límites configurados y advertencias de validez |
| `beam_optics.png` | Gráficas de avance de fase, ocupación, divergencia, alpha, anchura de fase, dispersión energética, tamaño y mismatch |

Para regenerar **solo estos diagnósticos**, sin ejecutar TRANSOPTR ni Warp:

```bash
conda activate idp
python -m rfq_transop_design.beam_diagnostics configs/mirfq1_sections.yaml
```

El comando usa la verificación guardada con más puntos por celda. Mantener juntos
los archivos de esa corrida (`history.csv`, `sections.json`, `vane_geometry.csv`
y `verification/`); no mezclar perfiles o configuraciones de diseños distintos.

Opciones bajo `beam_diagnostics` en el YAML:

```yaml
beam_diagnostics:
  envelope_multiplier: 1.0
  envelope_convention: unconfirmed
  max_occupancy: 1.0
  max_divergence_mrad: null
  max_pair_advance_deg: null
  matched_reference_csv: null
```

Los límites numéricos generan advertencias; no cambian ni optimizan el diseño.
Los límites de divergencia y avance quedan sin fijar hasta definir los requisitos
ópticos. La ocupación es `multiplicador * envolvente / apertura`, para un haz
centrado. Se excluye FF porque no hay metal. No se interpreta el centroide de
fase interno del modelo RFQ como una desviación física del haz.

La escala de las envolventes depende de la entrada de TRANSOPTR. No llamar a un
factor 3 «3 sigma» hasta confirmar que la entrada representa momentos RMS.
El máximo informado es el máximo **muestreado de la envolvente**, no el de todas
las partículas. La comparación con apertura no calcula transmisión ni halo.

El avance se lee de `fort.81` y `fort.83`: `VECTIVE.f` escribe la fase en vueltas,
que se convierten a grados. Se restan valores en los extremos de pares de celdas
regulares consecutivas, contados desde la entrada. Se omiten pares que crucen
secciones y TC/OM/UM/FF. Es avance de la óptica transportada; no es la fase
síncrona RF, ni una sintonía periódica calculada por autovalores. Los campos
`Sig0T_deg`/`Sig0L_deg` del historial no se usan como sustitutos.

Los Twiss proyectados se calculan con los segundos momentos y la emitancia
geométrica; la emitancia normalizada incluye el factor relativista beta*gamma.
Las conversiones longitudinales son locales y lineales: anchura temporal a partir
de longitud/velocidad y dispersión energética a partir de dp/p. Una envolvente
que abarca gran parte del ciclo RF requiere un modelo de partículas para estudiar
captura. Tampoco una emitancia conservada por un modelo lineal prueba ausencia
de crecimiento real por no linealidades.

Para calcular mismatch Bmag, proporcionar `matched_reference_csv` con columnas
`z_cm,x_alpha,x_beta_cm,y_alpha,y_beta_cm`, en las mismas coordenadas proyectadas.
La referencia debe ser una solución adaptada obtenida independientemente;
Bmag=1 representa coincidencia de Twiss en el caso desacoplado. No se extrapola
fuera de su rango. Sin referencia, el resultado queda vacío y se marca como no
evaluado. La depresión de sintonía queda pendiente de corridas compatibles con y
sin carga espacial. Los Twiss por plano requieren cautela si hay acoplamiento.

**Hallazgos de la corrida actual:** la divergencia inicial x es cero y produce
emitancia x inicial nula. Por ello sus Twiss y avance de fase se dejan vacíos,
aunque TRANSOPTR emita números por redondeo posteriormente. También se detecta
una envolvente longitudinal superior a 180 grados RF. No se modificó el haz de
entrada para ocultar estas condiciones: energía y sincronismo aprobados no
constituyen una validación de la óptica del haz completo.

## Comparación de anchos longitudinales de 180° y 90°

Se ejecutaron tres verificaciones con la **misma geometría**: la entrada original
(semilongitud 28 cm), ancho completo de 180° y ancho completo de 90°. Cada caso
usa 16, 32 y 64 nodos de campo por celda, con tolerancia RK diez veces menor en
el último. Los anchos son completos del elipsoide longitudinal nativo, no RMS
ni FWHM; a 30 keV y 162 MHz la conversión es:

`semilongitud_cm = beta * c_cm_s / frecuencia_hz * ancho_completo_deg / 720`.

| Entrada | Semilongitud [cm] | Energía final con 64 puntos [MeV] | Error máximo de fase [°] | Diferencia energética 32→64 [eV] |
|---|---:|---:|---:|---:|
| Original, ≈13623.35° completos | 28.000000 | 1.00419 | 0.144686 | 20 |
| 180° completos | 0.369953 | 1.00425 | 0.144129 | 20 |
| 90° completos | 0.184976 | 1.00425 | 0.144129 | 20 |

Los tres casos pasan los criterios existentes de energía y fase/refinamiento.
Entre 32 y 64 puntos, la diferencia máxima de cada envolvente, dividida por el
máximo de esa envolvente en la corrida fina, es menor de 0.5 %. La disminución
de las diferencias no es monótona en todas las coordenadas; estas tres
resoluciones no establecen un orden de convergencia.

Los resultados están en `outputs/mirfq1_width_comparison/`:

- `comparison.png`, `comparison.csv` y `comparison.json`: superposición de
  energía y envolventes, valores finales y condiciones del estudio.
- `refinement_comparison.csv`: diferencias de energía, fase y envolventes para
  cada resolución, comparadas con 64 puntos en posiciones comunes.
- `baseline/`, `width_180deg/` y `width_90deg/`: cada carpeta contiene `input.dat`,
  `config.yaml`, las tres verificaciones y los diagnósticos completos de óptica.

Para repetir la comparación desde la raíz:

```bash
conda activate idp
source /home/cvalerio/work1/transoptr/transoptr-master/activate_transoptr.sh
python -m rfq_transop_design.compare_bunch_widths configs/mirfq1_sections.yaml
```

El programa toma la geometría de `history.csv` y la entrada de la verificación
base guardada; no vuelve a optimizar las celdas. Conserva carga de entrada cero,
las dimensiones transversales originales y el sexto componente longitudinal.
Por tanto, **este estudio todavía no usa 1 mA ni 0.25 mm·mrad RMS normalizados**.
La emitancia x inicial sigue siendo cero. Reducir la longitud cambia también la
emitancia longitudinal, pues se mantiene el sexto componente y la correlación.

Con carga cero, el movimiento de referencia no depende físicamente del tamaño
del bunch en este modelo; los 60 eV entre la referencia y los casos reducidos
son una sensibilidad numérica a la integración de las distintas envolventes.
La prueba confirma una precisión numérica similar bajo estas condiciones;
no valida captura de un haz continuo, transporte con carga espacial ni estabilidad
de Hofmann.

## Ancho longitudinal predeterminado: 90° completos

El diseñador define `beam.longitudinal_full_width_deg: 90.0` por defecto.
`configs/mirfq1_sections.yaml` lo declara explícitamente y cada nueva entrada
TRANSOPTR calcula su semilongitud usando la energía inicial, masa y frecuencia.
Se respeta la escala longitudinal CONX del archivo de entrada.

Son ±45° respecto al centro del bunch, no una anchura RMS ni FWHM. Para 30 keV
y 162 MHz: semilongitud 0.184976481 cm, longitud completa 3.699530 mm.
El template `configs/mirfq1_transoptr.dat` también contiene ese valor de referencia.
El ancho del bunch no cambia la fase RF de su centro.

La dispersión de energía a lo largo de la RFQ sí puede cambiar al cambiar la
extensión longitudinal y su correlación con energía. Esta opción conserva el
sexto componente de entrada y las correlaciones del template; por ello también
cambia la emitancia longitudinal inicial. No introduce automáticamente la corriente
de 1 mA ni la emitancia transversal solicitada para el estudio posterior.

Para usar expresamente el ancho del template, definir
`longitudinal_full_width_deg: null`. El comparador fija por separado 180°, 90° y
la entrada histórica para que el nuevo predeterminado no altere su comparación.
Los resultados históricos no se regeneraron al cambiar este valor predeterminado.

## Corrección de inicialización: Twiss y emitancias del template

Las configuraciones MIRFQ ahora seleccionan:

```yaml
beam:
  initialization: template_twiss
  longitudinal_full_width_deg: 90.0
```

Antes, las líneas de parámetros alpha/beta/emitancia no se aplicaban. El cambio
180°→90° anterior modificaba solo la quinta dimensión del renglón 4: mantenía
el sexto componente y las correlaciones anteriores. Esa comparación permanece
como estudio histórico de aquella entrada, no de la entrada Twiss corregida.

Ahora `prepare_beam_input()` lee el bloque MIRFQ de 11 parámetros (voltaje, fase,
seguido de alpha/beta/emitancia para x, y, z) y genera las seis dimensiones y las
tres correlaciones antes de ejecutar TRANSOPTR. Usa las mismas conversiones que
`CIC3.f`, con las escalas CONX del template. Requiere parámetros fijos y un template
sin correlaciones previas: la inicialización es desacoplada entre planos.

Para cada plano: `size=sqrt(beta*emit)`,
`spread=sqrt((1+alpha²)*emit/beta)`, `correlation=-alpha/sqrt(1+alpha²)`.
Los valores de emitancia del archivo son geométricos, en la convención de la
envolvente TRANSOPTR; no se reinterpretan como RMS normalizados.

En longitudinal hay una restricción: no se pueden imponer simultáneamente beta_z,
emitancia_z y un ancho independiente. Con el ancho completo fijado a 90° se
conservan alpha_z=2 y emitancia_z=0.00049 cm·rad, y se calcula
`beta_z = semilongitud² / emitancia_z`:

| Propiedad | Valor efectivo |
|---|---:|
| Emitancia x / y | 0.003 / 0.003 cm·rad |
| Beta x / y | 31.5 / 33.6 cm |
| Beta z declarado | 47.27 cm |
| Beta z efectivo para 90° | 69.82918 cm |
| Semilongitud z | 0.184976481 cm |
| Sexto componente longitudinal | 0.00592331145 |
| Correlación z–sexto componente | −0.894427191 |

Para respetar también beta_z=47.27 cm exactamente, poner
`longitudinal_full_width_deg: null`: el ancho será el que resulte de beta_z y
emitancia_z. Para reproducir una entrada antigua basada en seis dimensiones,
usar `initialization: dimensions`. El comparador histórico selecciona este último
modo expresamente para conservar su experimento original.

Cada corrida guarda `beam_input.json` junto a `data.dat` con los valores declarados
y efectivos. La verificación final también lo copia al directorio principal.
Los valores del YAML de voltaje, fase RF y energía siguen controlando el diseño.

Se comprobó la corrección con TRANSOPTR sobre la geometría existente, a carga cero,
en `outputs/mirfq1_twiss_validation/`, usando 16/32/64 puntos por celda. La energía
final fina fue 1.00419 MeV y la diferencia energética 32→64 fue 70 eV. El máximo
error de fase fue 0.150812°: pasa el criterio de verificación de 0.5°, aunque queda
ligeramente por encima de los 0.15° usados al ajustar las celdas. No se reoptimizó
la geometría. La entrada efectiva ya muestra emitancias x/y no nulas y permite
exportar sus Twiss y avances de fase. Estos parámetros aún no equivalen a la
entrada futura de 1 mA y 0.25 mm·mrad RMS normalizados. Pasaron 29 pruebas.

## Generar entrada RMS y comparar 0 mA / 1 mA

El comando adicional toma la geometría existente y genera dos entradas con
**0.25 mm·mrad RMS normalizados en cada plano**, sin factor π. Conserva los Twiss
transversales y alpha_z/emitancia_z del template, y fija el ancho completo del
elipsoide longitudinal a 90°. No optimiza las celdas.

Solo preparar archivos, sin ejecutar TRANSOPTR:

```bash
python -m rfq_transop_design.compare_current configs/mirfq1_sections.yaml --prepare-only
```

Preparar y ejecutar la comparación, desde la raíz:

```bash
conda activate idp
source /home/cvalerio/work1/transoptr/transoptr-master/activate_transoptr.sh
python -m rfq_transop_design.compare_current configs/mirfq1_sections.yaml --current-ma 1 --emittance-n-rms 0.25
```

`--emittance-n-rms` está en mm·mrad por plano; `--current-ma` es la corriente
promedio en mA. La salida predeterminada es `outputs/mirfq1_current_comparison/`.
Cada carpeta `zero_current/` y `with_current/` contiene:

- `template.dat`: parámetros Twiss con emitancias convertidas.
- `prepared_data.dat`: dimensiones y correlaciones efectivas listas para modo 5.
- `input_parameters.json`: valores pedidos, carga y conversiones.
- `config.yaml`: configuración de cada caso.
- `verification/points_16`, `points_32`, `points_64`: los `data.dat` y
  `fort.envelope` efectivos, y la verificación de refinamiento.
- Los diagnósticos completos de fase, apertura y óptica y `beam_envelope_full.csv`.

En la raíz de la comparación: `comparison.csv`, `comparison.json`,
`refinement_comparison.csv` y `comparison.png`. Las gráficas comparativas de
este comando muestran tamaños y dispersión energética **RMS**; los CSV nativos de
óptica conservan su convención de envolvente, identificada como `sqrt5_rms`.

### Conversiones y alcance

Modo 5 espera carga por bunch: `Q = I_promedio / f`. Para 1 mA a 162 MHz,
Q=6.1728395 pC. Se supone un bunch por período y toda la corriente dentro de él.
Esto modela un haz ya agrupado de 90°, no la captura del haz continuo original.

En modo 5 las dimensiones nativas son sqrt(5) veces las RMS y las emitancias
nativas son cinco veces las RMS geométricas. Se usa:

`emitancia_nativa_cm_rad = 5 * emitancia_normalizada_RMS_mm_mrad * 1e-4 / (beta*gamma)`.

Los dos casos tienen exactamente la misma matriz inicial; únicamente cambia la
carga. Se verifican las emitancias iniciales directamente en `fort.envelope` y la
carga en el `data.dat` usado. La emitancia longitudinal sigue siendo la geométrica
nativa del template. No se vuelve a interpretar ese valor como normalizado RMS.

Se encontró en el RFQ.f instalado que QSC usaba BNCHARGE sin inicializar, en vez
de CURRENT. La copia local [RFQSC.f](transoptr/RFQSC.f) corrige esa referencia y
renombra las rutinas. Ambos casos la compilan dentro de su sy.f; el programa
TRANSOPTR instalado permanece intacto. Su procedencia y licencia están en
[transoptr/README.md](transoptr/README.md). Usar solo una carga distinta con la
rutina antigua no garantizaría que realmente se aplicara carga espacial.

### Resultados comprobados

| Magnitud, resolución fina | 0 mA | 1 mA |
|---|---:|---:|
| Energía de referencia [MeV] | 1.00420 | 1.00418 |
| Diferencia energética 32→64 [eV] | 10 | 40 |
| Dispersión energética RMS a la salida [keV] | 4.177 | 7.464 |
| Longitud RMS a la salida [mm] | 0.3820 | 0.3243 |
| Máxima ocupación nativa x | 1.800 | 1.795 |
| Máxima ocupación nativa y | 1.559 | 1.609 |

Pasaron los criterios de fase y refinamiento de energía. Las diferencias de
las envolventes 32→64 respecto al máximo de la corrida fina fueron inferiores
a 0.36 %. Las 31 pruebas del código pasaron y se detectó un efecto no nulo de
la corriente en las envolventes.

**La envolvente equivalente excede la apertura en ambos casos.** La integración
no elimina partículas al tocar las vanes y continúa propagando sus momentos:
los resultados de salida son del modelo de envolventes, no de un haz transmitido
validado. La matriz inicial de este estudio tiene emitancias diferentes de los
antiguos ejemplos de 0.003 cm·rad nativos. Esto requiere revisar adaptación y
contención antes de usar estos resultados como diseño aceptado.
