# Contratos de datos

## Entrada de diseño

La entrada canónica es una tabla separada por espacios con, como mínimo, las
columnas `Cell`, `V`, `Wsyn`, `A10`, `Phi`, `a`, `m`, `L` y `Z`. Las longitudes
de entrada están en centímetros. `Cell=0` es la condición inicial y `Z=0`.

## Eje electrostático exportado por Warp

CSV con encabezado `z_m,ez_vpm_per_v`, `z_m` creciente y campo dividido por la
diferencia intervane aplicada. Esta interfaz evita acoplar la comparación a una
versión concreta de la API interna de Warp.

## Mapa externo

HDF5 `rfq-fieldmap-v1`, ejes en metros, campo en V/m por volt intervane y orden
de los arrays `(x,y,z)`. Los nombres de datasets se documentan en el README.

## Tracking

La salida comprimida `particles_final.npz` contiene `x_m`, `y_m`, `z_m`,
`vx_mps`, `vy_mps` y `vz_mps`. El JSON homónimo registra los parámetros y el
número de supervivientes.

