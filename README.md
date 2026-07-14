# Schwinger Worm

Codigo en C para simular la representacion dual del modelo de Schwinger en
1+1 dimensiones mediante actualizaciones de plaqueta y algoritmos de gusano.

> Proyecto en desarrollo asociado a un Trabajo de Fin de Grado.

## Requisitos

- Windows 10 u 11.
- Visual Studio 2022 con soporte para C.
- CMake 3.20 o posterior.

## Estructura

```text
SchwingerWorm/
|-- CMakeLists.txt
|-- include/
|   `-- Lattice.h
|-- src/
|   |-- Main.c
|   |-- Initializer.c
|   |-- UpdCore.c
|   |-- Plaq_mod.c
|   |-- Worm_ferm.c
|   |-- Measures.c
|   |-- Visualizer.c
|   `-- Other.c
`-- Results/
```

## Compilacion

Desde una terminal de Visual Studio abierta en la carpeta del proyecto:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

El ejecutable se genera normalmente en:

```text
build/Release/SchwingerWorm.exe
```

## Uso

Ejecutar:

```powershell
.\build\Release\SchwingerWorm.exe
```

Al iniciar se muestra un menu con varias opciones:

1. Prueba rapida.
2. Barrido en theta.
3. Barrido en masa.
4. Plano masa-theta.
5. Punto unico.
6. Salir.

Debajo de cada opcion se muestran los parametros cargados, el numero de
series, los pasos de termalizacion, el numero de medidas y la separacion entre
medidas.

La opcion seleccionada debe confirmarse mediante `y` o cancelarse mediante
`n`.

## Paso Monte Carlo

Cada paso Monte Carlo contiene:

1. Un barrido secuencial completo de plaquetas.
2. `max(1, V / 10)` intentos de gusano.
3. Un segundo barrido secuencial completo de plaquetas.

## Resultados

Los resultados se guardan automaticamente en la carpeta:

```text
Results/
```

Los archivos `Series_*.dat` contienen:

- observables medidos;
- signo de la configuracion;
- productos entre el signo y cada observable para la reponderacion.

Los archivos `Thermalization_*.txt` contienen la evolucion de:

- la plaqueta media;
- la densidad de dimeros.

## Configuracion

Los parametros generales del reticulo se encuentran en:

```text
include/Lattice.h
```

Los preajustes mostrados en el menu se encuentran en:

```text
src/Main.c
```

## Estado

El codigo se encuentra en desarrollo. La memoria y la documentacion academica
completa se publicaran cuando finalice el proceso de evaluacion.
