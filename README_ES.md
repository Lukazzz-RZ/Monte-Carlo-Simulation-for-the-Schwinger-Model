# Monte Carlo Simulation of the Dual Schwinger Model

[English](README.md)

Implementacion en C para simular la representacion dual del modelo de
Schwinger en 1+1 dimensiones mediante actualizaciones de plaqueta y algoritmos
de gusano.

> Proyecto en desarrollo asociado a un Trabajo de Fin de Grado.

## Requisitos

- Windows 10 u 11.
- Visual Studio 2022 con soporte para C.
- CMake 3.20 o posterior.

## Estructura del proyecto

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

Abrir una terminal de Visual Studio en la carpeta del proyecto y ejecutar:

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

Al iniciar, el programa muestra un menu con las siguientes opciones:

1. Prueba rapida.
2. Barrido en theta.
3. Barrido en masa.
4. Plano masa-theta.
5. Punto unico.
6. Salir.

Cada opcion muestra los parametros cargados, el numero de series de
simulacion, los pasos de termalizacion, el numero de medidas y la separacion
entre medidas.

La opcion seleccionada debe confirmarse con `y` o cancelarse con `n`.

## Paso Monte Carlo

Cada paso Monte Carlo contiene:

1. Un barrido secuencial completo de plaquetas.
2. `max(1, V / 10)` intentos de gusano.
3. Un segundo barrido secuencial completo de plaquetas.

## Resultados

Los resultados se guardan automaticamente en:

```text
Results/
```

Los archivos `Series_*.dat` contienen:

- los observables medidos;
- el signo de la configuracion;
- los productos entre el signo y cada observable para la reponderacion.

Los archivos `Thermalization_*.txt` contienen la evolucion de:

- la plaqueta media;
- la densidad de dimeros.

## Configuracion

Los parametros generales del reticulo se definen en:

```text
include/Lattice.h
```

Los preajustes de simulacion mostrados en el menu se definen en:

```text
src/Main.c
```

## Estado

El codigo se encuentra actualmente en desarrollo. La memoria completa y la
documentacion academica se publicaran una vez finalizado el proceso de
evaluacion.
