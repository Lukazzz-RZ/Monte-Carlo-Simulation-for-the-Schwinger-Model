# Simulación Monte Carlo del modelo dual de Schwinger

[English](README.md)

Implementación en C de algoritmos Monte Carlo para la representación dual del **modelo de Schwinger en 1+1 dimensiones**, desarrollada para un Trabajo de Fin de Grado en Física (Universidad de Zaragoza, curso 2025–2026).

La formulación dual utiliza ocupaciones de plaqueta, lazos fermiónicos, dímeros y, en el caso masivo, monómeros. El estudio considera la dependencia de los observables con la masa fermiónica, el acoplamiento gauge y el ángulo topológico. En el caso masivo aparecen pesos con signo, cuya contribución se recupera mediante reponderación en el análisis.

> **Software de investigación:** poder compilar el programa no constituye, por sí solo, una validación científica de los resultados Monte Carlo.

## Funcionalidad

- Actualizaciones locales de plaqueta y algoritmos de gusano abierto y cerrado.
- Medidas de la plaqueta, densidad de dímeros, respuesta topológica, densidad de monómeros y estimador del condensado quiral.
- Registro del signo de cada configuración y de los observables multiplicados por el signo.
- Series de termalización y barridos de parámetros.

**Plataforma compatible en esta versión:** Windows con Visual Studio 2022 y MSVC. El código emplea rutinas específicas de Windows (`rand_s`, `_mkdir` y `Sleep`). No se anuncia compatibilidad con Linux o macOS.

## Compilación

Instala Visual Studio 2022 con «Desarrollo para el escritorio con C++» (compilador MSVC y SDK de Windows) y CMake 3.20 o superior.

Abre una terminal de desarrollador en la raíz del repositorio y ejecuta:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

Ejecuta el programa desde la raíz para generar los resultados en `Results/`:

```powershell
.\build\Release\SchwingerWorm.exe
```

**Importante:** el `src/Main.c` actualmente publicado **no tiene menú interactivo**: comienza directamente el barrido configurado. Los parámetros por defecto producen una simulación extensa.

## Parámetros de la versión actual

| Parámetro | Valor predeterminado |
| --- | --- |
| Tamaño del retículo | `Nx = Nt = 10`, en `include/Lattice.h` |
| Acoplamiento | `beta = 10` |
| Masa | `m/g = 1` |
| Ángulo | `theta/pi = 0, 0.25, ..., 2` |
| Termalización | 100 000 pasos Monte Carlo por punto |
| Medidas | 50 176 medidas por punto |
| Separación entre medidas | 10 pasos Monte Carlo |

Para una prueba pequeña, modifica los parámetros en `src/Main.c` y recompila. Un paso Monte Carlo consta de un barrido de plaquetas, `max(1, V/10)` intentos de gusano y otro barrido de plaquetas.

La semilla aleatoria se obtiene del reloj del sistema; esta versión no ofrece una opción de semilla fija. No se garantiza que dos ejecuciones produzcan resultados idénticos bit a bit.

## Resultados

El programa crea `Results/` en el directorio desde el que se ejecuta:

- `Series_*.dat`: metadatos, paso Monte Carlo, diez observables, signo de configuración y diez observables reponderados con el signo.
- `Thermalization_*.txt`: evolución de la plaqueta media y la densidad de dímeros durante la termalización.

Estos archivos contienen medidas brutas. El análisis jackknife y las incertidumbres estadísticas requieren un posprocesado apropiado. Los datos `q` y `q^2` no constituyen necesariamente un estimador completo de la susceptibilidad como curvatura de la energía libre; puede ser necesario incluir una derivada explícita.

## Estructura

```text
.
├── CMakeLists.txt
├── include/Lattice.h
├── src/
│   ├── Main.c
│   ├── Initializer.c
│   ├── UpdCore.c
│   ├── Plaq_mod.c
│   ├── Worm_ferm.c
│   ├── Measures.c
│   ├── Visualizer.c
│   └── Other.c
└── README.md
```

## Memoria

**Worm Algorithms for the Schwinger Model**, Trabajo de Fin de Grado en Física, Universidad de Zaragoza (2025–2026). Los enlaces a Zaguán o arXiv se incorporarán cuando estén disponibles. La formulación dual procede de trabajos previos; el objetivo de este repositorio es documentar la implementación y el estudio numérico.
