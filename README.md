# Monte Carlo Simulation of the Dual Schwinger Model

[Español](README_ES.md)

C implementation for simulating the dual representation of the Schwinger
model in 1+1 dimensions using plaquette updates and worm algorithms.

> Work in progress associated with a Bachelor's Thesis.

## Requirements

- Windows 10 or 11.
- Visual Studio 2022 with C support.
- CMake 3.20 or later.

## Project structure

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

## Build

Open a Visual Studio terminal in the project directory and run:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
```

The executable is usually generated at:

```text
build/Release/SchwingerWorm.exe
```

## Usage

Run:

```powershell
.\build\Release\SchwingerWorm.exe
```

At startup, the program displays a menu with the following options:

1. Quick test.
2. Theta sweep.
3. Mass sweep.
4. Mass-theta plane.
5. Single point.
6. Exit.

Each option shows the loaded parameters, number of simulation series,
thermalization steps, number of measurements, and spacing between
measurements.

The selected option must be confirmed with `y` or cancelled with `n`.

## Monte Carlo step

Each Monte Carlo step contains:

1. One complete sequential plaquette sweep.
2. `max(1, V / 10)` worm attempts.
3. A second complete sequential plaquette sweep.

## Output

Results are stored automatically in:

```text
Results/
```

Files named `Series_*.dat` contain:

- measured observables;
- the configuration sign;
- products between the sign and each observable for reweighting.

Files named `Thermalization_*.txt` contain the evolution of:

- the average plaquette;
- the dimer density.

## Configuration

General lattice parameters are defined in:

```text
include/Lattice.h
```

The simulation presets shown in the menu are defined in:

```text
src/Main.c
```

## Status

The code is currently under development. The full thesis and academic
documentation will be published after the evaluation process has been
completed.
