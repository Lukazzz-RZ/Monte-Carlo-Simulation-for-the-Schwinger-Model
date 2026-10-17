# Monte Carlo Simulation of the Schwinger Model

[Versión en español](README_ES.md) · [Source code](src/) · [Build and smoke test](.github/workflows/windows-build.yml)

**C implementation of Monte Carlo algorithms for the (1+1)-dimensional Schwinger model in its dual monomer–dimer–loop formulation.** Developed for a Bachelor's thesis in Physics at the University of Zaragoza (2025–2026), this project explores gauge-field and fermionic observables as functions of the mass, gauge coupling and topological angle.

## Representative result

![Heatmap of the topological response versus topological angle and fermion mass](docs/figures/topological_response_heatmap.svg)

*Topological-response estimator ⟨q⟩ as a function of θ/π and m/g, at β = 14.4 on a 12 × 12 lattice. Each coloured cell represents a sampled angle for one of four discrete masses; **no interpolation between masses** is implied. This overview is based on [rounded thesis data](docs/data/topological_response_preview.csv). Error bars are not shown here; the thesis should be consulted for the statistical analysis and physical interpretation. This is not the output of the repository's current default parameter sweep.*

## Key features

- **Dual lattice representation:** gauge plaquette occupation numbers, fermionic loops, dimers and monomers.
- **Monte Carlo updates:** local plaquette moves and open/closed worm algorithms.
- **Measurements:** plaquette, dimer and monomer densities, topological-response estimator, and a chiral-condensate estimator.
- **Sign reweighting:** configuration signs and sign-weighted observables saved for subsequent statistical analysis.
- **Build verification:** CMake build plus an automated Windows smoke test in GitHub Actions.

## Build and run

**Supported environment:** Windows, Visual Studio 2022 with *Desktop development with C++* (MSVC and Windows SDK), and CMake ≥ 3.20. The current source code uses Windows-specific functions such as `rand_s`, `_mkdir` and `Sleep`; Linux and macOS are not supported.

Open a Visual Studio developer terminal at the repository root:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
.\build\Release\SchwingerWorm.exe
```

**Note:** the current executable immediately runs the parameter sweep configured in `src/Main.c`; it does **not** display an interactive menu. The default workload is substantial.

### Current default parameters

| Setting | Value |
| --- | --- |
| Lattice size | `Nx = Nt = 10` |
| Gauge coupling | `β = 10` |
| Fermion mass | `m/g = 1` |
| Topological angles | `θ/π = 0, 0.25, …, 2` |
| Thermalisation | 100,000 Monte Carlo steps per point |
| Measurements | 50,176 per point |
| Steps between measurements | 10 |

The simulation settings and scan arrays live in [`src/Main.c`](src/Main.c); lattice dimensions are in [`include/Lattice.h`](include/Lattice.h). To run a shorter study, change those constants and rebuild. Each Monte Carlo step consists of one plaquette sweep, `max(1, V/10)` worm attempts, and another plaquette sweep.

## Simulation output

Running from the repository root creates a `Results/` directory containing:

| File | Contents |
| --- | --- |
| `Series_*.dat` | Parameters, Monte Carlo step, ten observables, configuration sign, and ten sign-weighted observables |
| `Thermalization_*.txt` | Thermalisation history of the mean plaquette and dimer density |

These are **raw measurements**, not final statistical estimates. In particular, `q` and `q²` alone must not be identified with the full topological susceptibility defined by the free-energy curvature; the corresponding explicit derivative contribution must also be accounted for.

## Validation and reproducibility

The [Windows CI workflow](.github/workflows/windows-build.yml) checks that the project builds with MSVC and that a **shortened simulation** runs on a disposable copy of the source, produces finite measurements and writes files with the expected structure. See [`scripts/smoke-test.ps1`](scripts/smoke-test.ps1).

This is an **engineering smoke test**, not a proof of ergodicity or of the physical correctness of the numerical results. The production executable currently seeds its pseudorandom generator from the clock and does not offer a user-selectable fixed seed.

## Project structure

```text
include/       Lattice declarations and parameters
src/           C source code: updates, measurements, initialisation, output
docs/data/     Rounded data for the README figure
docs/figures/  Scientific overview figures
scripts/       Isolated smoke test
CMakeLists.txt Windows/MSVC build
```

## Academic context

This repository accompanies **_Worm Algorithms for the Schwinger Model_**, a Bachelor's thesis in Physics at the University of Zaragoza (2025–2026). The dual formulation builds on previous lattice-field-theory literature; the purpose of this repository is to share the implementation and numerical study.

A public link to the complete thesis in the university repository or arXiv can be added after deposition.
