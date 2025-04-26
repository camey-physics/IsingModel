# 3D Ising Model Simulation

This project implements a three-dimensional (3D) Ising model with periodic boundary conditions and multiple Monte Carlo update methods. It includes unit tests and a basic simulated annealing script that measures the Binder cumulant as a function of temperature.

## Features

- **Ising Model**:
  - 3D cubic lattice with periodic boundary conditions
  - Spin values ±1
- **Monte Carlo Update Methods**:
  - Metropolis
  - Heat Bath
  - Wolff Cluster Algorithm
- **Observable Calculations**:
  - Energy per spin
  - Magnetization per spin
  - Binder cumulant (with bootstrap resampling for error estimation)
- **Example Simulation**:
  - Simulated annealing sweep across a temperature range
  - Measurement of Binder cumulant and statistical uncertainty
- **Unit Tests**:
  - Lattice initialization
  - Periodic boundary conditions
  - Neighbor table validation
  - Energy and magnetization calculations
  - Monte Carlo updates
  - Copying model state

## Requirements

- C++17 or later
- [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/)
- [GoogleTest](https://github.com/google/googletest) (for unit testing)

## Building

This project uses CMake for configuration.

Build steps:

```bash
mkdir build
cd build
cmake ..
make

This will generate:

- `simulatedAnnealing` — the basic simulated annealing program
- `IsingModelTest` — the unit tests executable
