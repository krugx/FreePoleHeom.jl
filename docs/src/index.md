# FreePoleHeom.jl

`FreePoleHeom.jl` is a Julia package for simulating Open Quantum Systems using the Free-Pole Hierarchy of Equations of Motion (FP-HEOM). This method allows for efficient simulation of system dynamics in structured environments by decomposing the bath correlation functions into a sum of exponentials (poles).

## Key Features

- **Efficient Decomposition**: Uses barycentric rational fitting (`bary_fit`) to optimally decompose spectral densities.
- **Sparse Algebra**: Utilizes sparse matrices for efficient memory usage and computation of the HEOM generator.
- **Stability Analysis**: Includes tools to check the convergence and stability of the hierarchy (`check_stability`).

## Usage Contexts

This package is designed to be versatile, supporting both:

1.  **Lab Frame**: Standard HEOM where the bath correlation function is decomposed directly.
2.  **Rotating Frame**: Extended HEOM formulation (FP-HEOM) that accounts for rotating frame transformations, splitting the correlation function into positive and negative frequency branches (see [Theory](theory.md)).

## Documentation

- [Theory](theory.md): Mathematical derivation and details on the Rotating Frame FP-HEOM.
- [API Reference](api.md): Documentation of types and functions.

## Installation

```julia
using Pkg
Pkg.add("FreePoleHeom")
```