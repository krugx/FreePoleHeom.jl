# FreePoleHeom.jl

`FreePoleHeom.jl` is a Julia package for simulating Open Quantum Systems using the Free-Pole Hierarchy of Equations of Motion (FP-HEOM). This method allows for efficient simulation of system dynamics in structured environments by decomposing the bath correlation functions into a sum of exponentials (poles).

## Key Features

- **Efficient Decomposition**: Uses barycentric rational fitting (`bary_fit`) to optimally decompose spectral densities.
- **Sparse Algebra**: Utilizes sparse matrices for efficient memory usage and computation of the HEOM generator.
- **Stability Analysis**: Includes tools to check the convergence and stability of the hierarchy (`check_stability`).


## Documentation

- [API Reference](api.md): Documentation of types and functions.
