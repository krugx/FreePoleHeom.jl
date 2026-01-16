# FreePoleHeom

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://krugx.github.io/FreePoleHeom.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://krugx.github.io/FreePoleHeom.jl/dev/)
[![Build Status](https://github.com/krugx/FreePoleHeom.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/krugx/FreePoleHeom.jl/actions/workflows/CI.yml?query=branch%3Amain)

`FreePoleHeom.jl` is a high-performance Julia package for simulating Open Quantum Systems using the **Free-Pole Hierarchy of Equations of Motion (FP-HEOM)**.

The package provides tools for:
- **Optimal Bath Decomposition**: Efficiently fitting spectral densities using barycentric rational approximations.
- **Versatile Dynamics**: Support for both standard **Lab Frame** HEOM and an extended formulation for the **Rotating Frame** (RWA-HEOM).
- **Scalable Computation**: Leveraging sparse matrices and optimized hierarchy structures for large-scale simulations.

For detailed theoretical background and API usage, please refer to the [Documentation](https://krugx.github.io/FreePoleHeom.jl/dev/).
