# FreePoleHeom.jl Setup & Troubleshooting

## Prerequisites

- **Julia**: Version 1.10 or higher.
- **OS**: Linux, macOS, or Windows (Linux recommended).

## Installation

1.  **Clone the repository**:
    ```bash
    git clone <repository_url>
    cd FreePoleHeom.jl
    ```

2.  **Instantiate the environment**:
    Launch Julia in the project directory:
    ```bash
    julia --project=.
    ```
    Enter the package manager by pressing `]`, then run:
    ```julia
    pkg> instantiate
    ```

## Running Tests

To verify the setup and current functionality, run the test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

This will run all tests defined in `test/runtests.jl`, including `test/qubit.jl` which demonstrates the GRAPE optimization usage.

## Troubleshooting: Zero Gradient with GRAPE + Krylov

### The Issue

When using `prop_method=:Krylov` with `method=GRAPE` in `QuantumControl.optimize`, the optimization may terminate early with a gradient norm of zero (`||∇J|| = 0.00e+00`).

**Root Cause Analysis:**
1.  **Backward Propagation Mismatch**: Standard GRAPE optimization relies on propagating a co-state (chi) backward in time.
2.  **GradVector Usage**: When `prop_method=:Krylov` is used, the generator is wrapped in a `GradGenerator`, causing the propagator to use `GradVector` (state + gradients).
3.  **Adjoint Physics**: The `Krylov` propagator (via `GradGenExt.jl`) implements the adjoint dynamics for `GradVector` as:
    $$ \dot{\phi} = H^\dagger \phi $$
    $$ \dot{g} = H^\dagger g + (\nabla H)^\dagger \phi $$
    In the backward pass for GRAPE, `phi` corresponds to the co-state $\chi$, and `g` corresponds to gradient accumulators. However, if the gradient accumulators `g` are initialized to zero (which they are), and the source term $(\nabla H)^\dagger \phi$ effectively couples into the state rather than accumulating a scalar overlap integral, the resulting gradient calculation becomes invalid or identically zero for standard functionals.
4.  **Forward vs. Backward**: `GradVector` is designed for **Forward Mode Differentiation (FMD)**, where gradients are propagated forward. Standard GRAPE uses the **Adjoint Method** (Backward). Using `GradVector` in a backward pass without a specific backward-compatible implementation (accumulating overlaps) leads to zero gradients.

### Solution Path 1: Enforce Forward Mode Differentiation (FMD)

The most direct fix, assuming `GradVector` logic is correct for forward propagation, is to configure the optimization to use Forward Mode Differentiation.

**Steps:**
1.  Modify the optimization call in your script (e.g., `test/qubit.jl`) to explicitly request forward gradients if the optimizer supports it.
    *   *Note*: `QuantumControl` / `GRAPE` might not expose a simple flag for this, requiring a change in how the problem is formulated (e.g., using `gradient_method=:forward` if supported, or using a different optimizer that defaults to FMD).

2.  Alternatively, ensure `GRAPE` detects that `Krylov` is an FMD-compatible propagator. This might involve defining specific traits for `Krylov` in `QuantumPropagators.Interfaces`.

### Solution Path 2: Implement Standard Adjoint Gradient Calculation

To support standard GRAPE (Backward Adjoint), the `KrylovPropagator` must handle standard `Vector{ComplexF64}` states (not `GradVector`) during the backward pass and compute the gradient via overlaps.

**Steps:**
1.  **Prevent `GradVector` Wrapping**: Ensure `GRAPE` passes standard vectors for the backward propagation. This usually happens automatically if the propagator doesn't claim to handle `GradVector`.
2.  **Manual Gradient Calculation**: If `GRAPE` relies on the propagator to calculate gradients, `KrylovPropagator` must expose an interface to compute $\langle \chi | \frac{\partial H}{\partial u} | \psi \rangle$ at each time step.
    *   This requires access to the **stored forward state** $\psi(t)$ during the backward pass.
    *   Currently, `KrylovPropagator` does not store the forward trajectory. You would need to implement `QuantumPropagators.Storage` or use `check=true` to store states, and then access them during the backward step.

### Summary
The current implementation of `evolve_step` for `GradVector` in `src/krylov.jl` enables the code to run (fixing the `MethodError`), but it does not mathematically support the backward gradient calculation expected by default GRAPE. Switching to a propagator that supports FMD (like `Cheby` seems to do) or re-implementing `Krylov` to support storage-based adjoint gradients is required for convergence.

## Debugging Notes (2026-02-03)

### Changes made for GRAPE + Krylov compatibility
- **`src/krylov.jl`**
  - Removed internal adjoint handling in `apply_H!` and `evolve_step` (backward propagation should use negative `dt` and an adjoint generator provided by GRAPE).
  - Implemented `QuantumPropagators.Interfaces.set_t!` and `set_state!` for proper reinitialization.
  - Fixed backward stepping logic in `prop_step!` (decrement `n` when `backward=true` and stop when out of valid range).
- **`ext/GradGenExt.jl`**
  - Aligned `GradVector` propagation with `GradGenerator` forward-mode semantics. Removed backward coupling that was mixing gradient states into the main state.
- **`test/qubit.jl`**
  - Added a Cheby vs Krylov gradient comparison block (relative diff ~5%, cosine similarity ~0.9998).
  - Experimented with tighter Krylov tolerances and Optim.jl optimizers.
- **Project deps / compatibility**
  - Added `Optim` to `Project.toml` and `Manifest.toml`.
  - Added `ext/OptimCompatExt.jl` to override GRAPE’s Optim integration for newer Optim (removed `ls_success` field).
  - Included this extension in `src/FreePoleHeom.jl`.

### Core issues encountered
- **GRAPE + Krylov (gradgen)**: Gradients are directionally correct but ~5% off in magnitude, enough to trigger LBFGSB line-search failure (`ABNORMAL_TERMINATION_IN_LNSRCH`).
- **GRAPE + Krylov (taylor)**: Terminates at 0 iterations (no progress).
- **Optim.jl integration**: GRAPE expects `ls_success` field (removed in Optim 1.13), requiring a compatibility override.
- **Early termination in Optim**: Even with `GradientDescent`, optimization stops after 1 iteration (no effective step).

### Suggested next steps
- **Isolate and validate gradient accuracy**:
  - Compare Krylov gradient vs finite differences for a handful of pulse parameters.
  - Increase Krylov accuracy further only if it changes gradient magnitude (not just direction).
- **Investigate GradVector propagation**:
  - Verify that the Krylov `GradVector` exponential is faithful to `GradGenerator` block structure.
  - Consider using a direct block-operator `LinearMap` to ensure exact block coupling.
- **Try GRAPE with `gradient_method=:taylor`**:
  - Tune `taylor_grad_max_order` and `taylor_grad_tolerance` to check if gradients stabilize.
- **Alternative route**:
  - Use a minimal custom GRAPE loop (forward/backward with overlap gradients) to validate Krylov propagation independently of GRAPE.jl’s abstractions.

### Revert plan (if desired)
- Remove `OptimCompatExt.jl` inclusion from `src/FreePoleHeom.jl` and delete `ext/OptimCompatExt.jl`.
- Remove `Optim` from `Project.toml`/`Manifest.toml`.
- Remove gradient comparison/debug blocks from `test/qubit.jl`.
