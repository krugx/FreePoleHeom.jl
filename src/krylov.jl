using QuantumPropagators
using QuantumPropagators.Controls
using QuantumPropagators.Generators: Generator
using QuantumPropagators.Interfaces: check_propagator
using QuantumGradientGenerators
using LinearAlgebra
using SparseArrays
using Expokit
using LinearMaps

"""
    Krylov(; m=30, tol=1e-7, anorm=1.0)

Configuration for the Krylov subspace propagation method.

# Arguments
- `m`: Dimension of the Krylov subspace.
- `tol`: Tolerance for the `expmv` method.
- `anorm`: Norm of the operator, if known.
"""
struct Krylov
  m::Int
  tol::Float64
  anorm::Float64
end

Krylov(; m=30, tol=1e-7, anorm=1.0) = Krylov(m, tol, anorm)

"""
    KrylovPropagator

Propagator for the Krylov subspace method.

# Fields
- `state`: Current state vector.
- `gen`: Generator of the dynamics.
- `t`: Current time.
- `n`: Current time step index.
- `tlist`: List of time points.
- `parameters`: Parameters for the propagator.
- `backward`: If true, propagates backward in time.
- `inplace`: If true, modifies `state` in-place.
- `method`: The `Krylov` configuration.
"""
mutable struct KrylovPropagator{T,G,M} <: QuantumPropagators.AbstractPropagator
  state::T
  gen::G
  t::Float64
  n::Int
  tlist::Vector{Float64}
  parameters::AbstractDict
  backward::Bool
  inplace::Bool
  method::M
end

"""
    apply_H!(y, x, u, gen; backwards::Bool=false)

Applies the time-dependent Hamiltonian/Liouvillian `gen` at control values `u`
to the state `x`, storing the result in `y`.

# Arguments
- `y`: Output vector.
- `x`: Input vector.
- `u`: Vector of control values.
- `gen`: Generator object.
- `backwards`: Unused. Backward propagation is handled via negative `dt` and
  an adjoint generator provided by the caller.
"""
function apply_H!(y, x, u, gen; backwards::Bool=false)
  A0 = gen.prop_0.mat
  # y = A0*x
  mul!(y, A0, x)
  # y += Σ u[i] * A_i * x
  for i in eachindex(gen.prop_C)
    Ai = gen.prop_C[i].mat
    mul!(y, Ai, x, u[i], 1.0)
  end
  return y
end

"""
    apply_H!(y, x, u, gen::Generator; backwards::Bool=false)

Specialized application of the generator for `QuantumPropagators.Generators.Generator`.
The `backwards` flag is unused; backward propagation is handled via negative
`dt` and an adjoint generator.
"""
function apply_H!(y, x, u, gen::Generator; backwards::Bool=false)
  ops = gen.ops
  amps = gen.amplitudes
  n_ops = length(ops)
  n_amps = length(amps)

  op_offset = 0

  # Handle drift term
  if n_ops == n_amps + 1
    # ops[1] is drift
    H0 = ops[1]
    mul!(y, H0, x)
    op_offset = 1
  elseif n_ops == n_amps
    # No drift, initialize y to 0
    fill!(y, zero(eltype(y)))
    op_offset = 0
  else
    error("Mismatch in ops and amplitudes")
  end

  # Apply control terms
  for i in 1:n_amps
    Hi = ops[i+op_offset]
    val = u[i]
    mul!(y, Hi, x, val, 1.0)
  end
  return y
end

"""
    evolve_step(phi, u, dt, gen, p::Krylov; backwards=false)

Evolves the state `phi` by one time step `dt` using the Krylov subspace method.

# Arguments
- `phi`: Input state vector.
- `u`: Control values.
- `dt`: Time step.
- `gen`: Generator.
- `p`: Krylov method configuration.
- `backwards`: Propagation direction.
"""
function evolve_step(phi, u, dt, gen, p::Krylov; backwards=false)
  # Build a square LinearMap A: v ↦ (prop_0 + Σ u[i] prop_C[i]) v
  n = length(phi)

  Af = x -> begin
    y = similar(x)
    apply_H!(y, x, u, gen; backwards=false)
  end
  Ab = Af

  A = LinearMap{ComplexF64}(Af, Ab, n, n)

  v_dense = phi
  w_dense = expmv(-1im * dt, A, v_dense; tol=p.tol, m=p.m, anorm=p.anorm)

  return w_dense
end

function evolve_step(phi::GradVector, u, dt, gen, p::Krylov; backwards=false)
  # Convert GradVector to dense vector
  v_flat = flatten(phi)
  n = length(v_flat)

  # Pre-allocate auxiliary structures
  aux_x = deepcopy(phi)
  aux_y = deepcopy(phi)

  Af! = (y_flat, x_flat) -> begin
    copyto!(aux_x, x_flat)
    apply_H!(aux_y, aux_x, u, gen; backwards=false)
    copyto!(y_flat, aux_y)
    return y_flat
  end

  Ab! = Af!

  A = LinearMap{ComplexF64}(Af!, Ab!, n, n; ismutating=true)

  w_flat = expmv(-1im * dt, A, v_flat; tol=p.tol, m=p.m, anorm=p.anorm)

  w = deepcopy(phi)
  copyto!(w, w_flat)
  return w
end

"""
    init_prop(state, gen, tlist, method::Val{:Krylov}; ...)

Initializes the Krylov propagator.

# Arguments
- `state`: Initial state.
- `gen`: Generator.
- `tlist`: Time grid.
- `method`: `Val(:Krylov)`.
- `inplace`: Whether to use in-place operations (default: `true`).
- `backward`: Whether to propagate backward (default: `false`).
- `m`: Krylov subspace dimension (default: 30).
- `tol`: Tolerance (default: 1e-7).
- `anorm`: Operator norm estimate (default: 1.0).
"""
function QuantumPropagators.init_prop(
  state,
  gen,
  tlist,
  method::Val{:Krylov};
  inplace=false,
  backward=false,
  verbose=false,
  parameters=nothing,
  m=30,
  tol=1e-7,
  anorm=1.0,
  _...
)
  t = tlist[1]
  parameters = isnothing(parameters) ? Dict() : parameters
  k_method = Krylov(m, tol, anorm)

  # Check consistency of state
  state_curr = inplace ? copy(state) : state

  return KrylovPropagator(
    state_curr,
    gen,
    t,
    1, # n (step index)
    tlist,
    parameters,
    backward,
    inplace,
    k_method
  )
end

# Alias for module name usage (if FreePoleHeom is used as method)
QuantumPropagators.init_prop(state, gen, tlist, method::Val{:FreePoleHeom}; kwargs...) =
  QuantumPropagators.init_prop(state, gen, tlist, Val(:Krylov); kwargs...)

"""
    prop_step!(p::KrylovPropagator)

Performs a single propagation step using the Krylov method.
"""
function QuantumPropagators.prop_step!(p::KrylovPropagator)
  N_T = length(p.tlist) - 1
  if p.n < 1 || p.n > N_T
    return nothing
  end

  t_start = p.tlist[p.n]
  t_end = p.tlist[p.n+1]
  dt = t_end - t_start
  if p.backward
    dt = -dt
  end

  t_mid = t_start + (t_end - t_start) / 2

  # Get controls and evaluate them
  controls = QuantumPropagators.Controls.get_controls(p.gen)
  vals = [QuantumPropagators.Controls.evaluate(c, p.tlist, p.n) for c in controls]

  new_state = evolve_step(p.state, vals, dt, p.gen, p.method; backwards=p.backward)
  setfield!(p, :state, new_state)

  setfield!(p, :t, p.backward ? t_start : t_end)
  setfield!(p, :n, p.n + (p.backward ? -1 : 1))
  return p.state
end

function flatten(v::AbstractVector)
  return Vector(v)
end


# set_t!, set_state! to support reinit_prop! and GRAPE workflows
function QuantumPropagators.Interfaces.set_t!(p::KrylovPropagator, t)
  setfield!(p, :t, t)
  # Update n: index of the interval starting at t
  n = searchsortedlast(p.tlist, t)
  if n < 1
    n = 1
  elseif n >= length(p.tlist)
    n = length(p.tlist) - 1
  end
  setfield!(p, :n, n)
end

function QuantumPropagators.Interfaces.set_state!(p::KrylovPropagator, state)
  if p.inplace && p.state isa AbstractArray && state isa AbstractArray && size(p.state) == size(state)
    copyto!(p.state, state)
    return p.state
  else
    setfield!(p, :state, state)
    return p.state
  end
end
