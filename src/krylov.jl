using QuantumPropagators
using QuantumPropagators.Controls
using QuantumPropagators.Generators: Generator
using QuantumPropagators.Interfaces: check_propagator
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
mutable struct KrylovPropagator{T, G, M} <: QuantumPropagators.AbstractPropagator
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

Applies the time-dependent Hamiltonian/Liouvillian `gen` at control values `u` to the state `x`, storing the result in `y`.

# Arguments
- `y`: Output vector.
- `x`: Input vector.
- `u`: Vector of control values.
- `gen`: Generator object.
- `backwards`: If true, applies the adjoint of the generator.
"""
function apply_H!(y, x, u, gen; backwards::Bool=false)
  A0 = gen.prop_0.mat
  # y = A0*x  
  if backwards
    mul!(y, adjoint(A0), x)
  else
    mul!(y, A0, x)
  end
  # y += Σ u[i] * A_i * x
  for i in eachindex(gen.prop_C)
    Ai = gen.prop_C[i].mat
    if backwards
      mul!(y, adjoint(Ai), x, u[i], 1.0)
    else
      mul!(y, Ai, x, u[i], 1.0)
    end
  end
  return y
end

"""
    apply_H!(y, x, u, gen::Generator; backwards::Bool=false)

Specialized application of the generator for `QuantumPropagators.Generators.Generator`.
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
        if backwards
            mul!(y, adjoint(H0), x)
        else
            mul!(y, H0, x)
        end
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
        if backwards
            mul!(y, adjoint(Hi), x, val, 1.0)
        else
            mul!(y, Hi, x, val, 1.0)
        end
    end
    return y
end

"""
    evolve_step!(phi_out, phi_in, u, dt, gen, p::Krylov; backwards=false)

In-place version of `evolve_step` for Krylov. Writes result to `phi_out`.

# Arguments
- `phi_out`: Output state vector.
- `phi_in`: Input state vector.
- `u`: Control values for the step.
- `dt`: Time step size.
- `gen`: Generator.
- `p`: Krylov method configuration.
- `backwards`: Propagation direction.
"""
function evolve_step!(phi_out::Vector, phi_in::Vector, u, dt, gen, p::Krylov; backwards=false)
  n = length(phi_in)
  Af = x -> begin
    y = similar(x)
    apply_H!(y, x, u, gen; backwards=false)
  end
  Ab = x -> begin
    y = similar(x)
    apply_H!(y, x, u, gen; backwards=true)
  end
  A = LinearMap{ComplexF64}(Af, Ab, n, n)

  # Use in-place expmv! which requires dense Vectors
  expmv!(phi_out, dt, (backwards ? adjoint(A) : A), phi_in; tol=p.tol, m=p.m, anorm=p.anorm)
  return phi_out
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
function evolve_step(phi::AbstractVector, u, dt, gen, p::Krylov; backwards=false)
  # Build a square LinearMap A: v ↦ (prop_0 + Σ u[i] prop_C[i]) v
  n = length(phi)

  Af = x -> begin
    y = similar(x)
    apply_H!(y, x, u, gen; backwards=false)
  end
  Ab = x -> begin
    y = similar(x)
    apply_H!(y, x, u, gen; backwards=true)
  end

  A = LinearMap{ComplexF64}(Af, Ab, n, n)

  v_dense = Vector(phi)
  w_dense = expmv(dt, (backwards ? adjoint(A) : A), v_dense; tol=p.tol, m=p.m, anorm=p.anorm)

  # Match output density to input density
  if phi isa SparseVector
    return sparse(w_dense)
  else # phi isa Vector
    return w_dense
  end
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
    inplace=true,
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
    if p.n >= length(p.tlist)
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
    
    if p.state isa Vector
        evolve_step!(p.state, p.state, vals, dt, p.gen, p.method; backwards=p.backward)
    else
        new_state = evolve_step(p.state, vals, dt, p.gen, p.method; backwards=p.backward)
        setfield!(p, :state, new_state)
    end
    
    setfield!(p, :t, t_end)
    setfield!(p, :n, p.n + 1)
    return p.state
end

# Optional: set_t!, set_state! (AbstractPropagator might handle these via fields, but good to be explicit)
function QuantumPropagators.Interfaces.set_t!(p::KrylovPropagator, t)
    setfield!(p, :t, t)
    # Update n: Find interval index
    # We assume t is in tlist. 
    # If t is not exactly in tlist, we find the interval it falls into.
    n = searchsortedlast(p.tlist, t)
    if n < 1
        n = 1
    elseif n >= length(p.tlist)
        n = length(p.tlist)
    end
    setfield!(p, :n, n)
end

function QuantumPropagators.Interfaces.set_state!(p::KrylovPropagator, state)
    if p.inplace && p.state isa AbstractArray && state isa AbstractArray && size(p.state) == size(state)
        copyto!(p.state, state)
    else
        setfield!(p, :state, state)
    end
end
