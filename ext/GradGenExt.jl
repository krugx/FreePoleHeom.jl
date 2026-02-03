module GradGenExt

using FreePoleHeom
using QuantumGradientGenerators
using LinearAlgebra
import Base: copyto!

function FreePoleHeom.flatten(v::GradVector)
  # Concatenate state and all grad_states
  return vcat(v.state, v.grad_states...)
end

function copyto!(dest::GradVector, src::AbstractVector)
  n_state = length(dest.state)
  n_grads = length(dest.grad_states)

  copyto!(dest.state, view(src, 1:n_state))

  offset = n_state
  for i in 1:n_grads
    copyto!(dest.grad_states[i], view(src, offset+1:offset+n_state))
    offset += n_state
  end
  return dest
end

function copyto!(dest::AbstractVector, src::GradVector)
  n_state = length(src.state)
  n_grads = length(src.grad_states)

  copyto!(view(dest, 1:n_state), src.state)

  offset = n_state
  for i in 1:n_grads
    copyto!(view(dest, offset+1:offset+n_state), src.grad_states[i])
    offset += n_state
  end
  return dest
end

function FreePoleHeom.apply_H!(y::GradVector, x::GradVector, u, gen::GradGenerator; backwards::Bool=false)
  # Diagonal parts (H on state, H on grad_states)
  FreePoleHeom.apply_H!(y.state, x.state, u, gen.G; backwards=backwards)
  for i in eachindex(y.grad_states)
    FreePoleHeom.apply_H!(y.grad_states[i], x.grad_states[i], u, gen.G; backwards=backwards)
  end

  # Off-diagonal parts
  # Forward:  d(grad)/dt = ... + H_i * state
  # Backward: d(state)/dt = ... + H_i' * grad
  for i in eachindex(gen.control_derivs)
    Hi = gen.control_derivs[i]
    if backwards
      mul!(y.state, adjoint(Hi), x.grad_states[i], 1.0, 1.0)
    else
      mul!(y.grad_states[i], Hi, x.state, 1.0, 1.0)
    end
  end
  return y
end

end # module
