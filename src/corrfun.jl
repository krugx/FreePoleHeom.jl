
function bary_fit(
  beta::Float64,
  alpha::Float64,
  omega_c::Float64;
  rtol=1e-3,
  K_max=Int(100),
  npoints=Int(1e4),
  grid_fac=4.0,
  cutoff_fun=x -> 1 / (1 + x^2)^2
)
  J(omega) = alpha * omega * cutoff_fun(omega / omega_c)
  S(omega) = J(omega) * 1 / (1 - exp(-beta * omega))
  d = ComplexF64[]
  gamma = ComplexF64[]
  if alpha == 0.0
    return [0.0], [0.0], 1
  end

  omega = collect(
    range(-grid_fac * omega_c, grid_fac * omega_c, length=npoints)
  )

  Ga = aaa(
    omega,
    real(S.(omega)),
    tol=rtol,
    mmax=2 * K_max + 1,
    verbose=false
  )

  poles = prz(Ga)[1]
  residuals = prz(Ga)[2]
  for i in eachindex(poles)
    if imag(poles[i]) <= 0
      push!(d, -2 * im * residuals[i])
      push!(gamma, im * poles[i])
    end
  end
  return d, gamma, length(d)
end
