
"""
    check_stability(prop::HEOMPropagator)

Computes the steady state of the HEOM propagator and checks for stability.

Returns the physical steady state `rho_ss` and the full ADO vector if the smallest magnitude eigenvalue is close to 0 (stable), otherwise returns -1.
"""
function check_stability(prop::HEOMPropagator)
  st = prop.structure
  dim = length(st.ados) * st.hild^2

  v = zeros(ComplexF64, dim)
  # initial guess for steady state in convergence
  v[1:st.hild^2] = reshape(1 / 2 * [1 1-1im; 1+1im 1], (st.hild^2, 1))
  lambda, v_lambda = eigs(
    prop.mat,
    nev=1;
    which=:LR,
    maxiter=4000,
    check=1,
    v0=v
  )

  lambda = real(lambda[1])
  v_lambda = sparsevec(v_lambda)

  if lambda <= 1e-6
    rho_ss = reshape(v_lambda[1:st.hild^2], st.hild, st.hild)
    phys_SS = rho_ss / tr(rho_ss)
    ado_ss = v_lambda / tr(rho_ss)

    return phys_SS, ado_ss
  else
    return -1
  end
end
