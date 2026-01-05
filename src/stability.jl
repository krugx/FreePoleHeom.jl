
function check_stability(prop::HEOMPropagator)
  st = prop.structure
  dim = length(st.ados) * st.hild^2

  v = zeros(ComplexF64, dim)
  # initial guess for steady state in convergence
  v[1:st.hild^2] = reshape(1 / 2 * [1 1-im; 1+im 1], (st.hild^2, 1))
  lambda, v_lambda = eigs(prop.mat, nev=1; which=:LR, maxiter=2000, check=1, v0=v)

  lambda = real(lambda[1])
  v_lambda = sparsevec(v_lambda)

  ## Check if stable
  println("max({Re(lambda_i)}) = $lambda \n")
  if lambda <= 1e-6
    rho_SS = reshape(v_lambda[1:st.hild^2], st.hild, st.hild)
    phys_SS = rho_SS / tr(rho_SS)
    ado_SS = v_lambda / tr(rho_SS)

    println("Physical steady state rho=")
    display(rho_SS)
    println("\n")
    println("Steady state of ADO state = ")
    display(v_lambda)
    return phys_SS, ado_SS
  else
    println("Instable HEOM dynamics!")
    return -1
  end
end
