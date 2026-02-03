module FreePoleHeom

using LinearAlgebra
using SparseArrays
using BaryRational
using Arpack

include("types.jl")
include("heom.jl")
include("corrfun.jl")
include("stability.jl")
include("krylov.jl")

include("../ext/GradGenExt.jl")
using .GradGenExt

export
  HEOMStructure,
  HEOMOperator,
  HEOMPropagator,
  bary_fit,
  build_heom_structure,
  check_stability,
  Krylov,
  KrylovPropagator

end
