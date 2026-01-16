
"""
    HEOMStructure

Data structure defining the hierarchy layout.
# Fields
- `K`: Number of bath modes (or pairs of poles in RWA).
- `T`: Hierarchy depth (truncation level).
- `hild`: Hilbert space dimension of the system.
- `useLtrunc`: If true, uses L-truncation (sum n_k <= T).
- `ados`: List of multi-indices.
- `idxmap`: Map from multi-index to linear index.
- `itvs`: UnitRanges for indexing the large matrix.
"""
struct HEOMStructure
  K::Int
  T::Int
  hild::Int
  useLtrunc::Bool
  ados::Vector{Vector{Int}}        # multi-index ADOs
  idxmap::Dict{Vector{Int},Int}   # multi-index → flat index
  itvs::Vector{UnitRange{Int}}
end

"""
    HEOMOperator

Contains the physical operators and bath parameters required to build the HEOM generator.
# Fields
- `structure`: The `HEOMStructure`.
- `H`: System Hamiltonian.
- `q`: Coupling operators (vector of matrices).
- `d`: Expansion coefficients (residues).
- `gamma`: Expansion frequencies (poles).
"""
struct HEOMOperator
  structure::HEOMStructure
  H::Matrix{ComplexF64}
  q::Vector{Matrix{ComplexF64}}
  d::Vector{ComplexF64}
  gamma::Vector{ComplexF64}
end

struct CachedOps
  Id_op::SparseMatrixCSC
  L_op::SparseMatrixCSC
  Q_op::Vector{SparseMatrixCSC}
  Q_L_op::Vector{SparseMatrixCSC}
  Q_R_op::Vector{SparseMatrixCSC}
  function CachedOps(H, q)
    id = sparse(diagm(ones(size(H, 1))))
    Id_op = kron(transpose(id), id)
    L_op = kron(transpose(id), H) - kron(transpose(H), id)
    Id_op = kron(transpose(id), id)
    Q_op = []
    Q_L_op = []
    Q_R_op = []
    for (i, qi) in enumerate(q)
      push!(Q_op, kron(transpose(id), qi) - kron(transpose(qi), id))
      push!(Q_L_op, kron(transpose(id), qi))
      push!(Q_R_op, kron(transpose(qi), id))
    end

    new(Id_op, L_op, Q_op, Q_L_op, Q_R_op)
  end
end

"""
    HEOMPropagator

Wraps the sparse Liouvillian matrix for the HEOM evolution.
# Fields
- `mat`: The sparse HEOM generator matrix.
- `structure`: The `HEOMStructure`.
- `op`: The `HEOMOperator`.
"""
struct HEOMPropagator
  mat::SparseMatrixCSC
  structure::HEOMStructure
  op::HEOMOperator
  function HEOMPropagator(op::HEOMOperator)
    new(build_matrix(op), op.structure, op)
  end
  function HEOMPropagator(structure, H, q, d, gamma)
    @assert size(H) == size(q[1])
    @assert length(q) == structure.K
    @assert length(d) == length(gamma)
    op = HEOMOperator(structure, H, q, d, gamma)
    new(build_matrix(op), structure, op)
  end
end
