
struct HEOMStructure
  K::Int
  T::Int
  hild::Int
  useLtrunc::Bool
  ados::Vector{Vector{Int}}        # multi-index ADOs
  idxmap::Dict{Vector{Int},Int}   # multi-index → flat index
  itvs::Vector{UnitRange{Int}}
end

struct HEOMOperator
  structure::HEOMStructure
  H::Matrix{ComplexF64}
  q::Matrix{ComplexF64}
  d::Vector{ComplexF64}
  gamma::Vector{ComplexF64}
end

struct CachedOps
  Id_op::SparseMatrixCSC
  L_op::SparseMatrixCSC
  Q_op::SparseMatrixCSC
  Q_L_op::SparseMatrixCSC
  Q_R_op::SparseMatrixCSC
  function CachedOps(H, q)
    id = sparse(diagm(ones(size(H, 1))))
    Id_op = kron(transpose(id), id)
    L_op = kron(transpose(id), H) - kron(transpose(H), id)
    Id_op = kron(transpose(id), id)
    Q_op = kron(transpose(id), q) - kron(transpose(q), id)
    Q_L_op = kron(transpose(id), q)
    Q_R_op = kron(transpose(q), id)

    new(Id_op, L_op, Q_op, Q_L_op, Q_R_op)
  end
end

struct HEOMPropagator
  mat::SparseMatrixCSC
  structure::HEOMStructure
  op::HEOMOperator
  function HEOMPropagator(op::HEOMOperator)
    new(build_matrix(op), op.structure, op)
  end
  function HEOMPropagator(structure, H, q, d, gamma)
    @assert size(H) == size(q)
    @assert length(d) == length(gamma)
    op = HEOMOperator(structure, H, q, d, gamma)
    new(build_matrix(op), structure, op)
  end
end
