


function build_heom_structure(K::Int, T::Int, hild::Int; useLtrunc::Bool=true)
  ados = Vector{Vector{Int}}()
  idxmap = Dict{Vector{Int},Int}()
  itvs = Vector{UnitRange{Int}}()

  counter = 0
  current = zeros(Int, 2K)

  function traverse(i::Int, sum::Int)
    if i > 2K
      r = copy(current)
      push!(ados, r)
      idxmap[r] = counter
      push!(itvs, counter*hild^2+1:(counter+1)*hild^2)
      counter += 1
      return
    end
    for val in 0:T-1
      if useLtrunc && (sum + val > T - 1)
        break
      end
      current[i] = val
      traverse(i + 1, sum + val)
    end
  end
  traverse(1, 0)

  return HEOMStructure(K, T, hild, useLtrunc, ados, idxmap, itvs)
end

function diag_ado(r, op::HEOMOperator, cops::CachedOps)
  st = op.structure

  ado_op = -im * cops.L_op
  ado_op += -sum(r[1:st.K] .* op.gamma[1:st.K] + r[st.K+1:2st.K] .* conj.(op.gamma)) * cops.Id_op

  return ado_op
end

function shift_ado(r, j::Int, offset::Int, direction::Int, op::HEOMOperator, cops::CachedOps)
  st = op.structure

  d_ext = vcat(op.d, conj.(op.d))

  if direction == -1
    return -im * sqrt(r[j+offset] * d_ext[j+offset]) * cops.Q_op
  elseif direction == +1
    if offset == 0
      return -im * sqrt((r[j+offset] + 1) * d_ext[j+offset]) * cops.Q_L_op
    elseif offset == st.K
      return +im * sqrt((r[j+offset] + 1) * d_ext[j+offset]) * cops.Q_R_op
    end
  end
end

function build_matrix(op::HEOMOperator)::SparseMatrixCSC
  st = op.structure

  offset_lst = [0, st.K]
  direction_lst = [-1, +1]

  dim = length(st.ados) * st.hild^2
  heom_matrix = spzeros(ComplexF64, dim, dim)

  cops = CachedOps(op.H, op.q)

  for i in eachindex(st.ados)
    r = st.ados[i]
    itv = st.itvs[i]

    ## Diag elements
    heom_matrix[itv, itv] = diag_ado(r, op, cops)

    ## Offdiagonal elements
    for j in 1:st.K
      for offset in offset_lst
        for direction in direction_lst
          r_shift = copy(r)
          r_shift[offset+j] += direction
          if haskey(st.idxmap, r_shift)
            itv_shift = st.itvs[st.idxmap[r_shift]+1]
            heom_matrix[itv_shift, itv] += shift_ado(r, j, offset, direction, op, cops)
          end
        end
      end
    end
  end

  return heom_matrix
end

