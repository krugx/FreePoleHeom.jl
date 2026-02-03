using FreePoleHeom
using Test
using LinearAlgebra
using SparseArrays

@testset "FreePoleHeom.jl" begin

  @testset "Barycentric Fit" begin
    beta = 1.0
    alpha = 0.1
    omega_c = 10.0
    d, gamma, K = bary_fit(beta, alpha, omega_c; rtol=1e-2)

    @test K > 0
    @test length(d) == K
    @test length(gamma) == K
    @test eltype(d) == ComplexF64
    @test eltype(gamma) == ComplexF64
  end

  @testset "HEOM Structure" begin
    K = 1
    T = 2
    hild = 2
    st = build_heom_structure(K, T, hild; useLtrunc=true)

    @test st.K == K
    @test st.T == T
    @test st.hild == hild
    @test !isempty(st.ados)
    @test length(st.itvs) == length(st.ados)
  end

  @testset "HEOM Propagator Construction" begin
    K = 1
    T = 2
    hild = 2
    st = build_heom_structure(K, T, hild)

    # Simple 2-level system
    H = [0.0+0im 0.0; 0.0 1.0]
    q = [[0.0+0im 1.0; 1.0 0.0]] # One coupling operator

    d = [0.1 + 0.0im]
    gamma = [0.5 + 0.0im]

    # q needs to be replicated for each pole if we treat them as independent modes in this specific API structure
    # Actually, heom.jl says: @assert length(q) == structure.K
    # If K=1, q should have length 1.

    prop = HEOMPropagator(st, H, q, d, gamma)

    @test isa(prop.mat, SparseMatrixCSC)
    dim = length(st.ados) * hild^2
    @test size(prop.mat) == (dim, dim)
  end

end

include("test_krylov.jl")
