module WeirdMatrices

using LinearAlgebra
import LinearAlgebra: Matrix, AbstractMatrix, Array,
    Adjoint, Transpose,
    transpose, adjoint, diag, mul!, _modify!, MulAddMul
import Base: size, similar,
    +, -, *, /, ==,
    conj, copy, real, imag,
    require_one_based_indexing, getindex

using BandedMatrices
using SparseArrays

include("skew_tridiagonal.jl")
include("one_block_matrix.jl")
include("block_sparse_diagonal.jl")
include("column_sparse_matrix.jl")

end
