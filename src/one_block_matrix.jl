struct OneBlockMatrix{T,Block<:AbstractMatrix{T}} <: AbstractMatrix{T}
    block::Block
    m::Int
    n::Int
    function OneBlockMatrix(block::Block, m, n) where {T,Block<:AbstractMatrix{T}}
        require_one_based_indexing(block)
        mb,nb = size(block)
        mb ≤ m && nb ≤ n ||
            throw(DimensionMismatch("Size of subblock ($mb,$nb) must be smaller than or equal to size of matrix ($m,$n)"))
        new{T,Block}(block, m, n)
    end
end

Base.size(A::OneBlockMatrix) = (A.m, A.n)
Base.size(A::OneBlockMatrix, i) = size(A)[i]

Base.similar(A::OneBlockMatrix, ::Type{T}, m::Integer, n::Integer) where T =
    OneBlockMatrix(similar(A.block, T), m, n)

Base.similar(A::OneBlockMatrix, ::Type{T}, (m,n)::Tuple{<:Integer,<:Tuple}) where T =
    similar(A, T, m, n)

Base.similar(A::OneBlockMatrix, ::Type{T}) where T =
    OneBlockMatrix(similar(A.block, T), A.m, A.n)

# Apparently ambiguous
Base.similar(A::OneBlockMatrix{T}, dims::Union{Integer, AbstractUnitRange}...) where T=
    similar(A, T, dims...)

function Base.getindex(A::OneBlockMatrix{T}, i, j) where T
    m,n = size(A.block)
    i ≤ m && j ≤ n && return A.block[i,j]
    zero(T)
end

function Base.replace_in_print_matrix(A::OneBlockMatrix,i::Integer,j::Integer,s::AbstractString)
    m,n = size(A.block)
    i ≤ m && j ≤ n && return s
    Base.replace_with_centered_mark(s)
end

for T in [Integer,Real]
    @eval begin
        LinearAlgebra.:(^)(A::OneBlockMatrix, p::$T) =
            OneBlockMatrix(A.block^p, A.m, A.n)
    end
end

Base.:(*)(α::Number, A::OneBlockMatrix) =
    OneBlockMatrix(α*A.block, A.m, A.n)

Base.:(*)(A::OneBlockMatrix, α::Number) =
    OneBlockMatrix(α*A.block, A.m, A.n)

Base.:(/)(A::OneBlockMatrix, α::Number) = A * inv(α)
Base.:(//)(A::OneBlockMatrix, α::Number) = A * (1//α)

for f in [:exp, :sin, :cos, :tan, :sqrt]
    @eval begin
        function Base.$f(A::OneBlockMatrix)
            M,N = size(A)
            M == N || throw(DimensionMismatch("matrix is not square: dimensions are ($(M), $(N))"))

            b = A.block
            m,n = size(b)

            block = if m == n
                $f(b)
            else
                ee = eigen(A)
                Q = first(ee.vectors.blocks)
                Q*Diagonal($f.(ee.values[1:m]))*inv(Q)
            end

            fz = $f(zero(eltype(block)))
            D = Diagonal(fill(fz, M))

            BlockSparseDiagonal(D, [block-fz*I], [1])
        end
    end
end

# * Matrix products

# ** Arbitrary matrices

function LinearAlgebra.mul!(y::AbstractMatrix, A::OneBlockMatrix, x::AbstractMatrix,
                            α::Number=true, β::Number=false)
    mY,nY = size(y)
    mA = size(A,1)
    nX = size(x,2)
    mY == mA && nY == nX ||
        throw(DimensionMismatch("Dimensions of y ($(mY),$(nY)) must agree with first dimension of A ($(mA)), and second dimension of x ($(nX)), respectively"))

    if iszero(β)
        y .= false
    elseif !isone(β)
        lmul!(β, y)
    end

    iszero(α) && return y

    m,n = size(A.block)
    mul!(view(y, 1:m, :), A.block, view(x, 1:n, :),
         α, true)

    y
end

# ** OneBlockMatrices

function LinearAlgebra.mul!(y::OneBlockMatrix, A::OneBlockMatrix, x::OneBlockMatrix,
                            α::Number=true, β::Number=false)
    mY,nY = size(y)
    mA = size(A,1)
    nX = size(x,2)
    mY == mA && nY == nX ||
        throw(DimensionMismatch("Dimensions of y ($(mY),$(nY)) must agree with first dimension of A ($(mA)), and second dimension of x ($(nX)), respectively"))

    mY,nY = size(y.block)
    mA,nA = size(A.block)
    mX,nX = size(x.block)

    mY == mA && nY == nX ||
        throw(DimensionMismatch("Dimensions of y subblock ($(mY),$(nY)) must agree with first dimension of A subblock ($(mA)), and second dimension of x subblock ($(nX)), respectively"))

    if iszero(β)
        y.block .= false
    elseif !isone(β)
        lmul!(β, y.block)
    end

    iszero(α) && return y

    interval = 1:min(nA,mX)
    mul!(y.block, view(A.block, :, interval), view(x.block, interval, :), α, true)

    y
end

function *(A::OneBlockMatrix, B::OneBlockMatrix)
    T = Base.promote_op(*, eltype(A), eltype(B))
    m = size(A.block,1)
    n = size(B.block,2)
    block = Matrix{T}(undef, m, n)
    C = OneBlockMatrix(block, size(A, 1), size(B, 2))
    mul!(C, A, B)
end

# ** Diagonal matrices

function mul!(y::OneBlockMatrix, A::OneBlockMatrix, D::Diagonal,
              α::Number=true, β::Number=false)
    mY,nY = size(y)
    mA = size(A,1)
    nD = size(D,2)
    mY == mA && nY == nD ||
        throw(DimensionMismatch("Dimensions of y ($(mY),$(nY)) must agree with first dimension of A ($(mA)), and second dimension of D ($(nD)), respectively"))

    m = size(A.block,1)
    mul!(y.block, A.block, Diagonal(view(D.diag, 1:m)), α, β)
    y
end

*(A::OneBlockMatrix, D::Diagonal) =
    mul!(similar(A), A, D)

function mul!(y::OneBlockMatrix, D::Diagonal, A::OneBlockMatrix,
              α::Number=true, β::Number=false)
    mY,nY = size(y)
    mD = size(D,1)
    nA = size(A,2)
    mY == mD && nY == nA ||
        throw(DimensionMismatch("Dimensions of y ($(mY),$(nY)) must agree with first dimension of D ($(mD)), and second dimension of A ($(nA)), respectively"))

    m = size(A.block,1)
    mul!(y.block, Diagonal(view(D.diag, 1:m)), A.block, α, β)
    y
end

*(D::Diagonal, A::OneBlockMatrix) =
    mul!(similar(A), D, A)

# ** Various banded matrices

const VariousBandedMatrices = Union{Tridiagonal,SymTridiagonal,Bidiagonal,BandedMatrix,SkewTridiagonal}

function _mul!(y::OneBlockMatrix, T::VariousBandedMatrices, A::OneBlockMatrix,
               α::Number=true, β::Number=false)
    mY,nY = size(y)
    mT = size(T,1)
    nA = size(A,2)
    mY == mT && nY == nA ||
        throw(DimensionMismatch("Dimensions of y ($(mY),$(nY)) must agree with first dimension of T ($(mT)), and second dimension of A ($(nA)), respectively"))

    b = y.block
    mY,nY = size(b)
    mA,nA = size(A.block)

    nrows = bandwidths(T)[1]
    mn = min(mA+nrows,size(T,1))
    mY == mn && nY == nA ||
        throw(DimensionMismatch("Dimensions of y subblock ($(mY),$(nY)) must agree with the dimensions of A subblock ($(mA)+nrows,$(nA))"))

    mul!(b, view(BandedMatrix(T), 1:mn, 1:mA), A.block)
    y
end

mul!(y::OneBlockMatrix, T::Union{Bidiagonal, SymTridiagonal, Tridiagonal}, A::OneBlockMatrix, args...) =
    _mul!(y, T, A, args...)

mul!(y::OneBlockMatrix, T::Union{BandedMatrix, SkewTridiagonal}, A::OneBlockMatrix, args...) =
    _mul!(y, T, A, args...)

function *(T::VariousBandedMatrices, A::OneBlockMatrix)
    m,n = size(A.block)
    U = Base.promote_op(*, eltype(T), eltype(A))
    M = size(T,1)
    nrows = bandwidths(T)[1]
    b = zeros(U, min(m+nrows,M), n)
    mul!(OneBlockMatrix(b, M, A.n), T, A)
end

function _mul!(y::OneBlockMatrix, A::OneBlockMatrix, T::VariousBandedMatrices,
               α::Number=true, β::Number=false)
    mY,nY = size(y)
    mA = size(A,1)
    nT = size(T,2)
    mY == mA && nY == nT ||
        throw(DimensionMismatch("Dimensions of y ($(mY),$(nY)) must agree with first dimension of A ($(mA)), and second dimension of T ($(nT)), respectively"))

    b = y.block
    mY,nY = size(b)
    mA,nA = size(A.block)

    ncols = bandwidths(T)[2]
    mn = min(nA+ncols,size(A,2))
    mY == mA && nY == mn ||
        throw(DimensionMismatch("Dimensions of y subblock ($(mY),$(nY)) must agree with the dimensions of A subblock ($(mA),$(nA)+1)"))

    mul!(b, A.block, view(BandedMatrix(T), 1:nA, 1:mn))
    y
end

mul!(y::OneBlockMatrix, A::OneBlockMatrix, T::Union{Bidiagonal, SymTridiagonal, Tridiagonal}, args...) =
    _mul!(y, A, T, args...)

mul!(y::OneBlockMatrix, A::OneBlockMatrix, T::Union{BandedMatrix, SkewTridiagonal}, args...) =
    _mul!(y, A, T, args...)

function *(A::OneBlockMatrix, T::VariousBandedMatrices)
    m,n = size(A.block)
    U = Base.promote_op(*, eltype(A), eltype(T))
    N = size(T,2)
    ncols = bandwidths(T)[2]
    b = zeros(U, m, min(n+ncols,N))
    mul!(OneBlockMatrix(b, A.m, N), A, T)
end

export OneBlockMatrix

# * Eigenfactorization

function LinearAlgebra.eigen(A::OneBlockMatrix{T}) where T
    M,N = size(A)
    M == N || throw(DimensionMismatch("matrix is not square: dimensions are ($(M), $(N))"))
    m,n = size(A.block)

    mn = max(m,n)
    Mmn = M-mn

    b = zeros(T, mn, mn)
    copyto!(view(b, 1:m, 1:n), A.block)

    ee = eigen(b)
    U = eltype(ee)
    λ = vcat(ee.values, zeros(U, Mmn))
    Φ = BlockSparseDiagonal(Diagonal(ones(U, M)), [ee.vectors - I], [1])

    Eigen(λ, Φ)
end
