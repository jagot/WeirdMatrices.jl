struct BlockSparseDiagonal{T,Diag<:Diagonal{T},Blocks} <: AbstractMatrix{T}
    diag::Diag
    blocks::Blocks
    block_starts::Vector{Int}
end

Base.size(A::BlockSparseDiagonal, args...) = size(A.diag, args...)

function find_block(A::BlockSparseDiagonal, i, j)
    k = findlast(≤(i), A.block_starts)
    isnothing(k) && return nothing
    bs = A.block_starts[k]
    interval = bs:bs+size(A.blocks[k],1)-1
    i ∈ interval && j ∈ interval || return nothing
    k
end

function Base.getindex(A::BlockSparseDiagonal{T}, i, j) where T
    v = zero(T)
    b = find_block(A, i, j)
    if !isnothing(b)
        bs = A.block_starts[b]
        v += A.blocks[b][i-bs+1,j-bs+1]
    end
    if i == j
        v += A.diag.diag[i]
    end
    v
end

Base.replace_in_print_matrix(A::BlockSparseDiagonal,i::Integer,j::Integer,s::AbstractString) =
    (i == j || !isnothing(find_block(A, i, j))
     ? s
     : Base.replace_with_centered_mark(s))

Base.:(*)(α::Number, A::BlockSparseDiagonal) =
    BlockSparseDiagonal(α*A.diag, α .* A.blocks, A.block_starts)

Base.:(*)(A::BlockSparseDiagonal, α::Number) =
    BlockSparseDiagonal(α*A.diag, α .* A.blocks, A.block_starts)

Base.:(/)(A::BlockSparseDiagonal, α::Number) = A * inv(α)

function LinearAlgebra.mul!(y::AbstractMatrix, A::BlockSparseDiagonal, x::AbstractMatrix,
                            α::Number=true, β::Number=false)
    mul!(y, A.diag, x, α, β)
    for (b,bs) in zip(A.blocks, A.block_starts)
        interval = bs:bs+size(b,2)-1
        mul!(view(y, interval, :), b, view(x, interval, :), α, true)
    end
    y
end

export BlockSparseDiagonal
