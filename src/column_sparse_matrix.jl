struct ColumnSparseMatrix{T,Columns<:AbstractMatrix{T},Indices} <: AbstractMatrix{T}
    m::Int
    n::Int
    columns::Columns
    col_indices::Indices
    function ColumnSparseMatrix(m::Int, n::Int, columns::Columns, col_indices::Indices) where {Columns,Indices}
        if !issorted(col_indices)
            p = sortperm(col_indices)
            col_indices = col_indices[p]
            columns = copy(columns)
            Base.permutecols!!(columns, p)
        end
        new{eltype(columns), Columns, Indices}(m, n, columns, col_indices)
    end
end

function ColumnSparseMatrix(A::AbstractMatrix)
    m,n = size(A)
    col_indices = [j for j in 1:n
                   if !all(iszero, view(A, :, j))]
    ColumnSparseMatrix(m, n, A[:, col_indices], col_indices)
end

Base.size(A::ColumnSparseMatrix) = (A.m,A.n)

function Base.getindex(A::ColumnSparseMatrix{T}, i, j) where T
    k = findfirst(isequal(j), A.col_indices)
    if isnothing(k)
        zero(T)
    else
        A.columns[i,k]
    end
end

function Base.replace_in_print_matrix(A::ColumnSparseMatrix, i::Integer, j::Integer, s::AbstractString)
    j ∈ A.col_indices ? s : Base.replace_with_centered_mark(s)
end

function LinearAlgebra.dot(A::ColumnSparseMatrix, B::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B))
    v = zero(T)
    for (j,k) in enumerate(A.col_indices)
        v += dot(view(A.columns, :, j), view(B, :, k))
    end
    v
end

LinearAlgebra.dot(A::AbstractMatrix, B::ColumnSparseMatrix) =
    conj(dot(B, A))

"""
    intersect_with_indices(u, v)

Intersect `u` and `v` and return the elements present in both, as well
as their indices in `u` and `v`, respectively.

Credits: @MasonProtter
"""
function intersect_with_indices(u, v)
    out = NamedTuple{(:i, :j, :x), Tuple{Int, Int, Base.promote_eltypeof(u, v)}}[]
    foreach(eachindex(u)) do i
        @inbounds ui = u[i]
        j = findfirst(==(ui), v)
        if j !== nothing
            @inbounds vj = v[j]
            push!(out, (i=i, j=j, x=ui))
        end
    end
    out
end

function LinearAlgebra.dot(A::ColumnSparseMatrix, B::ColumnSparseMatrix)
    T = promote_type(eltype(A), eltype(B))
    v = zero(T)
    for (i,j,_) in intersect_with_indices(A.col_indices, B.col_indices)
        v += dot(view(A.columns, :, i), view(B.columns, :, j))
    end
    v
end

function Base.:(*)(A::AbstractMatrix, B::ColumnSparseMatrix)
    ma,na = size(A)
    mb,nb = size(B)
    na == mb || throw(DimensionMismatch("Sizes of A $((ma,na)) and B $((mb,nb)) not compatible ($(na) != $(mb))"))

    T = promote_type(eltype(A), eltype(B))
    columns = A*B.columns

    ColumnSparseMatrix(ma, nb, columns, B.col_indices)
end

function Base.:(*)(A::Adjoint{<:Any,<:ColumnSparseMatrix}, B::ColumnSparseMatrix)
    ma,na = size(A)
    mb,nb = size(B)
    na == mb || throw(DimensionMismatch("Sizes of A $((ma,na)) and B $((mb,nb)) not compatible ($(na) != $(mb))"))

    T = promote_type(eltype(A), eltype(B))

    I = Int[]
    J = Int[]
    V = T[]

    Ap = parent(A)
    for (k,i) in enumerate(Ap.col_indices)
        u = view(Ap.columns, :, k)
        for (l,j) in enumerate(B.col_indices)
            v = view(B.columns, :, l)
            push!(I, i)
            push!(J, j)
            push!(V, u'v)
        end
    end

    sparse(I, J, V, ma, nb)
end

Base.:(*)(α::Number, A::ColumnSparseMatrix) =
    ColumnSparseMatrix(A.m, A.n, α*A.columns, A.col_indices)

Base.:(*)(A::ColumnSparseMatrix, α::Number) = α*A

Base.:(/)(A::ColumnSparseMatrix, α::Number) = inv(α)*A

export ColumnSparseMatrix
