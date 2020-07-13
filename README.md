# WeirdMatrices.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jagot.github.io/WeirdMatrices.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jagot.github.io/WeirdMatrices.jl/dev)
[![Build Status](https://github.com/jagot/WeirdMatrices.jl/workflows/CI/badge.svg)](https://github.com/jagot/WeirdMatrices.jl/actions)
[![Coverage](https://codecov.io/gh/jagot/WeirdMatrices.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jagot/WeirdMatrices.jl)

A home for weird matrices that do not fit in anywhere else.

```julia
julia> SkewTridiagonal(1:10)
11×11 SkewTridiagonal{Int64,UnitRange{Int64}}:
 ⋅  -1   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅    ⋅
 1   ⋅  -2   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅    ⋅
 ⋅   2   ⋅  -3   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅    ⋅
 ⋅   ⋅   3   ⋅  -4   ⋅   ⋅   ⋅   ⋅   ⋅    ⋅
 ⋅   ⋅   ⋅   4   ⋅  -5   ⋅   ⋅   ⋅   ⋅    ⋅
 ⋅   ⋅   ⋅   ⋅   5   ⋅  -6   ⋅   ⋅   ⋅    ⋅
 ⋅   ⋅   ⋅   ⋅   ⋅   6   ⋅  -7   ⋅   ⋅    ⋅
 ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   7   ⋅  -8   ⋅    ⋅
 ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   8   ⋅  -9    ⋅
 ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   9   ⋅  -10
 ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  10    ⋅

julia> OneBlockMatrix(Matrix(reshape(1:15, 5, 3)), 10, 20)
10×20 OneBlockMatrix{Int64,Array{Int64,2}}:
 1   6  11  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 2   7  12  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 3   8  13  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 4   9  14  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 5  10  15  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
 ⋅   ⋅   ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅  ⋅
```

`OneBlockMatrix` are also closed under some operations:

```julia
julia> N = 7
7

julia> A = OneBlockMatrix(Matrix(reshape(1:16, 4,4)), N, N)
7×7 OneBlockMatrix{Int64,Array{Int64,2}}:
 1  5   9  13  ⋅  ⋅  ⋅
 2  6  10  14  ⋅  ⋅  ⋅
 3  7  11  15  ⋅  ⋅  ⋅
 4  8  12  16  ⋅  ⋅  ⋅
 ⋅  ⋅   ⋅   ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅   ⋅   ⋅  ⋅  ⋅  ⋅
 ⋅  ⋅   ⋅   ⋅  ⋅  ⋅  ⋅

julia> o = ones(Int, N)
7-element Array{Int64,1}:
 1
 1
 1
 1
 1
 1
 1

julia> T = Tridiagonal(o[2:end], 2o, 3o[2:end])
7×7 Tridiagonal{Int64,Array{Int64,1}}:
 2  3  ⋅  ⋅  ⋅  ⋅  ⋅
 1  2  3  ⋅  ⋅  ⋅  ⋅
 ⋅  1  2  3  ⋅  ⋅  ⋅
 ⋅  ⋅  1  2  3  ⋅  ⋅
 ⋅  ⋅  ⋅  1  2  3  ⋅
 ⋅  ⋅  ⋅  ⋅  1  2  3
 ⋅  ⋅  ⋅  ⋅  ⋅  1  2

julia> T*A
7×7 OneBlockMatrix{Int64,Array{Int64,2}}:
  8  28  48  68  ⋅  ⋅  ⋅
 14  38  62  86  ⋅  ⋅  ⋅
 20  44  68  92  ⋅  ⋅  ⋅
 11  23  35  47  ⋅  ⋅  ⋅
  4   8  12  16  ⋅  ⋅  ⋅
  ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅
  ⋅   ⋅   ⋅   ⋅  ⋅  ⋅  ⋅

julia> A^2
7×7 OneBlockMatrix{Int64,Array{Int64,2}}:
  90  202  314  426  ⋅  ⋅  ⋅
 100  228  356  484  ⋅  ⋅  ⋅
 110  254  398  542  ⋅  ⋅  ⋅
 120  280  440  600  ⋅  ⋅  ⋅
   ⋅    ⋅    ⋅    ⋅  ⋅  ⋅  ⋅
   ⋅    ⋅    ⋅    ⋅  ⋅  ⋅  ⋅
   ⋅    ⋅    ⋅    ⋅  ⋅  ⋅  ⋅
```
