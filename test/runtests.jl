using WeirdMatrices
using LinearAlgebra
using BandedMatrices
using Test

@testset "WeirdMatrices.jl" begin
    @testset "OneBlockMatrices" begin
        N = 7
        A = OneBlockMatrix(Matrix(reshape(1:16, 4,4)), N, N)

        eye = Matrix{Int}(I, size(A))
        out = similar(eye)
        mul!(out, A, eye)
        @test out == Matrix(A)
        @test eye*A == Matrix(A)

        o = ones(Int, N)
        D = Diagonal(o)
        T = Tridiagonal(o[2:end], 2o, 3o[2:end])
        sT = SymTridiagonal(1o, 2o[2:end])
        S = SkewTridiagonal(1:N-1)
        B = BandedMatrix(2 => 1:N-2, 3=>1:N-3)
        Bu = Bidiagonal(1:N, 1:N-1, 'U')
        Bl = Bidiagonal(1:N, 1:N-1, 'L')

        function case(A,ref)
            @test A isa OneBlockMatrix
            @test A == ref
        end

        case(A, A*D)
        case(A, D*A)

        case(T*A, Matrix(T)*Matrix(A))
        @test Matrix(T)*A == Matrix(T)*Matrix(A)

        case(A*T, Matrix(A)*Matrix(T))
        @test A*Matrix(T) == Matrix(A)*Matrix(T)

        case(T*A*T, Matrix(T)*Matrix(A)*Matrix(T))

        case(A*(T*A), Matrix(A)*Matrix(T)*Matrix(A))
        case((A*T)*A, Matrix(A)*Matrix(T)*Matrix(A))

        case((A*T)*T, Matrix(A)*Matrix(T)*Matrix(T))

        case((A*T)*A, Matrix(A)*Matrix(T)*Matrix(A))

        case(T*(T*A), Matrix(T)*Matrix(T)*Matrix(A))

        case(sT*A, Matrix(sT)*Matrix(A))
        case(A*sT, Matrix(A)*Matrix(sT))
        case((sT*A)*sT, Matrix(sT)*Matrix(A)*Matrix(sT))
        case((A*sT)*sT, Matrix(A)*Matrix(sT)*Matrix(sT))
        case((sT*A)*T, Matrix(sT)*Matrix(A)*Matrix(T))
        case((A*sT)*T, Matrix(A)*Matrix(sT)*Matrix(T))

        case(S*A, Matrix(S)*Matrix(A))

        case(B*A, Matrix(B)*Matrix(A))
        case(B*B*A, Matrix(B)*Matrix(B)*Matrix(A))
        case(B*A*B, Matrix(B)*Matrix(A)*Matrix(B))

        case(Bl*A, Matrix(Bl)*Matrix(A))
        case(Bl*(Bl*A), Matrix(Bl)*Matrix(Bl)*Matrix(A))
        case(Bl*A*Bl, Matrix(Bl)*Matrix(A)*Matrix(Bl))

        case(Bu*A, Matrix(Bu)*Matrix(A))
        case(Bu*(Bu*A), Matrix(Bu)*Matrix(Bu)*Matrix(A))
        case(Bu*A*Bu, Matrix(Bu)*Matrix(A)*Matrix(Bu))
    end
end