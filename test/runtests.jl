using WeirdMatrices
using LinearAlgebra
using BandedMatrices
using SparseArrays
using Test

@testset "WeirdMatrices.jl" begin
    @testset "OneBlockMatrices" begin
        N = 7
        upper_block = Matrix(reshape(1:16, 4,4))
        A = OneBlockMatrix(upper_block, N, N)
        denseA = zeros(N, N)
        denseA[1:4,1:4] = upper_block

        for i = 1:N, j = 1:N
            @test A[i,j] == denseA[i,j]
        end
        @test A[1:5,3:7] == denseA[1:5,3:7]

        eye = Matrix{Int}(I, size(A))
        out = similar(eye)
        mul!(out, A, eye)
        @test out == Matrix(A)
        @test eye*A == Matrix(A)

        mA = -A
        @test mA isa OneBlockMatrix
        @test mA.block == -A.block

        @testset "Resultant OneBlockMatrix" begin
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

        @testset "Multiplication with arbitrary arrays" begin
            M = reshape(1.0:4N, N, 4)
            out = similar(M)
            mul!(out, A, M)
            @test out == Matrix(A)*M

            v = M[:,1]
            w = similar(v)
            mul!(w, A, v)
            @test w == Matrix(A)*v
        end
    end

    @testset "ColumnSparseMatrix" begin
        m,n = 10, 150
        u = rand(ComplexF64, m, 4)

        is = [2,4,10,6]
        A = ColumnSparseMatrix(m, n, u, is)
        @test issorted(A.col_indices)
        MA = Matrix(A)
        B = rand(ComplexF64, n, m)

        @test A*B ≈ MA*B

        C = rand(m, n)
        @test dot(A, C) ≈ dot(MA, C)
        @test dot(C, A) ≈ dot(C, MA)
        # display(@benchmark(dot(A,C)))
        # display(@benchmark(dot(MA,C)))

        @test dot(A,A) ≈ dot(MA,MA)
        # # Not as dramatic a difference as e.g. dot(A,C) vs dot(MA,C)
        # display(@benchmark(dot(A,A)))
        # display(@benchmark(dot(MA,MA)))

        CMA = ColumnSparseMatrix(MA)
        @test CMA isa ColumnSparseMatrix
        @test CMA == A

        BA = B*A
        @test BA isa ColumnSparseMatrix
        @test BA ≈ B*MA

        # display(@benchmark(B*A))
        # display(@benchmark(B*MA))AA = A'A

        AA = A'A
        @test AA isa SparseMatrixCSC
        @test AA ≈ MA'MA

        Ahalf = A/2
        @test Ahalf isa ColumnSparseMatrix
        @test Ahalf ≈ MA/2

        Atwice = A*2
        @test Atwice isa ColumnSparseMatrix
        @test Atwice ≈ MA*2
    end
end
