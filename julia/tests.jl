using Base.Test
using Cubature

# Fortran Wrappers

function test_2d_angular_integration(f, nθ, nϕ)
    fptr = cfunction(f, Float64, (Ref{Float64}, Ref{Float64},))
    ccall((:__test_grid_MOD_test_2d_angular_integration,
           "test_grid"),
          Float64, (Ptr{Void}, Ref{Int64}, Ref{Int64}), fptr, nθ, nϕ)
end

function test_angle_p_conversions(nθ, nϕ)
    ccall((:__test_grid_MOD_test_angle_p_conversions,
           "test_grid"),
          Bool, (Ref{Int64}, Ref{Int64}), nθ, nϕ)
end

# Julia Functions

function gengrid(nθ, nϕ)
    """generate grid. nϕ includes poles."""
    dθ = 2π/(nθ)
    dϕ = π/(nϕ-1)

    l = 1:nθ
    m = 1:nϕ
    θ = collect((l-1)*dθ)
    ϕ = collect((m-1)*dϕ)

    # edges between cells
    θe = θ+dθ/2
    ϕe = ϕ[1:nϕ-1]+dϕ/2
    #ϕe[1] = dϕ/2
    #ϕe[nϕ-1] = π - dϕ/2

    return dθ, θ, θe, dϕ, ϕ, ϕe
end

function angularquad(f, nθ, nϕ)
    dθ, θ, θe, dϕ, ϕ, ϕe = gengrid(nθ, nϕ)

    f_arr = [f(θ1, ϕ1) for θ1=θ, ϕ1=ϕ]

    # initialize with value only considering poles
    integ = 2π*(1-cos(dϕ/2))*(f_arr[1,1]+f_arr[1,nϕ])
    # exclude poles
    for l=1:nθ
        for m=2:nϕ-1
            integ += f_arr[l,m] * dθ * (cos(ϕe[m-1])-cos(ϕe[m]))
        end
    end

    return integ
end

function angularcubature(f, tol=1e-5)
    cubature, cube_err = hcubature(ω -> begin (θ, ϕ) = ω; sin(ϕ)*f(θ, ϕ) end, (0.,0.), (2π,π), reltol=tol, abstol=tol)
    return cubature
end

@testset "Grid Tests" begin
    @testset "Angular Integral" begin
        @testset "Julia Grid" begin
            f(θ,ϕ) = 1+.1*(sin(10*ϕ)+cos(10*θ))
            for nθ in [5,10,30]
                for nϕ in [5,10,30]
                    @test angularquad(f, nθ, nϕ) ≈ test_2d_angular_integration(f, nθ, nϕ) atol=1e-5
                end
            end
        end

        @testset "Cubature" begin
            f(θ,ϕ) = 1+.1*(sin(10*ϕ)+cos(10*θ))
            g(θ,ϕ) = ϕ^2 - log(ϕ+4.) - cos(5θ)
            h(θ,ϕ) = sin(2*θ)*sin(3*ϕ) + sin(θ)^2-sin(ϕ^2)^3

            for func in (f, g, h)
                @test angularcubature(func, 1e-5) ≈ angularquad(func, 100, 100) atol=1e-3
            end
        end
    end

    @testset "P-Conversions" begin
        for nθ=[5,10,30] nϕ=[5,10,30]
            @test test_angle_p_conversions(nθ, nϕ)
        end
    end
end

