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

function make_grid(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nθ, nϕ)
    dx = zeros(nx)
    dy = zeros(ny)
    dz = zeros(nz)
    # Julia doesn't pass by value, so have to wrap as array
    # Fortran still thinks it's a single number.
    dθ = [0.]
    dϕ = [0.]

    x = zeros(nx)
    y = zeros(ny)
    z = zeros(nz)
    θ = zeros(nθ)
    ϕ = zeros(nϕ)

    xe = zeros(nx)
    ye = zeros(ny)
    ze = zeros(nz)
    θe = zeros(nθ)
    ϕe = zeros(nϕ-1)

    ccall((:__test_grid_MOD_make_grid, "test_grid"), Void,
            (
                # Bounds
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},

                # Num
                Ref{Int64},
                Ref{Int64},
                Ref{Int64},
                Ref{Int64},
                Ref{Int64},

                # Spacing
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},

                # Vals
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},

                # Edges
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
                Ref{Float64},
            ),
          xmin, xmax, ymin, ymax, zmin, zmax,
          nx, ny, nz, nθ, nϕ,
          dx, dy, dz, dθ, dϕ,
          x, y, z, θ, ϕ,
          xe, ye, ze, θe, ϕe)
    cells = x, y, z, θ, ϕ
    edges = xe, ye, ze, θe, ϕe
    spacing = dx, dy, dz, dθ[1], dϕ[1]
    return cells, edges, spacing
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
    @testset "Cell Locations" for n=10:10:50
        xmin = -6
        ymin = -6
        zmin = -6
        xmax = 6
        ymax = 6
        zmax = 6
        lims = [xmin ymin zmin; xmax ymax zmax]
        nums = repeat([n], outer=5)
        cells, edges, spacing = make_grid(lims...,nums...)
        x, y, z, θ, ϕ = cells
        xe, ye, ze, θe, ϕe = edges
        dx, dy, dz, dθ, dϕ = spacing
        @testset "Even Spacing" for (cell, edge) in zip(cells, edges)
            dc = diff(cell)
            @test maximum(dc) ≈ minimum(dc) atol=1e-6
            de = diff(edge)
            @test maximum(de) ≈ minimum(de) atol=1e-6
        end
        @testset "Cell Endpoints" begin
            @test x[1] ≈ xmin+dx[1]/2 atol=1e-6
            @test y[1] ≈ ymin+dy[1]/2 atol=1e-6
            @test z[1] ≈ zmin+dz[1]/2 atol=1e-6
            @test θ[1] ≈ 0 atol=1e-6
            @test ϕ[1] ≈ 0 atol=1e-6

            @test x[n] ≈ xmax-dx[n]/2 atol=1e-6
            @test y[n] ≈ ymax-dy[n]/2 atol=1e-6
            @test z[n] ≈ zmax-dz[n]/2 atol=1e-6
            @test θ[n] ≈ 2π-dθ atol=1e-6
            @test ϕ[n] ≈ π atol=1e-6
        end
        @testset "Edge Endpoints" begin
            @test xe[1] ≈ xmin atol=1e-6
            @test ye[1] ≈ ymin atol=1e-6
            @test ze[1] ≈ zmin atol=1e-6
            @test θe[1] ≈ dθ/2 atol=1e-6
            @test ϕe[1] ≈ dϕ/2 atol=1e-6

            @test xe[n] ≈ xmax-dx[n] atol=1e-6
            @test ye[n] ≈ ymax-dy[n] atol=1e-6
            @test ze[n] ≈ zmax-dz[n] atol=1e-6
            @test θe[n] ≈ 2π-dθ/2 atol=1e-6
            @test ϕe[n-1] ≈ π-dϕ/2 atol=1e-6
        end
    end

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
