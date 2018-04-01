using Base.Test
include("test_definitions.jl")

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

@testset "Asymptotics" begin
    @testset "Traverse" begin
        xmin = -6
        ymin = -6
        zmin = -6
        xmax = 6
        ymax = 6
        zmax = 6
        lims = [xmin ymin zmin; xmax ymax zmax]
        nx = 10
        ny = 10
        nz = 10
        nθ = 10
        nϕ = 10
        nums = [nx ny nz nθ nϕ]
        # i, j, k, l, m
        indices = (3, 5, 4, 5, 3)
        s, ds, ã, gₙ, rad_scatter = test_traverse(lims..., nums..., indices...)

        println("ã = $ã")
        println("gₙ = $gₙ")
        println("s = $s")
        println("ds = $ds")
        # No NaN values
        @test all(.!isnan.(ã))
        @test all(.!isnan.(gₙ))
        @test all(.!isnan.(s))
        @test all(.!isnan.(ds))
        @test all(.!isnan.(rad_scatter))

        # s increasing
        @test all(diff(s).≥0)

    end
    @testset "Path Integration" for test_num=1:5
        # Grid vs analytical integral
        xmin = -6
        ymin = -6
        zmin = -6
        xmax = 6
        ymax = 6
        zmax = 6
        lims = [xmin ymin zmin; xmax ymax zmax]
        nx = 30
        ny = 30
        nz = 30
        nθ = 12
        nϕ = 12
        nums = [nx ny nz nθ nϕ]
        # i, j, k, l, m
        # Choose random indices
        indices = [rand(1:num) for num in nums]

        println("inds = $indices")

        pkelp_fun(x,y,z) = 1+sin(x) + 4*(2+sin(y)^2) + (1+cos(x*y+z))
        s, ds, ã, gₙ, rad_scatter = test_traverse(lims..., nums..., indices..., pkelp_fun)
        # println("ã = $ã")
        # println("gₙ = $gₙ")
        # println("s = $s")
        # println("ds = $ds")

        grid_integ = ds'*ã
        quad_integ = test_ray_integral(lims...,nums...,indices...,pkelp_fun)
        println("grid_integ = $grid_integ")
        println("quad_integ = $quad_integ")
        println()
        @test grid_integ > 0
        # Hopefully grid integration is accurate to 20%
        @test grid_integ ≈ quad_integ rtol=0.20
    end

    @testset "Full Scattering" begin
    end

end


