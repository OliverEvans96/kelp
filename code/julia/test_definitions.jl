module KelpTest

using Cubature
using Interpolations

# Single integral alias
∫(args...) = hquadrature(args...)[1]
export ∫

function test_fun(args...)
    println(args)
end
export test_fun

# Fortran Wrappers

function test_make_vsf(nθ, nϕ)
    θ = zeros(nθ)
    ϕ = zeros(nϕ)
    θe = zeros(nθ)
    ϕe = zeros(nϕ-1)
    dθ = [0.0]
    dϕ = [0.0]
    nω = nθ * (nϕ-2) + 2
    θₚ = zeros(nω)
    ϕₚ = zeros(nω)
    aₚ = zeros(nω)
    β = zeros(nω, nω)
    βᵢ = zeros(nω, nω)

    num_vsf = 55
    vsf_angles = zeros(num_vsf)
    vsf_vals = zeros(num_vsf)
    ccall((:__test_grid_MOD_make_vsf, "test_grid"),
          Void, (
              Ref{Int64},
              Ref{Int64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64}
          ),
          nθ, nϕ, θ, ϕ, θe, ϕe, dθ, dϕ,
          θₚ, ϕₚ, aₚ, β, βᵢ, vsf_angles, vsf_vals
    )
    return θ, ϕ, θe, ϕe, dθ[1], dϕ[1], θₚ, ϕₚ, aₚ, β, βᵢ, vsf_angles, vsf_vals
end
export test_make_vsf

function test_2d_angular_integration(f, nθ, nϕ)
    fptr = cfunction(f, Float64, (Ref{Float64}, Ref{Float64},))
    ccall((:__test_grid_MOD_test_2d_angular_integration,
           "test_grid"),
          Float64, (Ptr{Void}, Ref{Int64}, Ref{Int64}), fptr, nθ, nϕ)
end
export test_2d_angular_integration

function test_angle_p_conversions(nθ, nϕ)
    ccall((:__test_grid_MOD_test_angle_p_conversions,
           "test_grid"),
          Bool, (Ref{Int64}, Ref{Int64}), nθ, nϕ)
end
export test_angle_p_conversions

function calculate_max_cells(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nθ, nϕ)
    max_cells = ccall((:__test_asymptotics_MOD_test_max_cells, "test_asymptotics"),
        Int64,
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
        ),
        xmin, xmax, ymin, ymax, zmin, zmax,
        nx, ny, nz, nθ, nϕ)
    return max_cells
end
export calculate_max_cells

function test_ray_integral(xmin, xmax, ymin, ymax, zmin, zmax,
                           nx, ny, nz, nθ, nϕ,
                           i, j, k, l, m, pkelp_fun)

    cells, edges, spacing = make_grid(
        xmin, xmax, ymin, ymax, zmin, zmax,
        nx, ny, nz, nθ, nϕ
    )

    x, y, z, θ, ϕ = cells

    vec_x = [x[i],y[j],z[k]]
    vec_omega = [
        sin(ϕ[m])*cos(θ[l]),
        sin(ϕ[m])*sin(θ[l]),
        cos(ϕ[m])
    ]

    z0 = vec_omega[3] < 0 ? zmax : zmin
    s̃ = (vec_x[3] - z0) / vec_omega[3]
    vec_x0 = vec_x - s̃ * vec_omega

    println("vec_x = $vec_x")
    println("vec_omega = $vec_omega")
    println("z0 = $z0")
    println("s̃ = $s̃")
    println("vec_x0 = $vec_x0")

    # Path function
    # This is "l"
    path_fun(s) = vec_x0 + s/s̃ * (vec_x - vec_x0)

    shift_mod(x, xmin, xmax) = xmin + mod(x-xmin, xmax-xmin)
    # Make path function periodic
    function path_fun_per(s)
        path = path_fun(s)
        return [
            shift_mod(path[1], xmin, xmax),
            shift_mod(path[2], ymin, ymax),
            shift_mod(path[3], zmin, zmax),
        ]
    end

    # Assume that a_kelp = 1.0, a_water = 1.0,
    # so that a(x,y,z) = p_kelp(x,y,z)
    ã(s) = pkelp_fun(path_fun_per(s)...)

    # Integrate numerically
    tol = 1.0e-4
    integral, err = hquadrature(ã, 0.0, s̃, reltol=tol, abstol=tol)

    return integral
end
export test_ray_integral

function test_traverse(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nθ, nϕ, i, j, k, l, m, pkelp_jfun=((args...) -> 1.0))
    max_cells = calculate_max_cells(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nθ, nϕ)
    nω = nθ*(nϕ-2)+2

    println("max_cells = $max_cells ($(typeof(max_cells)))")
    s = zeros(max_cells)
    ds = zeros(max_cells)
    ã = zeros(max_cells)
    gₙ = zeros(max_cells)
    rad_scatter = zeros(nx, ny, nz, nω)
    num_cells = [0]

    # Convert to fortran-callable function
    # To determine pkelp. With abs_kelp = 1.0, abs_water = 0.0,
    # this function effectively determines abs_grid.
    pkelp_cfun = cfunction(
        pkelp_jfun, Float64,
        (Ref{Float64}, Ref{Float64}, Ref{Float64})
    )

    ccall((:__test_asymptotics_MOD_test_traverse,"test_asymptotics"),
          Void, (
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

              # Indices
              Ref{Int64},
              Ref{Int64},
              Ref{Int64},
              Ref{Int64},
              Ref{Int64},

              # Results
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Float64},
              Ref{Int64},

              # Kelp Function
              Ptr{Void}
          ),
          xmin, xmax, ymin, ymax, zmin, zmax,
          nx, ny, nz, nθ, nϕ,
          i, j, k, l, m,
          s, ds, ã, gₙ, rad_scatter,
          num_cells, pkelp_cfun)

    num_cells = num_cells[1]
    s = s[1:num_cells]
    ds = ds[1:num_cells]
    ã = ã[1:num_cells]
    gₙ = gₙ[1:num_cells]

    return s, ds, ã, gₙ, rad_scatter
end
export test_traverse

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
export make_grid

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
export gen_grid

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
export angularquad

function angularcubature(f, tol=1e-5)
    cubature, cube_err = hcubature(ω -> begin (θ, ϕ) = ω; sin(ϕ)*f(θ, ϕ) end, (0.,0.), (2π,π), reltol=tol, abstol=tol)
    return cubature
end
export angularcubature

function rte1d_exact(I₀, a, b, zmin, zmax, z)
    z₀, z₁ = zmin, zmax
    s = a/b + 1
    q = sqrt(a^2+2*a*b)/b
    c₂ = -(s+q)*I₀*exp(b*q*(2*z₁-z₀))/(s-q-(s+q)*exp(2b*q*(z₁-z₀)))
    c₁ = exp(-b*q*z₀) * (I₀ - c₂*exp(-b*q*z₀))
    L⁺ = c₁*exp.(b*q*z) + c₂*exp.(-b*q*z)
    L⁻ = c₁*(s+q)*exp.(b*q*z) + c₂*(s-q)*exp.(-b*q*z)
    return L⁺, L⁻
end
export rte1d_exact

function asymptotics1d_exact(I₀::Float64, a::Float64, b::Float64, βπ::Float64,
        zmin::Float64, zmax::Float64, z::Array{Float64,1}, num_scatters::Int)
    # Get interpolation index from z
    indx(z′) = (z′-z[1])/(z[end]-z[1]) * (length(z) - 1) + 1

    # % of scattered light remaining in angular grid cell
    β₀ = 1 - βπ

    # Radiance without scattering
    L₀⁺ = I₀ * exp.(-a*(z.-zmin))
    L₀⁻ = zeros(z)

    # Radiance to use for source calculation
    Lₙ₋₁⁺ = L₀⁺[:]
    Lₙ₋₁⁻ = L₀⁻[:]

    # Initialize final radiance
    L⁺ = L₀⁺[:]
    L⁻ = L₀⁻[:]

    # Debug variables
    Lₙ = zeros(length(z), 2, num_scatters+1)
    gₙ = zeros(length(z), 2, num_scatters+1)

    Lₙ[:,1,1] = L₀⁺
    Lₙ[:,2,1] = L₀⁻

    for n = 1:num_scatters
        # Calculate source and interpolate
        gₙ⁺ = z′ -> (βπ * interpolate(Lₙ₋₁⁻ - Lₙ₋₁⁺, BSpline(Linear()), OnCell())[indx(z′)])
        gₙ⁻ = z′ -> (βπ * interpolate(Lₙ₋₁⁺ - Lₙ₋₁⁻, BSpline(Linear()), OnCell())[indx(z′)])

        # for (k, zz) in enumerate(z)
        #     println("k = $k. βπ(L⁺ - L⁻) = $(βπ*(Lₙ₋₁⁺[k] - Lₙ₋₁⁻[k])) =? $(gₙ⁻(zz))")
        # end

        # Define nested anonymous functions and apply elementwize to z
        Lₙ⁺ = (z -> ∫(zp -> exp(-a*(z-zp)) * gₙ⁺(zp), zmin, z)).(z)
        Lₙ⁻ = (z -> ∫(zp -> exp(-a*(zp-z)) * gₙ⁻(zp), z, zmax)).(z)

        # Debug
        gₙ[:,1,n] = gₙ⁺.(z)
        gₙ[:,2,n] = gₙ⁻.(z)
        Lₙ[:,1,n+1] = Lₙ⁺
        Lₙ[:,2,n+1] = Lₙ⁻

        # Accumulate total radiance
        L⁺ += b^n * Lₙ⁺
        L⁻ += b^n * Lₙ⁻

        # Update
        Lₙ₋₁⁺[:] = Lₙ⁺[:]
        Lₙ₋₁⁻[:] = Lₙ⁻[:]
    end

    return L⁺, L⁻, Lₙ, gₙ
end
export asymptotics1d_exact

function asymptotics1d_grid(I₀, a, b, β̃, zmin, zmax, nz, num_scatters)
    z = zeros(nz)
    L⁺ = zeros(nz)
    L⁻ = zeros(nz)
    vsf_cfunc = cfunction(β̃, Float64, (Ref{Float64},))

    # This syntax allows updating library between calls
    # See https://discourse.julialang.org/t/unload-a-shared-library/5344/4
    funcsym = :__test_asymptotics_MOD_test_asymptotics_1d
    lib = Libdl.dlopen("../../include/test_asymptotics.so")
    sym = Libdl.dlsym(lib, funcsym)
    ccall(sym,
          Void,
          (Ref{Float64},
           Ref{Float64},
           Ref{Float64},
           Ptr{Void},
           Ref{Float64},
           Ref{Float64},
           Ref{Int64},
           Ref{Float64},
           Ref{Float64},
           Ref{Float64},
           Ref{Int64}),
          I₀, a, b, vsf_cfunc,
          zmin, zmax, nz, z,
          L⁺, L⁻, num_scatters
          )
    Libdl.dlclose(lib)

    return z, L⁺, L⁻

end
export asymptotics1d_grid

function p̂(l,m)
    if m == 1
        p = 1
    elseif m == nϕ
        p = nω
    else
        p = (m-2)*nθ + l + 1
    end
end
export p̂

function m̂(p)
    if p == 1
        m = 1
    elseif p == nω
        m = nϕ
    else
        m = ceil(Int, (p-1)/nθ) + 1
    end
end
export m̂

function l̂(p)
    if p == 1
        l = 1
    elseif p == nω
        l = 1
    else
        m = mod1(p-1, nθ)
    end
end
export l̂

end
