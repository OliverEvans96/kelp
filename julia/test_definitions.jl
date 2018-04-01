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

"""
No kelp. Uniform absorption & spacing over depth. Unit box.
"""
function calculate_light_no_kelp(a, b, num_scatters, nx, ny, nz, nθ, nϕ, I₀, θₛ, ϕₛ)

    num_si = 1
    si_ind = ones(nz)
    si_area = zeros(nz, num_si)

    abs_kelp = a
    scat_kelp = b
    abs_water = fill(a, nz)
    scat_water = fill(b, nz)

    rope_spacing = 1.0
    depth_spacing = fill(dz, nz)
    current_angles = zeros(nz)
    current_speeds = zeros(nz)

    num_vsf = 55
    vsf_file="data/vsf/nuc_vsf.txt"

    nω = nθ * (nϕ-2) + 2
    radiance = zeros(nx,ny,nz,nω)
    irradiance = zeros(nx,ny,nz)

    # Frond shape
    fs = 0.5
    fr = 2.0
    ft = 0.1

    ccall((:__light_interface_MOD_full_light_calculations, "light_interface"), Void,
          (
              # OPTICAL PROPERTIES
              # abs_kelp
              Ref{Float64},
              # scat_kelp
              Ref{Float64},
              # abs_water
              Ref{Float64},
              # scat_water
              Ref{Float64}
              # num_vsf,
              Ref{Int64},
              # vsf_file
              Ref{Cstring},
              # SUNLIGHT
              # solar_zenith
              Ref{Float64},
              # solar_azimuthal
              Ref{Float64},
              # surface_irrad
              Ref{Float64},
              # KELP
              # num_si
              Ref{Int64},
              # si_area
              Ref{Float64},
              # si_ind
              Ref{Int64},
              # frond_thickness
              Ref{Float64},
              # frond_aspect_ratio
              Ref{Float64},
              # frond_shape_ratio
              Ref{Float64},
              # WATER CURRENT
              # current_speeds
              Ref{Float64},
              # current_angles
              Ref{Float64},
              # SPACING
              # rope_spacing
              Ref{Float64},
              # depth_spacing
              Ref{Float64},
              # SOLVER PARAMETERS
              # nx
              Ref{Int64},
              # ny
              Ref{Int64},
              # nz
              Ref{Int64},
              # ntheta
              Ref{Int64},
              # nphi
              Ref{Int64},
              # num_scatters
              Ref{Int64},
              # ! FINAL RESULTS
              # available_light
              Ref{Float64},
              # avg_irrad
              Ref{Float64},
              # radiance
              Ref{Float64},
              # irradiance
              Ref{Float64}
          ),
          # OPTICAL PROPERTIES
          abs_kelp,
          scat_kelp,
          abs_water,
          scat_water,
          num_vsf,
          vsf_file,
          # SUNLIGHT
          ϕₛ,
          θₛ,
          I₀,
          # KELP
          num_si,
          si_area,
          si_ind,
          ft,
          fr,
          fs,
          # WATER CURRENT
          current_speeds,
          current_angles,
          # SPACING
          rope_spacing,
          depth_spacing,
          # SOLVER PARAMETERS
          nx,
          ny,
          nz,
          nθ,
          nϕ,
          num_scatters,
          # ! FINAL RESULTS
          available_light,
          avg_irrad,
          radiance,
          irradiance
        )

    return available_light, avg_irrad, radiance, irradiance
