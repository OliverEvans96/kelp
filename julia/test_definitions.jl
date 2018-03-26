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

function test_traverse(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nθ, nϕ, i, j, k, l, m)
    max_cells = calculate_max_cells(xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nθ, nϕ)
    nω = nθ*(nϕ-2)+2

    println("max_cells = $max_cells ($(typeof(max_cells)))")
    s = zeros(max_cells)
    ds = zeros(max_cells)
    ã = zeros(max_cells)
    gₙ = zeros(max_cells)
    rad_scatter = zeros(nx, ny, nz, nω)
    num_cells = [0]

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
              Ref{Int64}
          ),
          xmin, xmax, ymin, ymax, zmin, zmax,
          nx, ny, nz, nθ, nϕ,
          i, j, k, l, m,
          s, ds, ã, gₙ, rad_scatter, num_cells)

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

