# Step 1: Generate kelp from ind and area
module KelpModel

export ComputationalDomainPoints
export ComputationalDomainSpacing
export GetGrid

"""
3D Bounding box for domain
[xmin, xmax;
 ymin, ymax;
 zmin, zmaz]
"""
struct BoundingBox
    bounds::Matrix{<:Number}

    xmin::Number
    xmax::Number
    ymin::Number
    ymax::Number
    zmin::Number
    zmax::Number

    "Create box from array"
    function BoundingBox(bounds::Matrix{<:Number})
        if size(bounds) != (3,2)
            error("must be 3x2 Array")
        end

        xmin, xmax = bounds[1,:]
        ymin, ymax = bounds[2,:]
        zmin, zmax = bounds[3,:]

        verifyBounds(new(bounds, xmin, xmax, ymin, ymax, zmin, zmax))
    end

    "Create box from unnamed keyword parameters"
    function BoundingBox(xmin::F, xmax::F, ymin::F, ymax::F,
            zmin::F, zmax::F) where F<:Number
        bounds = [xmin xmax; ymin ymax; zmin zmax]

        verifyBounds(new(bounds, xmin, xmax, ymin, ymax, zmin, zmax))
    end

    "Create box from named keyword parameters"
    function BoundingBox(;xmin::F=0, xmax::F=1, ymin::F=0, ymax::F=1,
            zmin::F=0, zmax::F=1) where F<:Number
        bounds = [xmin xmax; ymin ymax; zmin zmax]

        verifyBounds(new(bounds, xmin, xmax, ymin, ymax, zmin, zmax))
    end

    "Disallow specifying both bounds and keyword parameters"
    function BoundingBox(kws::Vararg{Any,7})
        error("specify either array or keyword args, not both")
    end
end

"Verify that min < max for each dimension in BoundingBox"
function verifyBounds(box::BoundingBox)
    all(box.bounds[:,1] .< box.bounds[:,2]) ? box :
        error("Must have min < max in each dimension")
end

"""
Contains number of and spacing between points in a uniform grid.
GridResolution([nx ny nz] [dx dy dz])
"""
struct GridResolution
    npoints::Vector{<:Integer}
    nx::Integer
    ny::Integer
    nz::Integer

    spacing::Vector{<:Number}
    dx::Number
    dy::Number
    dz::Number

    function GridResolution(npoints::Vector{<:Integer}, spacing::Vector{<:Number})
        nx, ny, nz = npoints
        dx, dy, dz = spacing

        verifyGridResolution(new(npoints, nx, ny, nz, spacing, dx, dy, dz))
    end
end

"Verify that grid resolution is positive"
function verifyGridResolution(resolution::GridResolution)
    all([resolution.npoints resolution.spacing] .> 0) ? resolution :
        error("Grid resolution must be positive")
end


"""
3D uniform computationl domain, consisting of bounds,
number of points in each dimension, and spacing between
points in each dimension.

Either spacing or number of points may be specified.
"""
struct ComputationalDomain
    # Box bounds
    box::BoundingBox

    # Grid resolution
    mesh::GridResolution

    # Grid ranges
    x::Range{<:Number}
    y::Range{<:Number}
    z::Range{<:Number}

    # Whether endpoints are included
    endpoints::Bool
end

"Create domain from spacing, exclude upper endpoints by default"
function ComputationalDomainSpacing(bounds::Matrix{<:Number}, spacing::Vector{<:Number}, endpoints::Bool=false)
    if size(bounds) == (3,2) && length(spacing) == 3
        box::BoundingBox = BoundingBox(bounds)
        x::Range = box.xmin : spacing[1] : box.xmax
        y::Range = box.ymin : spacing[2] : box.ymax
        z::Range = box.zmin : spacing[3] : box.zmax

        # If we're removing endpoints, only remove the last point if it's
        # actually the maximum value (i.e., the endpoint)
        if !endpoints
            x = x[1:end - (x[end] == box.xmax)]
            y = y[1:end - (y[end] == box.ymax)]
            z = z[1:end - (z[end] == box.zmax)]
        end

        nx::Integer = length(x)
        ny::Integer = length(y)
        nz::Integer = length(z)
        npoints = [nx, ny, nz]

        mesh = GridResolution(npoints, spacing)

        ComputationalDomain(box, mesh, x, y, z, endpoints)
    else
        error("bounds must be a 3x2 matrix, spacing must be a length 3 vector")
    end
end

"Create domain from number of points, exclude upper endpoints by default"
function ComputationalDomainPoints(bounds::Matrix{<:Number}, npoints::Vector{<:Integer}, endpoints::Bool=false)
    if size(bounds) == (3,2) && length(npoints) == 3
        spacing::Vector{<:Number} = (diff(bounds, 2) ./ (npoints - endpoints))[:]

        ComputationalDomainSpacing(bounds, spacing, endpoints)
    else
        error("bounds must be a 3x2 matrix, npoints must be a length 3 vector")
    end
end

struct OpticalProperties
    vsf::Function
    absorption::Number
end

function createkelp(nz::Integer, ind::Array, area::Array)
end

function propagatelight()
end

function absorblight()
end

end
