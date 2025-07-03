using GeometryBasics

"""
    @enum CELLTYPE

Enumeration representing different types of cells.

## Fields
- `inner`: Represents an inner cell.
- `inlet`: Represents an inlet cell.
- `outlet`: Represents an outlet cell.
- `wall`: Represents a wall cell.
"""
@enum CELLTYPE begin
    inner = 1
    inlet = -1
    outlet = -2
    wall = -3
end

"""
    @enum Verbosity

Enumeration representing different levels of verbosity.

## Fields
- `silent`: Represents silent verbosity.
- `verbose`: Represents verbose verbosity.
"""
@enum Verbosity begin
    silent = 0
    verbose = 2
end


"""
    @enum ModelType

Enumeration representing different types of models.

## Fields
- `model_1`: Represents model 1.
- `model_2`: Represents model 2.
- `model_3`: Represents model 3.
"""
@enum ModelType begin
    model_1 = 1
    model_2 = 2
    model_3 = 3
end

"""
    NeighbouringCell(id::Int, face_area::Float64, face_normal::Point2{Float64}, toCenter::Point2{Float64}, transformation::Matrix{Float64})
    
A struct representing a neighbouring cell.

# Arguments
- `id::Int`: The ID of the neighbouring cell.
- `face_area::Float64`: The area of the face shared with the neighbouring cell.
- `face_normal::Point2{Float64}`: The normal vector of the face shared with the neighbouring cell.
- `toCenter::Point2{Float64}`: The vector from the cell center to the center of the neighbouring cell.
- `transformation::Matrix{Float64}`: The transformation matrix from the neighbouring cell to the current cell.

"""
struct NeighbouringCell
    id::Int
    face_area::Float64
    face_normal::Point2{Float64}
    toCenter::Point2{Float64}
    transformation::Matrix{Float64}
end


"""
    struct LcmCell

The `LcmCell` struct represents a cell in the LCM simulation.

# Fields
- `id::Int`: The ID of the cell.
- `vertex_ids::Tuple{Int, Int, Int}`: The IDs of the vertices that define the cell's geometry.
- `vertices::Vector{Point3{Float64}}`: The coordinates of the vertices that define the cell's geometry.
- `center::Point3{Float64}`: The center coordinate of the cell.
- `thickness::Float64`: The thickness of the cell.
- `area::Float64`: The area of the cell.
- `volume::Float64`: The volume of the cell.
- `permeability::Float64`: The permeability of the cell.
- `porosity::Float64`: The porosity of the cell.
- `permeability_Z::Float64`: The permeability in thickness direction the of the cell.
- `thickness_DM::Float64`: The thickness of the DM of the cell.
- `permeability_DM::Float64`: The permeability of the DM of the of the cell.
- `porosity_DM::Float64`: The porosity of the DM of the of the cell.
- `permeability_DM_Z::Float64`: The permeability in thickness direction of the DM of the of the cell.
- `ap::Float64`: The ap value of the cell.
- `bp::Float64`: The bp value of the cell.
- `cp::Float64`: The cp value of the cell.
- `reference_direction::Point3{Float64}`: The reference direction of the cell.
- `alpha::Float64`: The alpha value of the cell.
- `type::CELLTYPE`: The type of the cell.
- `part_id::Int`: The ID of the part that the cell belongs to.
- `num_neighbours::Int`: The number of neighbouring cells.
- `neighbour_ids::Vector{Int}`: The IDs of the neighbouring cells.
- `neighbours::Vector{NeighbouringCell}`: The neighbouring cells.
- `n_CK::Float64`: Cozeny-Karman parameter

# Constructors
- `LcmCell(id::Int, vertex_ids::Tuple{Int, Int, Int}, vertices::Vector{Point3{Float64}}, center::Point3{Float64}, part_id::Int, num_neighbours::Int, neighbour_ids::Vector{Int}, thickness::Float64, area::Float64, volume::Float64, permeability::Float64, porosity::Float64, permeability_Z::Float64, thickness_DM::Float64, permeability_DM::Float64, porosity_DM::Float64, permeability_DM_Z::Float64, ap::Float64, bp::Float64, cp::Float64, alpha::Float64, n_CK::Float64, reference_direction::Point3{Float64}, type::CELLTYPE)`: Constructs a `LcmCell` object with all the required fields.
- `LcmCell(id::Int, vertex_ids::Tuple{Int, Int, Int}, vertices::Vector{Point3{Float64}}, center::Point3{Float64}, part_id::Int, num_neighbours::Int, neighbour_ids::Vector{Int}, thickness::Float64, area::Float64, volume::Float64, permeability::Float64, porosity::Float64, permeability_Z::Float64, thickness_DM::Float64, permeability_DM::Float64, porosity_DM::Float64, permeability_DM_Z::Float64, ap::Float64, bp::Float64, cp::Float64, alpha::Float64, n_CK::Float64, reference_direction::Point3{Float64})`: Constructs a `LcmCell` object without the `type` field.

"""
struct LcmCell
    id::Int
    
    # basic cell geometry
    vertex_ids::Tuple{Int, Int, Int}
    vertices::Vector{Point3{Float64}}
    center::Point3{Float64}

    # cell parameters
    thickness::Float64
    area::Float64
    volume::Float64
    permeability::Float64
    porosity::Float64
    permeability_Z::Float64
    thickness_DM::Float64
    permeability_DM::Float64
    porosity_DM::Float64#
    permeability_DM_Z::Float64
    ap::Float64
    bp::Float64
    cp::Float64
    reference_direction::Point3{Float64}
    alpha::Float64
    n_CK::Float64
    type::CELLTYPE
    part_id::Int

    # neighbour information
    num_neighbours::Int
    neighbour_ids::Vector{Int}
    neighbours::Vector{NeighbouringCell}

    function LcmCell(
        id::Int,
        vertex_ids::Tuple{Int, Int, Int},
        vertices::Vector{Point3{Float64}},
        center::Point3{Float64},
        part_id::Int,
        num_neighbours::Int,
        neighbour_ids::Vector{Int},

        thickness::Float64,
        area::Float64,
        volume::Float64,
        permeability::Float64,
        porosity::Float64,
        permeability_Z::Float64,
        thickness_DM::Float64,
        permeability_DM::Float64,
        porosity_DM::Float64,
        permeability_DM_Z::Float64,
        ap::Float64,
        bp::Float64,
        cp::Float64,

        alpha::Float64,
        n_CK::Float64,
        reference_direction::Point3{Float64},
        type::CELLTYPE,

    )
        @assert num_neighbours == length(neighbour_ids) "Mismatch in num_neighbours and length(neighbour_ids)."

        neighbours = Vector{NeighbouringCell}(undef, num_neighbours)

        new(
            id,
            vertex_ids,
            vertices,
            center,
            thickness,
            area,
            volume,
            permeability,
            porosity,
            permeability_Z,
            thickness_DM,
            permeability_DM,
            porosity_DM,
            permeability_DM_Z,
            ap,
            bp,
            cp,
            reference_direction,
            alpha,
            n_CK,
            type,
            part_id,
            num_neighbours,
            neighbour_ids,
            neighbours
        )
    end

    function LcmCell(
        id::Int,
        vertex_ids::Tuple{Int, Int, Int},
        vertices::Vector{Point3{Float64}},
        center::Point3{Float64},
        part_id::Int,
        num_neighbours::Int,
        neighbour_ids::Vector{Int},

        thickness::Float64,
        area::Float64,
        volume::Float64,
        permeability::Float64,
        porosity::Float64,
        permeability_Z::Float64,
        thickness_DM::Float64,
        permeability_DM::Float64,
        porosity_DM::Float64,
        permeability_DM_Z::Float64,
        ap::Float64,
        bp::Float64,
        cp::Float64,

        alpha::Float64,
        n_CK::Float64,
        reference_direction::Point3{Float64}
    )
        @assert num_neighbours == length(neighbour_ids) "Mismatch in num_neighbours and length(neighbour_ids)."

        neighbours = Vector{NeighbouringCell}(undef, num_neighbours)

        new(
            id,
            vertex_ids,
            vertices,
            center,
            thickness,
            area,
            volume,
            permeability,
            porosity,
            permeability_Z,
            thickness_DM,
            permeability_DM,
            porosity_DM,
            permeability_DM_Z,
            ap,
            bp,
            cp,
            reference_direction,
            alpha,
            n_CK,
            inner::CELLTYPE,
            part_id,
            num_neighbours,
            neighbour_ids,
            neighbours
        )
    end
end

"""
    struct NamedPart

A struct representing a named part.

# Fields
- `name::String`: The name of the part.
- `cell_ids::Vector{Int}`: The cell IDs associated with the part.

"""
struct NamedPart
    name::String
    cell_ids::Vector{Int}
end

"""
    struct LcmMesh

A struct representing a mesh in the LCM simulation.
Note that the ordering of the cells in the `cells` field is important and should 
be consistent with the LcmCell IDs. Else the indices in the `NamedPart` structs 
and the `textile_cell_ids`, `inlet_cell_ids`, and `outlet_cell_ids` fields will 
not be correct, since they refer to the position of the cell in the `cells` field.

# Fields
- `N::Int`: The number of vertices in the mesh.
- `vertices::Vector{Point3{Float64}}`: The vertices of the mesh.
- `cells::Vector{LcmCell}`: The cells of the mesh.
- `textile_cell_ids::Vector{Int}`: The cell IDs of the textile cells in the mesh.
- `inlet_cell_ids::Vector{Int}`: The cell IDs of the inlet cells in the mesh.
- `outlet_cell_ids::Vector{Int}`: The cell IDs of the outlet cells in the mesh.
- `named_parts::Vector{NamedPart}`: The inlets and outlets of the mesh.
- `permeability_ratio::Float64`: The permeability ratio of the mesh.

"""
struct LcmMesh 
    N::Int
    vertices::Vector{Point3{Float64}}
    cells::Vector{LcmCell}
    textile_cell_ids::Vector{Int}
    inlet_cell_ids::Vector{Int}
    outlet_cell_ids::Vector{Int}
    named_parts::Vector{NamedPart}
    permeability_ratio::Float64
end

"""
    abstract type AbstractModel
Abstract type representing a generic model.
"""
abstract type AbstractModel end

"""
    struct Model_1 <: AbstractModel

The `Model_1` struct represents a specific model that implements the `AbstractModel` interface.

# Fields
- `p_a::Float64`: Parameter `p_a`
- `p_init::Float64`: Parameter `p_init`
- `p_ref::Float64`: Parameter `p_ref`
- `rho_a::Float64`: Parameter `rho_a`
- `rho_init::Float64`: Parameter `rho_init`
- `mu_resin::Float64`: Parameter `mu_resin`
- `betat2::Float64`: Parameter `betat2`
- `ap1::Float64`: Parameter `ap1`
- `ap2::Float64`: Parameter `ap2`
- `ap3::Float64`: Parameter `ap3`
- `kappa::Float64`: Parameter `kappa`
- `gamma::Float64`: Parameter `gamma`
"""
struct Model_1 <: AbstractModel
    p_a::Float64
    p_init::Float64
    p_ref::Float64
    rho_a::Float64
    rho_init::Float64
    mu_resin::Float64
    betat2::Float64
    ap1::Float64
    ap2::Float64
    ap3::Float64
    kappa::Float64
    gamma::Float64
end

"""
    abstract type Model_2_3 <: AbstractModel
"""
abstract type Model_2_3 <: AbstractModel end

"""
    struct Model_2 <: Model_2_3

The `Model_2` struct represents a specific model that extends `Model_2_3`.

Fields:
- `p_a::Float64`: Pressure of air
- `p_init::Float64`: Initial pressure
- `p_ref::Float64`: Reference pressure
- `rho_a::Float64`: Density of air
- `rho_init::Float64`: Initial density
- `mu_resin::Float64`: Resin viscosity
- `betat2::Float64`: Beta times squared
- `rho_0_air::Float64`: Reference density of air
- `rho_0_oil::Float64`: Reference density of oil
- `rho_ref::Float64`: Reference density
- `betat2_fac::Float64`: Beta times squared factor
- `exp_val::Float64`: Exponential value

"""
struct Model_2 <: Model_2_3
    p_a::Float64
    p_init::Float64
    p_ref::Float64
    rho_a::Float64
    rho_init::Float64
    mu_resin::Float64
    betat2::Float64
    rho_0_air::Float64
    rho_0_oil::Float64
    rho_ref::Float64
    betat2_fac::Float64
    exp_val::Float64
end

"""
    struct Model_3 <: Model_2_3

The `Model_3` struct represents a specific type of model that extends `Model_2_3`.
It contains the following fields:

- `p_a::Float64`: Parameter a
- `p_init::Float64`: Initial parameter
- `p_ref::Float64`: Reference parameter
- `rho_a::Float64`: Density a
- `rho_init::Float64`: Initial density
- `mu_resin::Float64`: Resin viscosity
- `betat2::Float64`: Beta times 2
- `rho_0_air::Float64`: Reference density of air
- `rho_0_oil::Float64`: Reference density of oil
- `rho_ref::Float64`: Reference density
- `betat2_fac::Float64`: Beta times 2 factor
- `exp_val::Float64`: Exponential value

"""
struct Model_3 <: Model_2_3
    p_a::Float64
    p_init::Float64
    p_ref::Float64
    rho_a::Float64
    rho_init::Float64
    mu_resin::Float64
    betat2::Float64
    rho_0_air::Float64
    rho_0_oil::Float64
    rho_ref::Float64
    betat2_fac::Float64
    exp_val::Float64
end

"""
    struct State

The `State` struct represents the state of a simulation at a given time step.

# Fields
- `t::Float64`: Current time
- `iter::Int`: Current iteration
- `deltat::Float64`: Time step size
- `p::Vector{Float64}`: Pressure values
- `gamma::Vector{Float64}`: Fluid saturation values
- `rho::Vector{Float64}`: Density values
- `u::Vector{Float64}`: Velocity in the x-direction
- `v::Vector{Float64}`: Velocity in the y-direction
- `viscosity::Vector{Float64}`: Fluid viscosity values
- `porosity_times_porosity::Vector{Float64}`: Product of porosity and porosity values
- `h::cavity thickness
- `porostiy::cavity porosity
- `p_DM::Vector{Float64}`: Pressure values in distribution medium
- `gamma_DM::Vector{Float64}`: Fluid saturation values in distribution medium
- `rho_DM::Vector{Float64}`: Density values in distribution medium
- `u_DM::Vector{Float64}`: Velocity in the x-direction in distribution medium
- `v_DM::Vector{Float64}`: Velocity in the y-direction in distribution medium
- `gamma_out::Vector{Float64}`: Combined fluid saturation values in distribution medium and fibrous preform
- `w::Vector{Float64}`: Velocity in the z-direction between distribution medium and fibrous preform
- `filled::Vector{Float64}`: Saves the maximum filling state
- `filled_DM::Vector{Float64}`: Saves the maximum filling state in DM

"""
struct State
    t::Float64
    iter::Int
    deltat::Float64
    p::Vector{Float64}
    gamma::Vector{Float64}
    rho::Vector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
    viscosity::Vector{Float64}
    porosity_times_porosity::Vector{Float64}
    h_out::Vector{Float64}
    porosity_out::Vector{Float64}
    p_DM::Vector{Float64}
    gamma_DM::Vector{Float64}
    rho_DM::Vector{Float64}
    u_DM::Vector{Float64}
    v_DM::Vector{Float64}
    gamma_out::Vector{Float64}    
    w::Vector{Float64}
    filled::Vector{Float64}
    filled_DM::Vector{Float64}
end

"""
    struct ScaledProperties

Describes scaled properties of one cell.

# Fields
- `thickness::Float64`: The thickness of the cell.
- `volume::Float64`: The volume of the cell.
- `face::Vector{Float64}`: The face of the cell.
- `porosity::Float64`: The porosity of the cell.
- `permeability::Float64`: The permeability of the cell.
- `viscosity::Float64`: The viscosity of the cell.
- `alpha::Float64`: The alpha value of the cell.

"""
struct ScaledProperties
    thickness::Float64  
    volume::Float64  
    face::Vector{Float64}
    porosity::Float64
    realporosity::Float64
    permeability::Float64
    viscosity::Float64
    alpha::Float64
end

"""
    struct LcmCase

A struct representing an LcmCase.

# Fields
- `mesh::LcmMesh`: The LcmMesh associated with the case.
- `model::AbstractModel`: The Model associated with the case.
- `state::State`: The State associated with the case.
"""
mutable struct LcmCase
    mesh::LcmMesh
    model::AbstractModel
    state::State
end