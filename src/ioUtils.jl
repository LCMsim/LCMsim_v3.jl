using HDF5

"""
    Writes a list of (name, data) tuples as attributes to 'parent'. 
    Caller is responsible to provide a valid parent and attributes.
"""
function write_attributes(parent::Union{HDF5.File,HDF5.Group}, attributes)::Nothing
    for (name, attr) in attributes
        write_attribute(parent, name, attr)
    end
end

"""
    write_datasets(parent::Union{HDF5.File,HDF5.Group}, datasets)::Nothing

    Writes a list of (name, data) tuples as datasets to 'parent'. 
    Caller is responsible to provide a valid parent and data.
"""
function write_datasets(parent::Union{HDF5.File,HDF5.Group}, datasets)::Nothing
    for (name, arr) in datasets
        write_dataset(parent, name, arr)
    end
end

"""
    saveLcmMesh(mesh::LcmMesh, filename::String)

    Saves a LcmMesh struct instance in hdf-format.
"""
function saveLcmMesh(mesh::LcmMesh, filename::String)
    # collect properties from lcm mesh struct
    neighbours = Matrix{Int}(undef, mesh.N, 3)
    centers = Matrix{Float64}(undef, mesh.N, 3)
    cellvolume = Vector{Float64}(undef, mesh.N)
    T11 = Matrix{Float64}(undef, mesh.N, 3)
    T12 = Matrix{Float64}(undef, mesh.N, 3)
    T21 = Matrix{Float64}(undef, mesh.N, 3)
    T22 = Matrix{Float64}(undef, mesh.N, 3)
    cellfacearea = Matrix{Float64}(undef, mesh.N, 3)
    cellfacenormal = Array{Float64}(undef, mesh.N, 3, 2)
    cellcentertocellcenter = Array{Float64}(undef, mesh.N, 3, 2)
    for (cid, cell) in enumerate(mesh.cells)
        neighbours[cid, 1:cell.num_neighbours] = cell.neighbour_ids
        centers[cid, :] = cell.center[1:3]

        for nid in 1:cell.num_neighbours
            T11[cid, nid] = cell.neighbours[nid].transformation[1, 1]
            T12[cid, nid] = cell.neighbours[nid].transformation[1, 2]
            T21[cid, nid] = cell.neighbours[nid].transformation[2, 1]
            T22[cid, nid] = cell.neighbours[nid].transformation[2, 2]
            cellfacearea[cid, nid] = cell.neighbours[nid].face_area
            cellfacenormal[cid, nid, :] = cell.neighbours[nid].face_normal
            cellcentertocellcenter[cid, nid, :] = cell.neighbours[nid].toCenter
        end
        cellvolume[cid] = cell.volume
    end
    # write data to structured format
    h5open(filename, "w") do fid
        mesh = create_group(fid, GROUP_MESH)

        aux = create_group(mesh, "aux")
        write_dataset(aux, "cellvolume", cellvolume)
        write_dataset(aux, "T11", T11)
        write_dataset(aux, "T12", T12)
        write_dataset(aux, "T21", T21)
        write_dataset(aux, "T22", T22)
        write_dataset(aux, "cellcenters", centers)
        write_dataset(aux, "cellcentertocellcenter", cellcentertocellcenter)
        write_dataset(aux, "cellfacearea", cellfacearea)
        write_dataset(aux, "cellfacenormal", cellfacenormal)
        write_dataset(aux, "cellneighboursarray", neighbours)



    end
end

"""
    save_plottable_mesh(mesh::LcmMesh, filename::String)::Nothing

    Saves the information needed to plot a mesh in hdf-format.
"""
function save_plottable_mesh(mesh::LcmMesh, filename::String)::Nothing
    cells = Matrix{Int}(undef, mesh.N, 3)

    # properties that could be of interest to plot
    type = Vector{Int}(undef, mesh.N)
    part_id = Vector{Int}(undef, mesh.N)
    thickness = Vector{Float64}(undef, mesh.N)
    permeability = Vector{Float64}(undef, mesh.N)
    porosity = Vector{Float64}(undef, mesh.N)
    volume = Vector{Float64}(undef, mesh.N)
    area = Vector{Float64}(undef, mesh.N)
    alpha = Vector{Float64}(undef, mesh.N)
    reference_direction = Matrix{Float64}(undef, mesh.N, 3)
    ap = Vector{Float64}(undef, mesh.N)
    bp = Vector{Float64}(undef, mesh.N)
    cp = Vector{Float64}(undef, mesh.N)

    for (cid, cell) in enumerate(mesh.cells)
        cells[cid, :] .= cell.vertex_ids
        type[cid] = Integer(cell.type)
        part_id[cid] = cell.part_id
        thickness[cid] = cell.thickness
        permeability[cid] = cell.permeability
        porosity[cid] = cell.porosity
        volume[cid] = cell.volume
        area[cid] = cell.area
        alpha[cid] = cell.alpha
        reference_direction[cid, :] = cell.reference_direction
        ap[cid] = cell.ap
        bp[cid] = cell.bp
        cp[cid] = cell.cp
    end

    M = length(mesh.vertices)
    vertices = Matrix{Float64}(undef, M, 3)
    for (vid, vertex) in enumerate(mesh.vertices)
        vertices[vid, :] = vertex
    end

    # collect name, part_id tuples
    parts = []
    for part in mesh.named_parts
        name = part.name
        part_id_ = mesh.cells[part.cell_ids[1]].part_id # not very good, data structure misses part_id
        push!(parts, (name, part_id_))
    end

    # write data to structured format
    h5open(filename, "w") do fid
        mesh = create_group(fid, "mesh")
        write_dataset(mesh, "vertices", vertices)
        write_dataset(mesh, "cells", cells)

        props = create_group(fid, "properties")
        write_dataset(props, "volume", volume)
        write_dataset(props, "area", area)
        write_dataset(props, "permeability", permeability)
        write_dataset(props, "porosity", porosity)
        write_dataset(props, "thickness", thickness)
        write_dataset(props, "type", type)
        write_dataset(props, "part_id", part_id)
        write_dataset(props, "alpha", alpha)
        write_dataset(props, "reference_direction", reference_direction)
        write_dataset(props, "ap", ap)
        write_dataset(props, "bp", bp)
        write_dataset(props, "cp", cp)

        for (name, part_id) in parts
            write_attribute(props, string(part_id), name)
        end
        # for base part
        write_attribute(props, "1", "base")
    end
end

# functions to load and save models

"""
    save_model(model::Union{AbstractModel,Model_2_3}, filename::String)::Nothing

    Abstract method, will throw error.
    Concrete implementations are provided for Model_1, Model_2 and Model_3.
"""
function save_model(model::Union{AbstractModel,Model_2_3}, filename::String)::Nothing
    error("This is an abstract function.")
end

"""
    load_model(filename::String)::AbstractModel

    Loads model from a hdf5 file.
    Decides which model to load based on the model_type attribute.
"""
function load_model(filename::String)::AbstractModel
    model = nothing

    h5open(filename, "r") do f
        type = read_attribute(f, "model_type")
        if type == 1
            group = f["model"]
            ap1 = read_attribute(group, "ap1")
            ap2 = read_attribute(group, "ap2")
            ap3 = read_attribute(group, "ap3")
            betat2 = read_attribute(group, "betat2")
            gamma = read_attribute(group, "gamma")
            kappa = read_attribute(group, "kappa")
            p_a = read_attribute(group, "p_a")
            p_init = read_attribute(group, "p_init")
            p_ref = read_attribute(group, "p_ref")
            rho_a = read_attribute(group, "rho_a")
            rho_init = read_attribute(group, "rho_init")
            mu_resin = read_attribute(group, "mu_resin")


            model = Model_1(
                p_a,
                p_init,
                p_ref,
                rho_a,
                rho_init,
                mu_resin,
                betat2,
                ap1,
                ap2,
                ap3,
                kappa,
                gamma
            )
        elseif type == 2
            group = f["model"]
            rho_0_air = read_attribute(group, "rho_0_air")
            rho_0_oil = read_attribute(group, "rho_0_oil")
            betat2_fac = read_attribute(group, "betat2_fac")
            betat2 = read_attribute(group, "betat2")
            exp_val = read_attribute(group, "exp_val")
            p_a = read_attribute(group, "p_a")
            p_init = read_attribute(group, "p_init")
            p_ref = read_attribute(group, "p_ref")
            rho_a = read_attribute(group, "rho_a")
            rho_init = read_attribute(group, "rho_init")
            rho_ref = read_attribute(group, "rho_ref")
            mu_resin = read_attribute(group, "mu_resin")

            model = Model_2(
                p_a,
                p_init,
                p_ref,
                rho_a,
                rho_init,
                mu_resin,
                betat2,
                rho_0_air,
                rho_0_oil,
                rho_ref,
                betat2_fac,
                exp_val
            )

        elseif type == 3
            group = f["model"]
            rho_0_air = read_attribute(group, "rho_0_air")
            rho_0_oil = read_attribute(group, "rho_0_oil")
            betat2_fac = read_attribute(group, "betat2_fac")
            betat2 = read_attribute(group, "betat2")
            exp_val = read_attribute(group, "exp_val")
            p_a = read_attribute(group, "p_a")
            p_init = read_attribute(group, "p_init")
            p_ref = read_attribute(group, "p_ref")
            rho_a = read_attribute(group, "rho_a")
            rho_init = read_attribute(group, "rho_init")
            rho_ref = read_attribute(group, "rho_ref")
            mu_resin = read_attribute(group, "mu_resin")

            model = Model_3(
                p_a,
                p_init,
                p_ref,
                rho_a,
                rho_init,
                mu_resin,
                betat2,
                rho_0_air,
                rho_0_oil,
                rho_ref,
                betat2_fac,
                exp_val
            )
        end
    end

    return model
end

"""
    save_model(model::Model_1, filename::String)::Nothing

    Saves a Model_1 object to a hdf5 file.
"""
function save_model(model::Model_1, filename::String)::Nothing
    h5open(filename, "r+") do f
        write_attribute(f, "model_type", 1)
        group = create_group(f, "model")
        write_attribute(group, "ap1", model.ap1)
        write_attribute(group, "ap2", model.ap2)
        write_attribute(group, "ap3", model.ap3)
        write_attribute(group, "betat2", model.betat2,)
        write_attribute(group, "gamma", model.gamma)
        write_attribute(group, "kappa", model.kappa)
        write_attribute(group, "p_init", model.p_init)
        write_attribute(group, "p_a", model.p_a,)
        write_attribute(group, "p_ref", model.p_ref)
        write_attribute(group, "rho_a", model.rho_a)
        write_attribute(group, "rho_init", model.rho_init)
        write_attribute(group, "mu_resin", model.mu_resin)
        write_dataset(group, "dummy", [1])
    end
end

"""
    save_model(model::Model_2, filename::String)::Nothing

    Saves a Model_2 object to a hdf5 file.
"""
function save_model(model::Model_2, filename::String)::Nothing
    h5open(filename, "r+") do f
        write_attribute(f, "model_type", 2) 
        group = create_group(f, "model")
        write_attribute(group, "rho_0_air", model.rho_0_air)
        write_attribute(group, "rho_0_oil", model.rho_0_oil)
        write_attribute(group, "betat2_fac", model.betat2_fac)
        write_attribute(group, "betat2", model.betat2)
        write_attribute(group, "exp_val", model.exp_val)
        write_attribute(group, "p_a", model.p_a)
        write_attribute(group, "p_init", model.p_init)
        write_attribute(group, "p_ref", model.p_ref)
        write_attribute(group, "rho_a", model.rho_a)
        write_attribute(group, "rho_init", model.rho_init)
        write_attribute(group, "rho_ref", model.rho_ref)
        write_attribute(group, "mu_resin", model.mu_resin)
    end
end

"""
    save_model(model::Model_3, filename::String)::Nothing

    Saves a Model_3 object to a hdf5 file.
"""
function save_model(model::Model_3, filename::String)::Nothing
    h5open(filename, "r+") do f
        write_attribute(f, "model_type", 3)
        group = create_group(f, "model")
        write_attribute(group, "rho_0_air", model.rho_0_air)
        write_attribute(group, "rho_0_oil", model.rho_0_oil)
        write_attribute(group, "betat2_fac", model.betat2_fac)
        write_attribute(group, "betat2", model.betat2)
        write_attribute(group, "exp_val", model.exp_val)
        write_attribute(group, "p_a", model.p_a)
        write_attribute(group, "p_init", model.p_init)
        write_attribute(group, "p_ref", model.p_ref)
        write_attribute(group, "rho_a", model.rho_a)
        write_attribute(group, "rho_init", model.rho_init)
        write_attribute(group, "rho_ref", model.rho_ref)
        write_attribute(group, "mu_resin", model.mu_resin)
    end
end

"""
    save_state(
        state::State,
        file::String,
    )::Nothing

    Requires an existing hdf5 file and a state struct.
    Creates a group named ("state%09i", iter) at meshfile["/"] and
    writes state attributes/ datasets to it.

    Example: iter = 1 -> "state000000001"
"""
function save_state(
    state::State,
    fid::HDF5.File,
)::Nothing
    statestring = @sprintf("state%09i", state.iter)
    state_group = create_group(fid, statestring)

    write_attributes(
        state_group,
        [
            (HDF_T, state.t),
            (HDF_ITER, state.iter),
            (HDF_DELTAT, state.deltat)
        ]
    )
    write_datasets(
        state_group,
        [
            (HDF_P, state.p),
            (HDF_RHO, state.rho),
            (HDF_U, state.u),
            (HDF_V, state.v),
            (HDF_GAMMA, state.gamma),
            (HDF_VISCOSITY, state.viscosity),
            (HDF_CELLPOROSITYTIMESCELLPOROSITY_FACTOR, state.porosity_times_porosity),
            (HDF_H_OUT, state.h_out),
            (HDF_POROSITY_OUT, state.porosity_out),
            (HDF_P_DM, state.p_DM),
            (HDF_RHO_DM, state.rho_DM),
            (HDF_U_DM, state.u_DM),
            (HDF_V_DM, state.v_DM),
            (HDF_GAMMA_DM, state.gamma_DM),
            (HDF_GAMMA_OUT, state.gamma_out),
            (HDF_W, state.w),
            (HDF_FILLED, state.filled),
            (HDF_FILLED_DM, state.filled_DM)
        ]
    )
end


"""
    save_state(
    state::State,
    filename::String
)::Nothing

    Saves a state to an existing hdf5 file.
"""
function save_state(
    state::State,
    filename::String
)::Nothing
    h5open(filename, "r+") do f
        save_state(state, f)
    end
end

"""
    save_states(
        states::Vector{State},
        filename::String
    )::Nothing

    Saves a vector of states to an existing hdf5 file.
"""
function save_states(
    states::Vector{State},
    filename::String
)
    h5open(filename, "r+") do f
        for state in states
            save_state(state, f)
        end
    end
end

# function to load all states in a hdf file
function load_states(filename::String)::Vector{State}
    states = []
    h5open(filename, "r") do f
        for key in keys(f["/"])
            g = f[key]
            if occursin("state", String(key))
                t = read_attribute(g, "t")

                state = State(
                    t,
                    read_attribute(g, "iter"),
                    read_attribute(g, "deltat"),
                    read_dataset(g, "p"),
                    read_dataset(g, "gamma"),
                    read_dataset(g, "rho"),
                    read_dataset(g, "u"),
                    read_dataset(g, "v"),
                    read_dataset(g, "viscosity"),
                    read_dataset(g, "cellporositytimesporosityfactor")
                )
                push!(states, state)
            end
        end

        sort!(states, by=x -> x.t)
    end
    return states
end

"""
    check_path(path::String)::Tuple{String, String}

    Checks if the given path exists and returns the path to the hdf5 and jld2 file.
"""
function check_path(path::String)::Tuple{String,String}
    # assert that path exists
    @assert isdir(path) "The given path does not exist."

    # check if path contains trailing slash
    if path[end] == '/'
        return path * "data.h5", path * "data.jld2"
    else
        return path * "/data.h5", path * "/data.jld2"
    end
end

"""
    load_case(path::String)::LcmCase

    Loads a LcmCase object from a jld2 file.
"""
function load_case(path::String)::LcmCase
    # assert that path exists
    @assert isfile(path) "The given path does not exist."  #COb: path already includes the filename

    # load the jld2 file
    jld2file = JLD2.load(path)

    # extract the LcmCase object
    return jld2file["LcmCase"]
end

"""
    save_case(case::LcmCase, path::String)::Nothing

    Saves a LcmCase object to a jld2 file.
"""
function save_case(case::LcmCase, path::String)::Nothing
    ## assert that path exists
    #@assert isdir(path) "The given path does not exist."  
    ## save the LcmCase object
    #JLD2.save(path * "/data.jld2", "LcmCase", case)

    #COb: path already includes the filename
    @assert isfile(path) "The given path does not exist." 
    JLD2.save(path, "LcmCase", case)
end

# funtion to log license and version
function log_license()
    @info """
    LCMsim v2 
    
    LCMsim v2 is Julia code which simulates the mold filling in Liquid Composite Molding (LCM) 
    manufacturing process. 
    
    Contributors:
    - Christof Obertscheider (FHWN): Conceptual design
    - Leonard Heber (ISSE): Software engineering
    - Ewald Fauster (MUL): Scientific advisor

    This program is freely available and open source, licensed under the GNU General Public 
    License version 3 (GPL v3).
    
    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    See the GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License along with this program. If not, 
    see http://www.gnu.org/licenses/. 

    This software is free of charge and may be used for commercial and academic purposes. Please mention 
    the use of this software at an appropriate place in your work.

    Submit bug reports to christof.obertscheider@fhwn.ac.at
    """
end

"""
    load_project(
        filename::String
    )::Tuple{LcmCase,Vector{State}}

    Loads a project from a hdf5 file and returns a LcmCase object and a vector with all states.
    The LcmCase object contains the last state of the state vector.
"""
function load_project(filename::String)::Tuple{LcmCase,Vector{State}}
    model = load_model(filename)
    mesh = nothing
    h5open(filename, "r") do f
        mesh = create_LcmMesh(f)
    end
    states = load_states(filename)

    case = LcmCase(mesh, model, last(states))

    return case, states
end

"""
    save_project(
    case::LcmCase,
    states::Vector{State},
    filename::String
)

    Saves project ro hdf-file. Assumes the state of the case to be included in the states vector.
"""
function save_project(
    case::LcmCase,
    states::Vector{State},
    filename::String
)
    save_plottable_mesh(case.mesh, filename)
    save_model(case.model, filename)

    save_states(states, filename)

end