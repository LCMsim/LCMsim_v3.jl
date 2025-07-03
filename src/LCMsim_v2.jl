module LCMsim_v2

using CSV
using HDF5
using JLD2
using ProgressMeter

include( "constants.jl")
include("structs.jl")
include("preparation/transform_bdf.jl")
include("preparation/parse_mesh.jl")
include("preparation/meshUtils.jl")
include("ioUtils.jl")
include("preparation/simUtils.jl")
include("models/abstract_model.jl")
include("models/model_1.jl")
include("models/model_2_3.jl")
include("models/model_2.jl")
include("models/model_3.jl")

export create, 
create_and_solve,
continue_and_solve,
Verbosity, 
apply_pressure,
solve, 
solve!,
LcmCase, 
State, 
Mesh, 
AbstractModel, 
ModelType, 
save_plottable_mesh, 
save_state, 
save_case, 
load_case,
save_project,
load_project,
create_LcmMesh,
State,
LcmMesh

"""
    create(
        meshfile::String,
        partfile::String,
        simfile::String,
        i_model::ModelType,
        i_advanced::Int64
    )::LcmCase

    Creates a mesh, a model and an initial state from the given files.
    The model is chosen according to the given i_model.

    Returns a LcmCase struct.
"""
function create(
    meshfile::String,
    partfile::String,
    simfile::String,
    i_model::ModelType,
    i_advanced::Int64
)::LcmCase

mesh = create_LcmMesh(meshfile, partfile)
    model = create_SimParameters(mesh, simfile, i_model,i_advanced) 
    state = create_initial_state(mesh, model)

    case = LcmCase(
        mesh,
        model,
        state
    )

    return case
end


"""
    create(
        meshfile::String,
        partfile::String,
        simfile::String,
        i_model::ModelType,
        i_advanced::Int64,
        save_path::String,
        save_binary::Bool=true,
        save_hdf::Bool=true
    )::LcmCase

    Creates a mesh, a model and an initial state from the given files.
    The model is chosen according to the given i_model.
    Saves the mesh and the initial state in hdf format.
    Saves the LcmCase as binary in jld2 format.

    Returns the LcmCase struct.
"""
function create(
    meshfile::String,
    partfile::String,
    simfile::String,
    i_model::ModelType,
    i_advanced::Int64,
    save_path::String,
    save_binary::Bool=true,
    save_hdf::Bool=true
)::LcmCase

    # check if path exists and get paths to data.h5 and data.jld2
    hdf_path, jld2_path = check_path(save_path)

    case = create(
        meshfile,
        partfile,
        simfile,
        i_model,
        i_advanced
    )

    if save_hdf
        save_project(case, [case.state], hdf_path)
    end

    if save_binary
        save_case(case, jld2_path)
    end

    return case
end


"""
    create_and_solve(
        save_path::String,
        meshfile::String,
        partfile::String,
        simfile::String,
        i_model::ModelType,
        i_advanced::Int64,
        t_max::Float64,
        t_step::Float64,
        verbosity=verbose::Verbosity,
        save_binary::Bool=true,
        save_hdf::Bool=true
    )::Nothing

    Creates a mesh, a model and an initial state from the given files.
    The model is chosen according to the given i_model.
    Then solves the problem up to the specified end time.
    Saves results in hdf5 format at the specified time steps.
    If t_step is not given, it is set to t_max, aka the results are only saved at the end.
    Additionally saves the resulting LcmCase as binary in jld2 format.
"""
function create_and_solve(
    save_path::String,
    meshfile::String,
    partfile::String,
    simfile::String,
    i_model::ModelType,
    i_advanced::Int64,
    t_max::Float64,
    t_step=Float64(0.0),
    verbosity=verbose::Verbosity,
    save_binary::Bool=true,
    save_hdf::Bool=true  
)::Nothing

    # check if path exists and get paths to data.h5 and data.jld2
    hdf_path, jld2_path = check_path(save_path)

    # create and save case
    case = create(
        meshfile,
        partfile,
        simfile,
        i_model,
        i_advanced,
        save_path,
        save_binary,
        save_hdf
    )
    # retrieve initial state
    state = case.state

    if save_hdf
        save_project(case, [case.state], hdf_path)
    end

    # if t_step is not given, set it to t_max
    if t_step <= 0.0
        t_step = t_max
    end

    if verbosity == verbose::Verbosity
        log_license()
        # TODO log simulation parameters
    end

    t_next = 0.0
    while t_next < t_max
        t_next += t_step
        state = solve(case.model, case.mesh, state, t_next, verbosity)

        if save_hdf
            # append state to previously created hdf file
            save_state(state, hdf_path)
        end
    end

    # save end_case as jld2
    end_case = LcmCase(
        case.mesh,
        case.model,
        state
    )
    if save_binary
        save_case(end_case, jld2_path)
    end
end


"""
    continue_and_solve(
	    source_path::String,
        save_path::String,
        meshfile::String,
        partfile::String,
        simfile::String,
        i_model::ModelType,
        t_max::Float64,
        output_intervals::Int,
        verbosity=verbose::Verbosity,
        save_binary::Bool=true,
        save_hdf::Bool=true
    )::Nothing

    This function allows to continue a simulation. 
	A new LcmCase is created but the initial values are from the a saved LcmCase
    Creates a mesh, a model and an initial state from the given files.
    The model is chosen according to the given i_model.
    Then solves the problem up to the specified end time.
    Saves results in hdf5 format at the specified output intervals.
    If t_step is not given, it is set to t_max, aka the results are only saved at the end.
    Additionally saves the resulting LcmCase as binary in jld2 format.
"""

function continue_and_solve(
    source_path::String,
	save_path::String,
    meshfile::String,
    partfile::String,
    simfile::String,
    i_model::ModelType,
    i_advanced::Int64,
    t_max::Float64,
    output_intervals=Int(4),
    verbosity=verbose::Verbosity,
    save_binary::Bool=true,
    save_hdf::Bool=true  
)::Nothing

    case_old = load_case(source_path)
    # retrieve old state
    state_old = case_old.state

    # check if path exists and get paths to data.h5 and data.jld2
    hdf_path, jld2_path = check_path(save_path)

    # create and save case
    case = create(
        meshfile,
        partfile,
        simfile,
        i_model,
        i_advanced,
        save_path
    )
    # retrieve initial state
    state = case.state

    # write relevant old values into new state
    state=State(
        state_old.t,
        state_old.iter,
        state_old.deltat,
        state_old.p,
        state_old.gamma,
        state_old.rho,
        state_old.u,
        state_old.v,
        state.viscosity,
        state.porosity_times_porosity,
        state_old.h_out,
        state_old.porosity_out,
        state_old.p_DM,
        state_old.gamma_DM,
        state_old.rho_DM,
        state_old.u_DM,
        state_old.v_DM,
        state_old.gamma_out,
        state_old.w,
        state_old.filled,
        state_old.filled_DM
    )

    #COb: Apply new boundary conditions if changed
    state.p[case.mesh.inlet_cell_ids] .= case.model.p_a
    state.u[case.mesh.inlet_cell_ids] .= U_A
    state.v[case.mesh.inlet_cell_ids] .= V_A
    state.rho[case.mesh.inlet_cell_ids] .= case.model.rho_a
    state.gamma[case.mesh.inlet_cell_ids] .= GAMMA_A
    state.p[case.mesh.outlet_cell_ids] .= case.model.p_init
    state.u[case.mesh.outlet_cell_ids] .= U_INIT
    state.v[case.mesh.outlet_cell_ids] .= V_INIT
    state.rho[case.mesh.outlet_cell_ids] .= case.model.rho_init
    state.gamma[case.mesh.outlet_cell_ids] .= GAMMA_INIT
    state.p_DM[case.mesh.inlet_cell_ids] .= case.model.p_a
    state.u_DM[case.mesh.inlet_cell_ids] .= U_A
    state.v_DM[case.mesh.inlet_cell_ids] .= V_A
    state.rho_DM[case.mesh.inlet_cell_ids] .= case.model.rho_a
    state.gamma_DM[case.mesh.inlet_cell_ids] .= GAMMA_A
    state.p_DM[case.mesh.outlet_cell_ids] .= case.model.p_init
    state.u_DM[case.mesh.outlet_cell_ids] .= U_INIT
    state.v_DM[case.mesh.outlet_cell_ids] .= V_INIT
    state.rho_DM[case.mesh.outlet_cell_ids] .= case.model.rho_init
    state.gamma_DM[case.mesh.outlet_cell_ids] .= GAMMA_INIT

    if save_hdf
        save_plottable_mesh(case.mesh, hdf_path)
        save_state(state, hdf_path)
    end

    # if t_step is not given, set it to t_max
    if output_intervals > 0
        t_step = (ceil(t_max)-floor(state_old.t))/output_intervals
    else
        t_step=t_max
    end

    if verbosity == verbose::Verbosity
        log_license()
        # TODO log simulation parameters
    end

    t_next = state_old.t+0.0
    while t_next < t_max
        t_next += t_step
        state = solve(case.model, case.mesh, state, t_next, verbosity)

        if save_hdf
            # append state to previously created hdf file
            save_state(state, hdf_path)
        end
    end

    # save end_case as jld2
    end_case = LcmCase(
        case.mesh,
        case.model,
        state
    )
    if save_binary
        save_case(end_case, jld2_path)
    end
end

"""
    create_and_solve_gui(
        save_path::String,
        meshfile::String,
        partfile::String,
        simfile::String,
        i_model::ModelType,
        i_advanced::Int64,
        t_max::Float64,
        t_step::Float64,
        verbosity=verbose::Verbosity,
        save_binary::Bool=true,
        save_hdf::Bool=true
    )::Nothing

    Wrapper around create_and_solve
"""
function create_and_solve_gui(
    save_path::String,
    meshfile::String,
    partfile::String,
    simfile::String,
    i_model::ModelType,
    i_advanced::Int64,
    t_max::Float64,
    output_intervals::Int,
    verbosity=verbose::Verbosity,
    save_binary::Bool=true,
    save_hdf::Bool=true  
)::Nothing
    t_step=t_max/output_intervals
    create_and_solve(save_path,meshfile,partfile,simfile,i_model,i_advanced,t_max,t_step,verbosity,save_binary,save_hdf)
end

"""
    continue_and_solve_gui(
	    source_path::String,
        save_path::String,
        meshfile::String,
        partfile::String,
        simfile::String,
        i_model::ModelType,
        i_advanced::Int64,
        t_max::Float64,
        output_intervals::Int,
        verbosity=verbose::Verbosity,
        save_binary::Bool=true,
        save_hdf::Bool=true
    )::Nothing

    Wrapper around continue_and_solve
"""
function continue_and_solve_gui(
    source_path::String,
	save_path::String,
    meshfile::String,
    partfile::String,
    simfile::String,
    i_model::ModelType,
    i_advanced::Int64,
    t_max::Float64,
    output_intervals=Int(4),
    verbosity=verbose::Verbosity,
    save_binary::Bool=true,
    save_hdf::Bool=true  
)::Nothing
    continue_and_solve(source_path,save_path,meshfile,partfile,simfile,i_model,i_advanced,t_max,output_intervals,verbosity,save_binary,save_hdf)
end


"""
    solve(
        source_path::String,
        save_path::String,
        t_max::Float64,
        t_step::Float64,
        verbosity=verbose::Verbosity,
        save_binary::Bool=true,
        save_hdf::Bool=true
    )::Nothing

    This function allows to continue a simulation from a previously in binary saved LcmCase.
    Solves the problem up to the specified end time.
    Saves results in hdf5 format at the specified time steps.
    If t_step is given to be zero, it is set to t_max, aka the results are only saved at the end.
    Additionally saves the resulting LcmCase as binary in jld2 format.
"""
function solve(
    source_path::String,
    save_path::String,
    t_max::Float64,
    t_step=Float64(0.0),
    verbosity=verbose::Verbosity,
    save_binary::Bool=true,
    save_hdf::Bool=true
)::Nothing
    # check if path exists and get paths to data.h5 and data.jld2
    hdf_path, jld2_path = check_path(save_path)

    case = load_case(source_path)
    # retrieve initial state
    state = case.state

    if save_hdf
        save_project(case, [case.state], hdf_path)
    end

    # if t_step is not given, set it to t_max
    if t_step <= 0.0
        t_step = t_max
    end

    t_next = state.t
    while t_next < t_max
        t_next += t_step
        state = solve(case.model, case.mesh, state, t_next, verbosity)

        if save_hdf
            # append state to previously created hdf file
            save_state(state, hdf_path)
        end
    end

    # save end_case as jld2
    end_case = LcmCase(
        case.mesh,
        case.model,
        state
    )
    if save_binary
        save_case(end_case, jld2_path)
    end
end

"""
    solve(
        case::LcmCase,
        t_max::Float64
    )::State

    Solves the problem for the given LcmCase up to the specified end time.
    Returns the new State.
"""
function solve(
    case::LcmCase,
    t_max::Float64,
    verbosity=silent::Verbosity
)::State
    return solve(case.model, case.mesh, case.state, t_max, verbosity)
end

"""
    apply_pressure!(
        case::LcmCase,
        actions::Vector{Tuple{String, Float64}}
    )::LcmCase

    Applies actions consisting of a part name and a pressure value to the given LcmCase. 
    One action is a tuple of a part name and a pressure value.
    For every existing inlet and outlet, an action needs to be provided.
    Mutates the given LcmCase by replacing the state, but does not mutate the old state itself.
    Returns the mutated LcmCase.
"""
function apply_pressure!(
    case::LcmCase,
    actions::Vector{Tuple{String, Float64}}
)::LcmCase
    given_names = sort([action[1] for action in actions])
    existing_names = sort([part.name for part in case.mesh.named_parts])

    @assert given_names == existing_names "Given names do not match existing names."

    state = deepcopy(case.state)

    for action in actions
        part_index = findfirst(part -> part.name == action[1], case.mesh.named_parts)
        part = case.mesh.named_parts[part_index]
        for cell in case.mesh.cells[part.cell_ids]
            state.p[cell.id] = action[2]
        end
    end

    case.state = state

    return case
end

"""
    solve!(
        case::LcmCase,
        t_max::Float64,
        actions::Vector{Tuple{String, Float64}},
        verbosity=silent::Verbosity
    )::State

    Applies pressure actions and solves the problem for the given LcmCase up to the specified end time.
    One action is a tuple of a part name and a pressure value.
    For every existing inlet and outlet, an action needs to be provided.
    Mutates the given LcmCase by replacing the state, but does not mutate the old state itself.
    Returns the new state.
"""
function solve!(
    case::LcmCase,
    t_max::Float64,
    actions::Vector{Tuple{String, Float64}},
    verbosity=silent::Verbosity
)::State
    case = apply_pressure!(case, actions)
    case.state = solve(case, t_max, verbosity)
    return case.state
end

##################################################################################################
# INTERNAL FUNCTIONS
##################################################################################################

function numerical_gradient(i_method::Int, cell::LcmCell, p_old::Vector{Float64})::Point2{Float64}
    if i_method == 1
        #least square solution to determine gradient
        bvec = Vector{Float64}(undef, cell.num_neighbours)
        Amat = Array{Float64}(undef, cell.num_neighbours, 2)
        for (i_neighbour, neighbour) in enumerate(cell.neighbours)
            i_P = cell.id
            i_A = neighbour.id
            Amat[i_neighbour, :] = neighbour.toCenter
            bvec[i_neighbour] = p_old[i_A] - p_old[i_P]
        end

        if cell.num_neighbours > 1
            ∇p = Amat \ bvec
        else
            ∇p = [0.; 0.]
        end
    elseif i_method == 2
        #least square solution to determine gradient with limiter
        bvec = Vector{Float64}(undef, cell.num_neighbours)
        Amat = Array{Float64}(undef, cell.num_neighbours, 2)
        for (i_neighbour, neighbour) in enumerate(cell.neighbours)
            i_P = cell.id
            i_A = neighbour.id
            exp_limiter = 2
            wi = 1 / norm(neighbour.toCenter) ^ exp_limiter
            Amat[i_neighbour, :] = wi .* neighbour.toCenter
            bvec[i_neighbour] = wi * (p_old[i_A] - p_old[i_P])
        end

        if cell.num_neighbours > 1
            ∇p = Amat \ bvec
        else
            ∇p = [0.; 0.]
        end
    elseif i_method == 3
        #least square solution to determine gradient - runtime optimized
        bvec = Vector{Float64}(undef, cell.num_neighbours)
        Amat = Array{Float64}(undef, cell.num_neighbours, 2)
        for (i_neighbour, neighbour) in enumerate(cell.neighbours)
            i_P = cell.id
            i_A = neighbour.id
            Amat[i_neighbour, :] = neighbour.toCenter
            bvec[i_neighbour] = p_old[i_A] - p_old[i_P]
        end

        if cell.num_neighbours > 1
            Aplus = transpose(Amat) * Amat
            a = Aplus[1, 1]
            b = Aplus[1, 2]
            c = Aplus[2, 1]
            d = Aplus[2, 2]
            bvec_mod = transpose(Amat) * bvec
            inv = 1 / (a * d - b * c)
            # 1 / (ad -bc) * [d -b; -c a]
            dpdx = inv * d * bvec_mod[1] - inv * b * bvec_mod[2]
            dpdy = -inv * c * bvec_mod[1] + inv * a * bvec_mod[2]
            ∇p = [dpdx; dpdy]
        else
            ∇p = [0.; 0.]
        end

    end

    return Point2{Float64}(∇p)
end

function numerical_flux_function(i_method, vars_P, vars_A, face_normal, A)::Vector{Float64}
    if i_method == 1
        #first order upwinding
        rho_P = vars_P[1]
        u_P = vars_P[2]
        v_P = vars_P[3]
        gamma_P = vars_P[4]

        rho_A = vars_A[1]
        u_A = vars_A[2]
        v_A = vars_A[3]
        gamma_A = vars_A[4]

        n_dot_rhou1 = dot(face_normal, 0.5 * (rho_P + rho_A) * [0.5 * (u_P + u_A); 0.5 * (v_P + v_A)])
        n_dot_rhou = n_dot_rhou1
        phi = 1
        F_rho_num_add = n_dot_rhou * phi * A
        if n_dot_rhou1 >= 0
            phi = u_P
        else
            phi = u_A
        end
        F_u_num_add = n_dot_rhou * phi * A
        if n_dot_rhou1 >= 0
            phi = v_P
        else
            phi = v_A
        end
        F_v_num_add = n_dot_rhou * phi * A
        n_dot_u = dot(face_normal, [0.5 * (u_P + u_A); 0.5 * (v_P + v_A)])
        if n_dot_rhou1 >= 0
            phi = gamma_P
        else
            phi = gamma_A
        end
        F_gamma_num_add = n_dot_u * phi * A
        phi = 1
        F_gamma_num1_add = n_dot_u * phi * A
    end
    return [F_rho_num_add, F_u_num_add, F_v_num_add, F_gamma_num_add, F_gamma_num1_add]
end

function numerical_flux_function_boundary(i_method, vars_P, vars_A, A, n_dot_u)::Vector{Float64}
    if i_method == 1
        #first order upwinding
        rho_P = vars_P[1]
        u_P = vars_P[2]
        v_P = vars_P[3]
        gamma_P = vars_P[4]
        rho_A = vars_A[1]
        u_A = vars_A[2]
        v_A = vars_A[3]
        gamma_A = vars_A[4]
        n_dot_rhou = n_dot_u * 0.5 * (rho_A + rho_P)
        phi = 1
        F_rho_num_add = n_dot_rhou * phi * A
        if n_dot_u <= 0
            phi = u_A
        else
            phi = u_P
        end
        F_u_num_add = n_dot_rhou * phi * A
        if n_dot_u <= 0
            phi = v_A
        else
            phi = v_P
        end
        F_v_num_add = n_dot_rhou * phi * A
        if n_dot_u <= 0
            phi = gamma_A
        else
            phi = gamma_P
        end
        F_gamma_num_add = n_dot_u * phi * A
        phi = 1
        F_gamma_num1_add = n_dot_u * phi * A
    end
    return [F_rho_num_add, F_u_num_add, F_v_num_add, F_gamma_num_add, F_gamma_num1_add]
end

"""
    aggregate_neighour_flux(
        cell::LcmCell,
        scaled_properties::ScaledProperties,
        mesh::LcmMesh,
        ∇p::Point2{Float64},
        rho_old::Vector{Float64},
        u_old::Vector{Float64},
        v_old::Vector{Float64},
        gamma_old::Vector{Float64}
    )::Tuple{Float64, Float64, Float64, Float64, Float64}

    Aggregates the neighbour flux numerically for the given cell and scaled properties.
    Returns a tuple of the aggregated fluxes:
        F_rho_num, 
        F_u_num, 
        F_v_num, 
        F_gamma_num, 
        F_gamma_num1
"""
function aggregate_neighour_flux(
    cell::LcmCell,
    scaled_properties::ScaledProperties,
    mesh::LcmMesh,
    ∇p::Point2{Float64},
    rho_old::Vector{Float64},
    u_old::Vector{Float64},
    v_old::Vector{Float64},
    gamma_old::Vector{Float64},
)::Tuple{Float64, Float64, Float64, Float64, Float64} 

    # corresponds to [F_rho_num, F_u_num, F_rho_num, F_gamma_num, F_gamma_num1]
    numerical_flux = zeros(Float64, 5)

    for (i_neighbour, neighbour) in enumerate(cell.neighbours)
        # transform old u, v of neighbour into this cell's coordinates
        uvec = neighbour.transformation * [u_old[neighbour.id]; v_old[neighbour.id]]
        u_A = uvec[1]
        v_A = uvec[2]
    
        vars_P = [
            rho_old[cell.id], 
            u_old[cell.id], 
            v_old[cell.id], 
            gamma_old[cell.id]
        ]
        
        vars_A = [
            rho_old[neighbour.id], 
            u_A,
            v_A, 
            gamma_old[neighbour.id]
        ]

        # scaled face area of i-th neighbour
        A = scaled_properties.face[i_neighbour]

        neighbour_cell = mesh.cells[neighbour.id]
        neighbour_type = neighbour_cell.type
        if neighbour_type == inner::CELLTYPE || neighbour_type == wall::CELLTYPE

            numerical_flux .+= numerical_flux_function(1, vars_P, vars_A, neighbour.face_normal, A)
           
        elseif neighbour_type == inlet::CELLTYPE || neighbour_type == outlet::CELLTYPE

            A = A * cell.thickness / (0.5 * (cell.thickness + neighbour_cell.thickness)) # TODO thickness Factor?

            if neighbour_type == outlet::CELLTYPE

                n_dot_u = dot(neighbour.face_normal, [u_old[cell.id]; v_old[cell.id]])

            elseif neighbour_type == inlet::CELLTYPE

                n_dot_u = min(
                    0, 
                    -1 / scaled_properties.viscosity * dot(
                        [scaled_properties.permeability 0; 0 scaled_properties.permeability * cell.alpha] * ∇p, # <=> [k1; k2] .* ∇p
                        neighbour.face_normal
                    )
                )  #inflow according to Darcy's law and no backflow possible
            end
            numerical_flux .+= numerical_flux_function_boundary(1, vars_P, vars_A, A, n_dot_u)
        end
    end

    F_rho_num = numerical_flux[1]
    F_u_num = numerical_flux[2]
    F_v_num = numerical_flux[3]
    F_gamma_num = numerical_flux[4]
    F_gamma_num1 = numerical_flux[5]

    return F_rho_num, F_u_num, F_v_num, F_gamma_num, F_gamma_num1
end

"""
    aggregate_neighour_flux_DM(
        cell::LcmCell,
        scaled_properties::ScaledProperties,
        mesh::LcmMesh,
        ∇p::Point2{Float64},
        rho_old::Vector{Float64},
        u_old::Vector{Float64},
        v_old::Vector{Float64},
        gamma_old::Vector{Float64}
    )::Tuple{Float64, Float64, Float64, Float64, Float64}

    Aggregates the neighbour flux numerically for the given cell and scaled properties for the distribution medium (DM).
    Returns a tuple of the aggregated fluxes:
        F_rho_num, 
        F_u_num, 
        F_v_num, 
        F_gamma_num, 
        F_gamma_num1
"""
function aggregate_neighour_flux_DM(
    cell::LcmCell,
    scaled_properties::ScaledProperties,
    mesh::LcmMesh,
    ∇p::Point2{Float64},
    rho_old::Vector{Float64},
    u_old::Vector{Float64},
    v_old::Vector{Float64},
    gamma_old::Vector{Float64},
)::Tuple{Float64, Float64, Float64, Float64, Float64} 

    # corresponds to [F_rho_num, F_u_num, F_rho_num, F_gamma_num, F_gamma_num1]
    numerical_flux = zeros(Float64, 5)

#TO DO: Move this to models

    #Only if thickness_DM>0.
    if cell.thickness_DM>0.

        for (i_neighbour, neighbour) in enumerate(cell.neighbours)
            # transform old u, v of neighbour into this cell's coordinates
            uvec = neighbour.transformation * [u_old[neighbour.id]; v_old[neighbour.id]]
            u_A = uvec[1]
            v_A = uvec[2]
        
            vars_P = [
                rho_old[cell.id], 
                u_old[cell.id], 
                v_old[cell.id], 
                gamma_old[cell.id]
            ]
            
            vars_A = [
                rho_old[neighbour.id], 
                u_A,
                v_A, 
                gamma_old[neighbour.id]
            ]

            # scaled face area of i-th neighbour
            #A = scaled_properties.face[i_neighbour]  
            A = cell.neighbours[i_neighbour].face_area/cell.thickness*cell.thickness_DM  #DM: change face to thickness_DM


            neighbour_cell = mesh.cells[neighbour.id]
            neighbour_type = neighbour_cell.type
            if neighbour_type == inner::CELLTYPE || neighbour_type == wall::CELLTYPE

                if neighbour_cell.thickness>0.
                    numerical_flux .+= numerical_flux_function(1, vars_P, vars_A, neighbour.face_normal, A)
                end
            
            elseif neighbour_type == inlet::CELLTYPE || neighbour_type == outlet::CELLTYPE

                A = cell.neighbours[i_neighbour].face_area * cell.thickness / (0.5 * (cell.thickness + neighbour_cell.thickness)) 
                A = A / cell.thickness * cell.thickness_DM

                if neighbour_type == outlet::CELLTYPE

                    n_dot_u = dot(neighbour.face_normal, [u_old[cell.id]; v_old[cell.id]])

                elseif neighbour_type == inlet::CELLTYPE
                    #n_dot_u = min(
                    #    0, 
                    #    -1 / scaled_properties.viscosity * dot(
                    #        [scaled_properties.permeability 0; 0 scaled_properties.permeability * cell.alpha] * ∇p,
                    #        neighbour.face_normal
                    #    )
                    #)  #inflow according to Darcy's law and no backflow possible
                    n_dot_u = min(
                        0, 
                        -1 / scaled_properties.viscosity * dot(         #DM: has to be changed if viscosity is different in DM and FP
                            [cell.permeability_DM 0; 0 cell.permeability_DM] * ∇p,    #DM: change twice to cell.permeability_DM
                            neighbour.face_normal
                        )
                    )  #inflow according to Darcy's law and no backflow possible
                end
                numerical_flux .+= numerical_flux_function_boundary(1, vars_P, vars_A, A, n_dot_u)
            end
        end
    end

    F_rho_num = numerical_flux[1]
    F_u_num = numerical_flux[2]
    F_v_num = numerical_flux[3]
    F_gamma_num = numerical_flux[4]
    F_gamma_num1 = numerical_flux[5]

    return F_rho_num, F_u_num, F_v_num, F_gamma_num, F_gamma_num1
end

"""
    adaptive_timestep(
        model::AbstractModel,
        cells::Vector{LcmCell},
        u_new::Vector{Float64},
        v_new::Vector{Float64},
        u_DM_new::Vector{Float64},
        v_DM_new::Vector{Float64},
        w_new::Vector{Float64},
        Δt
    )::Float64

    Calculates the adaptive time step for the given model and cells.
    Takes all cells given in account.

"""
function adaptive_timestep(
    model::AbstractModel,
    cells::Vector{LcmCell},
    u_new::Vector{Float64},
    v_new::Vector{Float64},
    u_DM_new::Vector{Float64},
    v_DM_new::Vector{Float64},
    w_new::Vector{Float64},
    Δt
)::Float64      
    α =0.1;
    min_val = Inf64
    for cell in cells
        if cell.thickness_DM>0.0
            min_val = min(
                min_val, 
                cell.area / (u_DM_new[cell.id]^2 + v_DM_new[cell.id]^2)
            ) 
            min_val = min(
                min_val, 
                cell.area / (u_new[cell.id]^2 + v_new[cell.id]^2)
            ) 
            min_val = min(
                min_val, 
                (0.5*min(cell.thickness,cell.thickness_DM))^2 / (w_new[cell.id]^2)  #transverse velocity based on through-thickness-distance and rate of under-relaxation
            ) 
        else
            min_val = min(
                min_val, 
                cell.area / (u_new[cell.id]^2 + v_new[cell.id]^2)
            ) 
        end
    end
    # removed sqrt from loop, since it should changed the 
    # index of the minimum and only apply it afterwards
    Δt = (1 - α) * Δt + α * model.betat2 * sqrt(min_val)
    return Δt
end

"""
    solve(
        model::AbstractModel,
        mesh::LcmMesh,
        old_state::State,
        t_next::Float64,
        debug_frequency=nothing, 
        fixed_deltat=nothing,
        verbosity=verbose::Verbosity
    )::State

    Internal solve function. Not necessesarily meant for end users, use convenience interface instead.
    Solves the problem for the given model, mesh, old state and next time.
    Returns the new state.
    Allows to set a fixed time step, which will overwrite the adaptive time step.
    Allows to set a debug frequency, which will return the state after the given number of iterations.
"""
function solve(
    model::AbstractModel,
    mesh::LcmMesh,
    old_state::State,
    t_next::Float64,
    verbosity=silent::Verbosity,
    debug_frequency=nothing, 
    fixed_deltat=nothing
)::State

    # retrieve info from old state
    iter = old_state.iter
    t = old_state.t
    if isnothing(fixed_deltat)
        Δt = old_state.deltat
    else
        Δt = fixed_deltat
    end

    # vectors for old state, do deep copy to not mutate old state
    p_old = deepcopy(old_state.p)
    gamma_old = deepcopy(old_state.gamma)
    rho_old = deepcopy(old_state.rho)
    u_old = deepcopy(old_state.u)
    v_old = deepcopy(old_state.v)
    viscosity = deepcopy(old_state.viscosity)
    cellporositytimesporosityfactor_old = deepcopy(old_state.porosity_times_porosity)
    h_old = deepcopy(old_state.h_out)
    porosity_old = deepcopy(old_state.porosity_out)
    p_DM_old = deepcopy(old_state.p_DM)
    gamma_DM_old = deepcopy(old_state.gamma_DM)
    rho_DM_old = deepcopy(old_state.rho_DM)
    u_DM_old = deepcopy(old_state.u_DM)
    v_DM_old = deepcopy(old_state.v_DM)
    gamma_out_old = deepcopy(old_state.gamma_out)
    w_old = deepcopy(old_state.w)
    filled_old = deepcopy(old_state.filled)
    filled_DM_old = deepcopy(old_state.filled_DM)

    # vectors for new state; already here the boundary conditions are set because not changed inside t loop
    p_new = deepcopy(p_old)
    rho_new = deepcopy(rho_old)
    u_new = deepcopy(u_old)
    v_new = deepcopy(v_old)
    gamma_new = deepcopy(gamma_old)
    h_new = deepcopy(h_old)
    porosity_new = deepcopy(porosity_old)
    p_DM_new = deepcopy(p_DM_old)
    rho_DM_new = deepcopy(rho_DM_old)
    u_DM_new = deepcopy(u_DM_old)
    v_DM_new = deepcopy(v_DM_old)
    gamma_DM_new = deepcopy(gamma_DM_old)
    gamma_out_new = deepcopy(gamma_out_old)
    w_new = deepcopy(w_old)
    filled_new = deepcopy(filled_old)
    filled_DM_new = deepcopy(filled_DM_old)

    if verbosity == verbose::Verbosity
        @info "Start solving at t = $t."
        progress = ProgressMeter.Progress(100; desc="Solving", showspeed=true)
    end

    # time evolution
    while t <= t_next
        
        # iterate over cells, that are neither an inlet nor an outlet
        for cell in mesh.cells[mesh.textile_cell_ids]

            # calculate factors for this cell and scale properties (only applies to model 3)
            scaled_properties = scale_properties(
                model, 
                cell, 
                max(min(p_old[cell.id],model.p_a),model.p_init),  #p_old[cell.id], 
                cellporositytimesporosityfactor_old[cell.id], 
                viscosity[cell.id],
                iter
            )

            #Only for output
            h_new[cell.id]=scaled_properties.thickness
            porosity_new[cell.id]=scaled_properties.realporosity

            # calculate pressure gradient
            ∇p = numerical_gradient(3, cell, p_old)  #numerical_gradient(2, cell, p_old)

            # aggregate neighbour flux numerically
            F_rho_num, F_u_num, F_v_num, F_gamma_num, F_gamma_num1 = aggregate_neighour_flux(
                cell,
                scaled_properties,
                mesh,
                ∇p,
                rho_old,
                u_old,
                v_old,
                gamma_old
            )

            # calculate pressure gradient
            ∇p_DM = numerical_gradient(3, cell, p_DM_old)  #numerical_gradient(2, cell, p_DM_old)

            # aggregate neighbour flux numerically
            F_rho_DM_num, F_u_DM_num, F_v_DM_num, F_gamma_DM_num, F_gamma_DM_num1 = aggregate_neighour_flux_DM(
                cell,
                scaled_properties,
                mesh,
                ∇p_DM,
                rho_DM_old,
                u_DM_old,
                v_DM_old,
                gamma_DM_old
            )


            # aggregate thickness flux numerically
            if cell.thickness_DM>0.
                alpha_w=0.1
                t_FP=scaled_properties.thickness
                t_DM=cell.thickness_DM
                Kz_DM=cell.permeability_DM_Z
                if cell.n_CK>0.
                    Kz_FP=cell.permeability_Z*(1 - (1 - scaled_properties.realporosity)) ^ (cell.n_CK+1) / (1 - scaled_properties.realporosity) ^ cell.n_CK
                else
                    Kz_FP=cell.permeability_Z
                end
                Kz=1/((0.5*t_FP)/(0.5*t_FP+0.5*t_DM)*1/Kz_FP + (0.5*t_DM)/(0.5*t_FP+0.5*t_DM)*1/Kz_DM)
                
                dp=p_old[cell.id]-p_DM_old[cell.id]               
                dz=0.25*t_FP   #working solution
                #dz=(0.2+0.6*gamma_old[cell.id])*t_FP  #alternative with variable distance dz
                dpdz=dp/dz
                w_new[cell.id]=(1-alpha_w)*w_new[cell.id]  +(alpha_w)*(-Kz/scaled_properties.viscosity*dpdz)   #with under-relaxation for smoothing
                if w_new[cell.id]>=0.
                    F_rho_num=F_rho_num-(rho_DM_old[cell.id]*w_new[cell.id]*cell.area)  #negative sign because in the update function there is also a negative sign infront of the flux
                    F_rho_DM_num=F_rho_DM_num+(rho_DM_old[cell.id]*w_new[cell.id]*cell.area)
                else
                    F_rho_num=F_rho_num-(rho_old[cell.id]*w_new[cell.id]*cell.area)
                    F_rho_DM_num=F_rho_DM_num+(rho_old[cell.id]*w_new[cell.id]*cell.area)
                end
                
                #Contribution to F_u_num, F_v_num, F_u_DM_num, F_v_DM_num are neglected because convective terms are not relevant in momentum equations
            else
                w_new[cell.id]=0.
            end
            

            # calculate new state
            rho_new[cell.id] = update_rho(
                model,
                Δt,
                scaled_properties,
                rho_old[cell.id],
                F_rho_num,
                cellporositytimesporosityfactor_old[cell.id]
            )

            cellporositytimesporosityfactor_old[cell.id] = update_porosity_times_porosity(
                model,
                scaled_properties
            )

            u_new[cell.id] = update_u(
                model,
                Δt,
                scaled_properties,
                ∇p,
                rho_old[cell.id],
                rho_new[cell.id],
                F_u_num,
                u_old[cell.id]
            )

            v_new[cell.id] = update_v(
                model,
                Δt,
                scaled_properties,
                ∇p,
                rho_old[cell.id],
                rho_new[cell.id],
                F_v_num,
                v_old[cell.id]
            )

            p_new[cell.id] = update_p(
                model,
                mesh, 
                rho_new[cell.id],
                filled_new[cell.id],
                scaled_properties
            )

            gamma_new[cell.id] = update_gamma(
                model,
                cell,
                Δt,
                scaled_properties,
                gamma_old[cell.id],
                rho_new[cell.id],
                F_gamma_num,
                F_gamma_num1
            )

            viscosity[cell.id] = update_viscosity(
                model,
                scaled_properties
            )

            #filled_new[cell.id]=max(filled_new[cell.id],gamma_new[cell.id])
            filled_new[cell.id]=max(filled_DM_new[cell.id],rho_DM_new[cell.id]/model.rho_0_oil)



            if cell.thickness_DM>0.
                # calculate new state for DM
                rho_DM_new[cell.id] = update_rho_DM(
                    model,
                    Δt,
                    rho_DM_old[cell.id],
                    F_rho_DM_num,
                    cell.porosity_DM,
                    cell.area*cell.thickness_DM
                )

                #This is a workaround because the pressure in the DM cells with thickness_DM=0. is set to model.p_a and results during filling in a too high pressure gradient there
                #Make this a pre-calculated flag
                fact_gradp=1
                for (i_neighbour, neighbour) in enumerate(cell.neighbours)
                    neighbour_cell = mesh.cells[neighbour.id]
                    if neighbour_cell.thickness_DM<0.00001
                        fact_gradp=0.
                    end
                end 
                u_DM_new[cell.id] = update_u_DM(
                    model,
                    Δt,
                    fact_gradp*∇p_DM,
                    rho_DM_old[cell.id],
                    rho_DM_new[cell.id],
                    F_u_DM_num,
                    u_DM_old[cell.id],
                    cell.area*cell.thickness_DM,
                    cell.permeability_DM,
                    scaled_properties.viscosity
                )
                v_DM_new[cell.id] = update_v_DM(
                    model,
                    Δt,
                    fact_gradp*∇p_DM,
                    rho_DM_old[cell.id],
                    rho_DM_new[cell.id],
                    F_v_DM_num,
                    v_DM_old[cell.id],
                    cell.area*cell.thickness_DM,
                    cell.permeability_DM,
                    scaled_properties.viscosity
                )
                p_DM_new[cell.id] = update_p(
                    model,
                    mesh, 
                    rho_DM_new[cell.id],
                    filled_DM_new[cell.id],
                    scaled_properties
                )
                gamma_DM_new[cell.id] = update_gamma_DM(
                    model,
                    cell,
                    Δt,
                    scaled_properties,
                    gamma_DM_old[cell.id],
                    rho_DM_new[cell.id],
                    F_gamma_num,
                    F_gamma_num1
                )
                #filled_DM_new[cell.id]=max(filled_DM_new[cell.id],gamma_DM_new[cell.id])
                filled_DM_new[cell.id]=max(filled_DM_new[cell.id],rho_DM_new[cell.id]/model.rho_0_oil)
            else
                rho_DM_new[cell.id]=0.0
                u_DM_new[cell.id]=0.0
                v_DM_new[cell.id]=0.0
                gamma_DM_new[cell.id]=0.0
                p_DM_new[cell.id]=model.p_a  #model.p_init  #p_new[cell.id]  #because is used for DM pressure gradient calculation where DM ends
                filled_DM_new[cell.id]=0.0
            end


            if cell.thickness_DM>0.
                #1.0..DM+FP filled, 0.5..DM filled, 0.0..not filled
                gamma_out_new[cell.id] = 0.5*(gamma_new[cell.id]+gamma_DM_new[cell.id])       
            else
                gamma_out_new[cell.id] = gamma_new[cell.id]
            end


            # update old state of this cell
            u_old[cell.id] = u_new[cell.id] + 0.0
            v_old[cell.id] = v_new[cell.id] + 0.0
            rho_old[cell.id] = rho_new[cell.id] + 0.0
            p_old[cell.id] = p_new[cell.id] + 0.0
            gamma_old[cell.id] = gamma_new[cell.id] + 0.0

            u_DM_old[cell.id] = u_DM_new[cell.id] + 0.0
            v_DM_old[cell.id] = v_DM_new[cell.id] + 0.0
            rho_DM_old[cell.id] = rho_DM_new[cell.id] + 0.0
            p_DM_old[cell.id] = p_DM_new[cell.id] + 0.0
            gamma_DM_old[cell.id] = gamma_DM_new[cell.id] + 0.0
            w_old[cell.id] = w_new[cell.id] + 0.0
        end

        if verbosity == verbose::Verbosity
            percent = round(Int, 100 * (t - old_state.t) / (t_next - old_state.t))
            if percent < 100
                ProgressMeter.update!(progress, percent)
            end
        end

        # Δt calculation
        if isnothing(fixed_deltat)
            Δt = adaptive_timestep(
                model,
                mesh.cells[mesh.textile_cell_ids],
                u_new,
                v_new,
                u_DM_new,
                v_DM_new,
                w_new,
                Δt
            )
            if isnan(Δt) || isinf(Δt)
                @error "Δt is NaN or Inf. Stopping simulation."
                break
            end
        end
        
        iter = iter + 1
        t = t + Δt

        # for debugging
        if !isnothing(debug_frequency) && (iter - 1) % debug_frequency == 0
            break
        end
    end

    if verbosity == verbose::Verbosity
        ProgressMeter.finish!(progress)
        @info "Solving finished at t = $t."
        #@info "h_new = $h_new"
    end

    # create state to return
    return State(
        t,
        iter,
        Δt,
        p_new,
        gamma_new,
        rho_new,
        u_new,
        v_new,
        viscosity,
        cellporositytimesporosityfactor_old,
        h_new,
        porosity_new,
        p_DM_new,
        gamma_DM_new,
        rho_DM_new,
        u_DM_new,
        v_DM_new,
        gamma_out_new,
        w_new,
        filled_new,
        filled_DM_new
    )
end
end # module LcmSimCore
