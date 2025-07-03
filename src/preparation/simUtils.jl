using CSV
using Printf

"""
    Calculates physical parameters according to method 1.

    Returns (all Float64)
        p_a
        p_init
        p_ref
        rho_a
        rho_init
        kappa 
        ap1 
        ap2 
        ap3
"""
function calculate_physical_parameters_method_1(
    p_a::Float64, 
    p_init::Float64,
    p_ref::Float64, 
    rho_ref::Float64, 
    gamma::Float64
)
    #Normalization for Delta p: p->p-p_init
    p_eps = Float64(0.001e5) #Float64(0.000e5);  #
    p_a = p_a - p_init + p_eps
    p_init = p_init - p_init + p_eps
    p_ref = p_ref  #p_ref-p_init+p_eps;
    kappa = p_ref / (rho_ref^gamma)
    #Lookuptable for adiabatic law (required for stability)
    p_int1 = Float64(0.0e5)
    rho_int1 = (p_int1 / kappa)^(1 / gamma)
    p_int3 = Float64(0.5e5)
    rho_int3 = (p_int3 / kappa)^(1 / gamma)
    p_int4 = Float64(1.0e5)
    rho_int4 = (p_int4 / kappa)^(1 / gamma)
    A = [rho_int1^2 rho_int1 Float64(1.0); rho_int3^2 rho_int3 Float64(1.0); rho_int4^2 rho_int4 Float64(1.0)]
    b = [p_int1; p_int3; p_int4]
    apvals = A \ b
    ap1 = apvals[1]
    ap2 = apvals[2]
    ap3 = apvals[3]
    rho_a = (p_a / kappa)^(Float64(1) / gamma)
    rho_init = (p_init / kappa)^(Float64(1) / gamma)

    return p_a, p_init, p_ref, rho_a, rho_init, kappa, ap1, ap2, ap3
end

"""
    Calculates physical parameters according to methods 2 & 3.

    Returns (all Float64)
        rho_a
        rho_init
"""
function calculate_physical_parameters_method_2_3(
    rho_0_air::Float64, 
    rho_0_oil::Float64
)
    return rho_0_oil, rho_0_air
end

function parse_paramter_file(filename::String)
    sim_params = CSV.File(filename)
    column_names = ["p_ref", "rho_ref", "gamma", "mu_resin", "p_a", "p_init", "rho_0_air", "rho_0_oil"]
    @assert issetequal(String.(sim_params.names), column_names) "Invalid column name in parameter file!"
    @assert sim_params.rows == 1 "More than 1 row in parameter file not allowed."

    p_ref = sim_params["p_ref"][1]
    rho_ref = sim_params["rho_ref"][1]
    gamma = sim_params["gamma"][1]
    mu_resin = sim_params["mu_resin"][1]
    p_a = sim_params["p_a"][1]
    p_init = sim_params["p_init"][1]
    rho_0_air = float(sim_params["rho_0_air"][1])
    rho_0_oil = float(sim_params["rho_0_oil"][1])

    return (p_a,
    p_init,
    p_ref,
    rho_ref,
    rho_0_air,
    rho_0_oil,
    gamma,
    mu_resin)
end

"""
    Creates a Model from parameters given via 'parameter_file' and 
    according to 'i_model'.
"""
function create_SimParameters(
    mesh::LcmMesh, 
    parameter_file::String, 
    i_model::ModelType,
    i_advanced::Int64
)::AbstractModel
    (p_a,
    p_init,
    p_ref,
    rho_ref,
    rho_0_air,
    rho_0_oil,
    gamma,
    mu_resin) = parse_paramter_file(parameter_file)

    if i_model == model_1::ModelType
        betat2 = 0.1
        p_a, p_init, p_ref, rho_a, rho_init, kappa, ap1, ap2, ap3 = calculate_physical_parameters_method_1(p_a, p_init, p_ref, rho_ref, gamma)
        return Model_1(
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
    
    else
        rho_a, rho_init = calculate_physical_parameters_method_2_3(rho_0_air, rho_0_oil)

        betat2_fac = 1.
        exp_val = 1.

        # TODO kann permeability_ratio > 100 sein? nicht > 1.0?
        #if mesh.permeability_ratio >= 100
        #    betat2_fac = 0.01  #1.0  #0.25  #
        #    exp_val = 25  #4;  #10;  #
        #else
        #    betat2_fac = 1.0  #0.1  #0.25  #
        #    exp_val = 4  #25;  #10;  #
        #end
        #COb: Linear functions of permeabiliity ratio

        #exp_val=floor(4+(25-4)/(1250-1)*(mesh.permeability_ratio-1))
        #exp_val=25  #4  #
        #betat2_fac=0.1  #1.0+(0.01-1.0)/(1250-1)*(mesh.permeability_ratio-1)

        if i_advanced==0
            exp_val=4 
            betat2_fac=1.0
        elseif i_advanced==1
            exp_val=25
            betat2_fac=0.1
        elseif i_advanced==2
            exp_val=25
            betat2_fac=0.01
        end

        #@info "betat2_fac = $betat2_fac"  #COb
        #@info "exp_val = $exp_val"  #COb 
        betat2= 0.1 * betat2_fac 

        if i_model == model_2::ModelType
            return Model_2(
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
        else
            return Model_3(
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
end

function calculate_initial_timestep(
    mesh::LcmMesh,
    model::AbstractModel
)::Float64

    min_area = Inf64
    max_velocity = -Inf64

    for cell in mesh.cells
        min_area = min(
            min_area, 
            cell.area
        )
        # currently this calculation uses a global viscosity (mu_resin)
        max_velocity = max(
            max_velocity, 
            cell.permeability / model.mu_resin, 
            cell.alpha * cell.permeability / model.mu_resin
        ) 
    end
    # TODO Sollte man hier nicht area und velocity einer Zelle nehmen?

    max_velocity *= (model.p_a - model.p_init) / min_area
        
    Δt = sqrt(min_area) / max_velocity

    return Δt
end

function create_initial_state(
    mesh::LcmMesh, 
    model::AbstractModel
)::State
    # read some mesh properties needed to create arrays
    inlet_cells = mesh.inlet_cell_ids

    # initialize 
    p = zeros(Float64, mesh.N) .+ model.p_init
    rho = zeros(Float64, mesh.N) .+ model.rho_init
    u = zeros(Float64, mesh.N) .+ U_INIT
    v = zeros(Float64, mesh.N) .+ V_INIT
    gamma = zeros(Float64, mesh.N) .+ GAMMA_INIT
    viscosity = zeros(Float64, mesh.N) .+ model.mu_resin
    h_out = zeros(Float64, mesh.N)  #is only used as output
    porosity_out = zeros(Float64, mesh.N)  #is only used as output
    p_DM = zeros(Float64, mesh.N) .+ model.p_init
    rho_DM = zeros(Float64, mesh.N) .+ model.rho_init
    u_DM = zeros(Float64, mesh.N) .+ U_INIT
    v_DM = zeros(Float64, mesh.N) .+ V_INIT
    gamma_DM = zeros(Float64, mesh.N) .+ GAMMA_INIT
    gamma_out = zeros(Float64, mesh.N)  #is only used as output
    w = zeros(Float64, mesh.N)  #is only used as output
    filled = zeros(Float64, mesh.N)  #is only used as output
    filled_DM = zeros(Float64, mesh.N)  #is only used as output

    # set boundary conditions of inlet cells
    p[inlet_cells] .= model.p_a
    rho[inlet_cells] .= model.rho_a
    u[inlet_cells] .= U_A
    v[inlet_cells] .= V_A
    gamma[inlet_cells] .= GAMMA_A
    p_DM[inlet_cells] .= model.p_a
    rho_DM[inlet_cells] .= model.rho_a
    u_DM[inlet_cells] .= U_A
    v_DM[inlet_cells] .= V_A
    gamma_DM[inlet_cells] .= GAMMA_A
    

    cellporositytimesporosityfactor = Vector{Float64}(undef, mesh.N)
    for (cid, cell) in enumerate(mesh.cells)
        porosity_val = cell.ap + cell.bp * (p[cid]) + cell.cp * (p[cid] ^ 2)
        target_porosity_val = (1 - cell.porosity) / (1 - porosity_val) * porosity_val
        cellporositytimesporosityfactor[cid] = target_porosity_val
    end
    t = 0.0
    iter = 1
    deltat = calculate_initial_timestep(mesh, model)

    return State(
        t,
        iter,
        deltat,
        p,
        gamma,
        rho,
        u,
        v,
        viscosity,
        cellporositytimesporosityfactor,
        h_out,
        porosity_out,
        p_DM,
        gamma_DM,
        rho_DM,
        u_DM,
        v_DM,
        gamma_out,
        w,
        filled,
        filled_DM
    )
end
