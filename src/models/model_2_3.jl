function update_rho(
    model::Model_2_3,
    Δt::Float64, 
    props::ScaledProperties,
    rho_old::Float64, 
    F_rho_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    error("This is an abstract function.")
end

function update_rho_DM(
    model::Model_2_3,
    Δt::Float64, 
    rho_old::Float64, 
    F_rho_num::Float64, 
    porosity::Float64,
    volume::Float64
    )::Float64
    
    error("This is an abstract function.")
end 

function update_porosity_times_porosity(
    model::Model_2_3,
    props::ScaledProperties
)
    error("This is an abstract function.")
end

function update_p(
    model::Model_2_3,
    mesh::LcmMesh,
    rho_new::Float64,
    filled_new::Float64,
    props::ScaledProperties
)::Float64

    #if filled_new>0.8
    #    a_val = model.p_init
    #    c_val = (model.p_a - model.p_init) / (model.rho_0_oil - 0.5*model.rho_0_oil)^model.exp_val
    #    p_new = a_val + c_val * (rho_new - 0.5*model.rho_0_oil)^model.exp_val
    #    p_new = min(1.0 * model.p_a, p_new)
    #    p_new = max(model.p_init, p_new)
    #else
        a_val = model.p_init
        c_val = (model.p_a - model.p_init) / (model.rho_0_oil - model.rho_0_air)^model.exp_val
        p_new = a_val + c_val * (rho_new - model.rho_0_air)^model.exp_val
        p_new = min(1.0 * model.p_a, p_new)
        p_new = max(model.p_init, p_new)
    #end

    return p_new
end

function update_gamma(
    model::Model_2_3,
    cell::LcmCell,
    Δt::Float64,
    props::ScaledProperties,
    gamma_old::Float64,
    rho_new::Float64,
    F_gamma_num::Float64,
    F_gamma_num1::Float64
)
    if cell.thickness_DM>0.
        if rho_new>= 0.9 * model.rho_0_oil  #account for filling gradient in thickness direction in main preform for VARI with DM
            gamma_new = 1.
        else
            gamma_new = 0.
        end
    else
        if rho_new>= 0.5 * model.rho_0_oil
            gamma_new = 1.
        else
            gamma_new = 0.
        end
    end
    #gamma_new=min(max(rho_new/model.rho_0_oil,0.0),1.0)

    return gamma_new
end

function update_gamma_DM(
    model::Model_2_3,
    cell::LcmCell,
    Δt::Float64,
    props::ScaledProperties,
    gamma_old::Float64,
    rho_new::Float64,
    F_gamma_num::Float64,
    F_gamma_num1::Float64
)

    if rho_new>= 0.5 * model.rho_0_oil
        gamma_new = 1.
    else
        gamma_new = 0.
    end
    #gamma_new=min(max(rho_new/model.rho_0_oil,0.0),1.0)

    return gamma_new
end

function scale_properties(
    model::Model_2_3,
    cell::LcmCell,
    p_old::Float64,
    porosity_times_porosity_old::Float64,
    viscosity::Float64,
    iter::Int
    )::ScaledProperties

    error("This is an abstract function.")
end