function update_rho(
    model::AbstractModel,
    Δt::Float64, 
    props::ScaledProperties,
    rho_old::Float64, 
    F_rho_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64
    
    error("This is an abstract function.")
end 

function update_rho_DM(
    model::AbstractModel,
    Δt::Float64, 
    rho_old::Float64, 
    F_rho_num::Float64, 
    porosity::Float64,
    volume::Float64
    )::Float64
    
    error("This is an abstract function.")
end 

function update_porosity_times_porosity(
    model::AbstractModel,
    props::ScaledProperties
)
    error("This is an abstract function.")
end

function update_p(
    model::AbstractModel,
    mesh::LcmMesh,
    rho_new::Float64,
    filled_new::Float64,
    props::ScaledProperties
    )::Float64
    error("This is an abstract function.")
end

function update_u(
    model::AbstractModel, 
    Δt::Float64,
    props::ScaledProperties,
    ∇p::Point2{Float64},
    rho_old::Float64, 
    rho_new::Float64, 
    F_u_num,
    u_old
    )::Float64
    S_u = -∇p[1]
    return (rho_old * u_old - Δt * F_u_num / props.volume + S_u * Δt) / (rho_new + props.viscosity / props.permeability * Δt)
end

function update_u_DM(
    model::AbstractModel, 
    Δt::Float64,
    ∇p::Point2{Float64},
    rho_old::Float64, 
    rho_new::Float64, 
    F_u_num,
    u_old,
    volume,
    permeability,
    viscosity
    )::Float64
    S_u = -∇p[1]
    return (rho_old * u_old - Δt * F_u_num / volume + S_u * Δt) / (rho_new + viscosity / permeability * Δt)
end

function update_v(
    model::AbstractModel, 
    Δt::Float64,
    props::ScaledProperties,
    ∇p::Point2{Float64},
    rho_old::Float64, 
    rho_new::Float64, 
    F_v_num::Float64, 
    v_old::Float64
    )::Float64
    S_v = -∇p[2]
    return (rho_old * v_old - Δt * F_v_num / props.volume + S_v * Δt) / (rho_new + props.viscosity / (props.permeability * props.alpha) * Δt)
end

function update_v_DM(
    model::AbstractModel, 
    Δt::Float64,
    ∇p::Point2{Float64},
    rho_old::Float64, 
    rho_new::Float64, 
    F_v_num,
    v_old,
    volume,
    permeability,
    viscosity
    )::Float64
    S_v = -∇p[2]
    return (rho_old * v_old - Δt * F_v_num / volume + S_v * Δt) / (rho_new + viscosity / permeability * Δt)
end

function update_gamma(
    model::AbstractModel,
    cell::LcmCell,
    Δt::Float64,
    props::ScaledProperties,
    gamma_old::Float64,
    rho_new::Float64,
    F_gamma_num::Float64,
    F_gamma_num1::Float64
    )::Float64

    error("This is an abstract function.")
end

function update_gamma_DM(
    model::AbstractModel,
    Δt::Float64,
    props::ScaledProperties,
    gamma_old::Float64,
    rho_new::Float64,
    F_gamma_num::Float64,
    F_gamma_num1::Float64
    )::Float64

    error("This is an abstract function.")
end

function update_viscosity(
    model::AbstractModel,
    scaled_properties::ScaledProperties
)
    return scaled_properties.viscosity
end

function scale_properties(
    model::AbstractModel,
    cell::LcmCell,
    p_old::Float64,
    porosity_times_porosity_old::Float64,
    viscosity::Float64,
    iter::Int
    )::ScaledProperties

    error("This is an abstract function.")
end