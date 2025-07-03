function update_rho(
    model::Model_1,
    Δt::Float64, 
    props::ScaledProperties,
    rho_old::Float64, 
    F_rho_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    rho_new =  rho_old - Δt * F_rho_num / props.volume

    return max(rho_new, 0.0)
end

function update_rho_DM(
    model::Model_1,
    Δt::Float64, 
    rho_old::Float64, 
    F_rho_num::Float64, 
    porosity::Float64,
    volume::Float64
    )::Float64

    rho_new = 0.0

    return max(rho_new, 0.0)
end

function update_porosity_times_porosity(
    model::Model_1,
    props::ScaledProperties
)

    return 1.0
end

function update_p(
    model::Model_1,
    mesh::LcmMesh,
    rho_new::Float64,
    filled_new::Float64,
    props::ScaledProperties
    )::Float64

    #EOS:
    if model.gamma > 1.01
        return model.ap1 * rho_new^2 + model.ap2 * rho_new + model.ap3
    else
        return model.kappa * rho_new^model.gamma
    end
end


function update_gamma(
    model::Model_1,
    cell::LcmCell,
    Δt::Float64,
    props::ScaledProperties,
    gamma_old::Float64,
    rho_new::Float64,
    F_gamma_num::Float64,
    F_gamma_num1::Float64
    )::Float64

    gamma_new = (props.porosity * gamma_old - Δt * (F_gamma_num - gamma_old * F_gamma_num1) / props.volume) / props.porosity

    return clamp(gamma_new, 0., 1.)
end

function update_gamma_DM(
    model::Model_1,
    cell::LcmCell,
    Δt::Float64,
    props::ScaledProperties,
    gamma_old::Float64,
    rho_new::Float64,
    F_gamma_num::Float64,
    F_gamma_num1::Float64
    )::Float64

    gamma_new = (props.porosity * gamma_old - Δt * (F_gamma_num - gamma_old * F_gamma_num1) / props.volume) / props.porosity

    return clamp(gamma_new, 0., 1.)
end


function scale_properties(
    model::Model_1,
    cell::LcmCell,
    p_old::Float64,
    porosity_times_porosity_old::Float64,
    viscosity::Float64,
    iter::Int
    )::ScaledProperties

    faces = [x.face_area for x in cell.neighbours]

    return ScaledProperties(
        cell.thickness,
        cell.volume,
        faces,
        cell.porosity,
        cell.porosity,
        cell.permeability,
        viscosity,
        cell.alpha
    )
end