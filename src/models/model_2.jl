function update_rho(
    model::Model_2,
    Δt::Float64, 
    props::ScaledProperties,
    rho_old::Float64, 
    F_rho_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64,
    aux::Float64
    )::Float64

    #val=props.porosity
    #@info "porosity = $val."

    rho_new = (props.porosity * rho_old - Δt * F_rho_num / props.volume) / props.porosity

    return max(rho_new, 0.0)
end

function update_rho_air(
    model::Model_2,
    Δt::Float64, 
    props::ScaledProperties,
    rho_air_old::Float64, 
    F_rho_air_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    rho_air_new = 0.

    return rho_air_new
end

function update_rho_oil(
    model::Model_2,
    Δt::Float64, 
    props::ScaledProperties,
    rho_oil_old::Float64, 
    F_rho_oil_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    rho_oil_new = 0.

    return rho_oil_new
end

function update_rho_DM(
    model::Model_2,
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
    model::Model_2,
    props::ScaledProperties
)

    return 1.0
end

function scale_properties(
    model::Model_2,
    cell::LcmCell,
    p_old::Float64,
    rho_old::Float64,
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