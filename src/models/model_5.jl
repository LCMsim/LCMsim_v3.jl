function update_rho(
    model::Model_5,
    Δt::Float64, 
    props::ScaledProperties,
    rho_old::Float64, 
    F_rho_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64,
    aux::Float64
    )::Float64

    #val=props.porosity
    #@info "porosity = $val."

    rho_new = aux

    return max(rho_new, 0.0)
end

function update_rho_air(
    model::Model_5,
    Δt::Float64, 
    props::ScaledProperties,
    rho_air_old::Float64, 
    F_rho_air_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    rho_air_new = (props.porosity * rho_air_old - Δt * F_rho_air_num / props.volume) / props.porosity

    return max(rho_air_new, 0.0)
end

function update_rho_oil(
    model::Model_5,
    Δt::Float64, 
    props::ScaledProperties,
    rho_oil_old::Float64, 
    F_rho_oil_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    rho_oil_new = (props.porosity * rho_oil_old - Δt * F_rho_oil_num / props.volume) / props.porosity

    return max(rho_oil_new, 0.0)
end

function update_rho_DM(
    model::Model_5,
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
    model::Model_5,
    props::ScaledProperties
)

    return 1.0
end

function scale_properties(
    model::Model_5,
    cell::LcmCell,
    p_old::Float64,
    rho_old::Float64,
    porosity_times_porosity_old::Float64,
    viscosity::Float64,
    iter::Int
    )::ScaledProperties

    visc_factor=0.01  #air viscosity is two orders of magnitude smaller than oil viscosity
    viscosity_val=model.mu_resin*visc_factor + model.mu_resin*(1-visc_factor)*rho_old/model.rho_0_oil

    faces = [x.face_area for x in cell.neighbours]

    return ScaledProperties(
        cell.thickness,
        cell.volume,
        faces,
        cell.porosity,
        cell.porosity,
        cell.permeability,
        viscosity_val,
        cell.alpha
    )
end

function update_p(
    model::Model_5,
    mesh::LcmMesh,
    rho_new::Float64,
    filled_new::Float64,
    props::ScaledProperties,
    aux_vec::Vector{Float64}
    )::Float64

    E_val=aux_vec[1]
    kappa_air_val=aux_vec[2]
    m_0_oil_val=aux_vec[3]
    m_air_val=max(aux_vec[4],1e-12)
    m_oil_val=max(aux_vec[5],1e-12)

    a_val=-m_oil_val/m_0_oil_val*1/E_val
    b_val=m_oil_val/m_0_oil_val-1+model.p_a/E_val*m_oil_val/m_0_oil_val
    c_val=m_air_val*kappa_air_val/(props.porosity*props.volume)
    
    p_1_val=(-b_val+sqrt(b_val^2-4*a_val*c_val))/(2*a_val)
    p_2_val=(-b_val-sqrt(b_val^2-4*a_val*c_val))/(2*a_val)
    if p_1_val>=0
        p_new = p_1_val
    else
        p_new = p_2_val
    end 

    return p_new
end

function update_gamma(
    model::Model_5,
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