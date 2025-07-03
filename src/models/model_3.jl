function update_rho(
    model::Model_3,
    Δt::Float64, 
    props::ScaledProperties,
    rho_old::Float64, 
    F_rho_num::Float64, 
    cellporositytimescellporosityfactor_old::Float64
    )::Float64

    #Modified equation for variable porosity
    rho_new = (props.porosity * rho_old - Δt * F_rho_num / props.volume - rho_old * (props.porosity - cellporositytimescellporosityfactor_old)) / props.porosity

    return max(rho_new, 0.0)
end



function update_rho_DM(
    model::Model_3,
    Δt::Float64, 
    rho_old::Float64, 
    F_rho_num::Float64, 
    porosity::Float64,
    volume::Float64
    )::Float64

    #rho_new = 0.0
    rho_new = (porosity * rho_old - Δt * F_rho_num / volume) / porosity

    return max(rho_new, 0.0)
end


function update_porosity_times_porosity(
    model::Model_3,
    props::ScaledProperties
)

    return props.porosity
end

function scale_properties(
    model::Model_3,
    cell::LcmCell,
    p_old::Float64,
    porosity_times_porosity_old::Float64,
    viscosity::Float64,
    iter::Int
    )::ScaledProperties

    n_CK=cell.n_CK

    if n_CK>0.
        porosity_val = cell.ap + cell.bp * (p_old) + cell.cp * (p_old ^ 2)
        target_porosity_val = (1 - cell.porosity) / (1 - porosity_val) * porosity_val
        if iter == 1 # TODO move this to initialization
            target_porosity_val = target_porosity_val
        else
            target_porosity_val = porosity_times_porosity_old + 0.01 * (target_porosity_val - porosity_times_porosity_old)
        end
        porosity_factor = target_porosity_val / cell.porosity
    else
        porosity_val=cell.porosity
        porosity_factor=1.0
    end

    if n_CK>0.
        permeability_factor = (1 - (1 - porosity_val)) ^ (n_CK+1) / (1 - porosity_val) ^ n_CK
    else
        permeability_factor=1.0
    end
    
    thickness_factor = (1-cell.porosity)/(1-porosity_val)


    volume_factor = 1.0
    face_factor = 1.0
    viscosity_factor = 1.0

    faces = [x.face_area for x in cell.neighbours]

    return ScaledProperties(
        cell.thickness * thickness_factor,
        cell.volume * volume_factor,
        faces .* face_factor,
        cell.porosity * porosity_factor,
        porosity_val,
        cell.permeability * permeability_factor,
        viscosity * viscosity_factor,
        cell.alpha
    )
end





