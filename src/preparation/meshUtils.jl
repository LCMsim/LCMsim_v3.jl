using CSV
using HDF5
using LinearAlgebra
using GeometryBasics
using StaticArraysCore
using Random

rng = MersenneTwister(1234)

function seed!(seed::Int)
    global rng = MersenneTwister(seed)
end

function create_LcmMesh(file::HDF5.File)
    cellgridid = read_dataset(file["/mesh"], "cells")
    vertices = read_dataset(file["/mesh"], "vertices")

    part_ids = read_dataset(file["/properties"], "part_id")
    type = read_dataset(file["/properties"], "type")
    thickness = read_dataset(file["/properties"], "thickness")
    porosity = read_dataset(file["/properties"], "porosity")
    permeability = read_dataset(file["/properties"], "permeability")
    alpha = read_dataset(file["/properties"], "alpha")
    n_CK = read_dataset(file["/properties"], "n_CK")
    reference_direction = read_dataset(file["/properties"], "reference_direction")
    volume = read_dataset(file["/properties"], "volume")
    area = read_dataset(file["/properties"], "area")
    ap = read_dataset(file["/properties"], "ap")
    bp = read_dataset(file["/properties"], "bp")
    cp = read_dataset(file["/properties"], "cp")
    unique_part_ids = unique(part_ids)

    names = Dict()
    for pid in unique_part_ids
        try
            name = read_attribute(file["/properties"], string(pid))
            names[pid] = name
        catch e

        end
    end

    N = size(cellgridid)[1]

    # create part_id -> cell_id vector for each part
    partID = Dict()
    for part in unique_part_ids
        partID[part] = Vector{Int32}()
    end
    for i in 1:N
        part_id = part_ids[i]
        push!(partID[part_id], Int32(i))
    end

    # calculate neighbours and centers
    each_cells_neighbours, _ = __get_neighbours(cellgridid)
    vertices = [Point3{Float64}(vertices[i, :]) for i in 1:size(vertices)[1]]
    centers = __calculate_cellcenters(cellgridid, vertices)

    # loop to create all cells with basic information
    cells = Vector{LcmCell}(undef, 0)
    for cid in 1:N
        # get geometry and neighbour information for this cell
        vertex_ids = (cellgridid[cid, 1], cellgridid[cid, 2], cellgridid[cid, 3])
        center = centers[cid]
        neighbours = each_cells_neighbours[cid]
        num_neighbours = Int(length(neighbours))

        # get part information for this cell
        part_id = part_ids[cid]

        this_permeability = permeability[cid]
        this_porosity = porosity[cid]
        this_permeability_Z = permeability_Z[cid]
        this_thickness_DM = thickness_DM[cid]
        this_permeability_DM = permeability_DM[cid]
        this_porosity_DM = porosity_DM[cid]
        this_permeability_DM_Z = permeability_DM_Z[cid]
        this_thickness = thickness[cid]
        this_alpha = alpha[cid]
        this_n_CK = n_CK[cid]
        this_reference_direction = reference_direction[cid]
        this_volume = volume[cid]
        this_type = CELLTYPE(type[cid])
        this_area = area[cid]

        this_ap = ap[cid]
        this_bp = bp[cid]
        this_cp = cp[cid]

        cell_vertices = Vector{Point3{Float64}}(undef, 3)
        cell_vertices .= vertices[cellgridid[cid, :]]
        # create cell struct via basic constructor
        cell = LcmCell(
            cid,
            vertex_ids,
            cell_vertices,
            center, 
            part_id,
            num_neighbours,
            neighbours,
            this_thickness,
            this_area,
            this_volume,
            this_permeability,
            this_porosity,
            this_permeability_Z,
            this_thickness_DM,
            this_permeability_DM,
            this_porosity_DM,
            this_permeability_DM_Z,
            this_ap,
            this_cp,
            this_alpha,
            this_n_CK,
            Point3{Float64}(this_reference_direction),
            this_type
        )

        push!(cells, cell)
    end

    theta = Vector{Float64}(undef, N)
    temp_Tmat = Vector{Matrix{Float64}}(undef, N)
    local_verts = Vector{Vector{Point3{Float64}}}(undef, N)
    for cell in cells     
        _theta, 
        _Tmat, 
        _local_verts = __calculate_local_coordinates(cell)

        theta[cell.id] = _theta
        temp_Tmat[cell.id] = _Tmat
        local_verts[cell.id] = _local_verts
    end

    # loop over cells and add neighbour information
    for cell in cells
        
        for (neighbour_count, neighbour) in enumerate(cells[cell.neighbour_ids])

            face_normal, 
            center_to_center, 
            transformation = __calculate_transformation(
                cell,
                neighbour,
                local_verts[cell.id],
                theta[neighbour.id], 
                temp_Tmat[cell.id],
                vertices
            )

            # calculate face area
            face_area = __calculate_face_area(
                cell,
                neighbour
            )

            # create NeighbouringCell
            neighbouring_cell = NeighbouringCell(
                neighbour.id,
                face_area,
                face_normal,
                center_to_center,
                transformation
            )
            cell.neighbours[neighbour_count] = neighbouring_cell
        end

    end

    # create fields for easy access to inlet/ outlet/ textile cells
    inlet_cells = Vector{Int}(undef, 0)
    outlet_cells = Vector{Int}(undef, 0)
    textile_cells = Vector{Int}(undef, 0)
    named_parts = Vector{NamedPart}(undef, 0)
    for (pid, name) in names
        cell_ids = partID[pid]
        example_cell = cells[cell_ids[1]]
        type = example_cell.type
        
        if type == inlet::CELLTYPE
            append!(inlet_cells, cell_ids)
            push!(named_parts, NamedPart(name, cell_ids))
        elseif type == outlet::CELLTYPE
            append!(outlet_cells, cell_ids)
            push!(named_parts, NamedPart(name, cell_ids))
        else
            append!(textile_cells, cell_ids)
        end
    end

    # calculate permeability ration, which will be added as a general information
    permeability_ratio = __calculate_permeability_ratio(cells, textile_cells)

    mesh = LcmMesh(
        N,
        vertices,
        cells,
        textile_cells,
        inlet_cells,
        outlet_cells,
        named_parts,
        permeability_ratio
    )

    return mesh
end

function create_LcmMesh(meshfile::String, partfile::String)
    # read mesh file, extract basic geometric properties + information about physical groups (base, inlets, outlets, patches)
    cellgridid::Matrix{Int}, vertices::Vector{Point3{Float64}}, part_ids::Vector{Int}, N::Int = parse_mesh(meshfile)
    allowed_part_ids = unique(part_ids)

    # read csv part file, pass part_ids found in the meshfile to ensure the part descriptions match
    parts = parse_partfile(partfile, allowed_part_ids)

    # create part_id -> cell_id vector for each part
    partID = Dict()
    for part in unique(part_ids)
        partID[part] = Vector{Int32}()
    end
    for i in 1:N
        part_id = part_ids[i]
        push!(partID[part_id], Int32(i))
    end

    # calculate neighbours and centers
    each_cells_neighbours, type_wall = __get_neighbours(cellgridid)
    centers = __calculate_cellcenters(cellgridid, vertices)

    # loop to create all cells with basic information
    cells = Vector{LcmCell}(undef, 0)
    for cid in 1:N
        # get geometry and neighbour information for this cell
        vertex_ids = (cellgridid[cid, 1], cellgridid[cid, 2], cellgridid[cid, 3])
        center = centers[cid]
        neighbours = each_cells_neighbours[cid]
        num_neighbours = Int(length(neighbours))

        # get part information for this cell
        part_id = part_ids[cid]
        part_parameters = parts[part_id] 

        permeability = __calculate_permeability(part_parameters[KEY_PERMEABILITY], part_parameters[KEY_PERMEABILITY_NOISE])
        porosity = __calculate_porosity(part_parameters[KEY_POROSITY], part_parameters[KEY_POROSITY_NOISE])
        thickness = part_parameters[KEY_THICKNESS]
        alpha = part_parameters[KEY_ALPHA]
        n_CK = part_parameters[KEY_N_CK]
        reference_direction = part_parameters[KEY_REFERENCE_DIRECTION]

        type = part_parameters[KEY_TYPE]
        if type != inlet::CELLTYPE && type != outlet::CELLTYPE && type_wall[cid] == wall::CELLTYPE
            type = wall::CELLTYPE
        end

        name = part_parameters["name"]
        porosity_1 = part_parameters[KEY_POROSITY_1]
        p_1 = part_parameters[KEY_P_1]
        porosity_2 = part_parameters[KEY_POROSITY_2]
        p_2 = part_parameters[KEY_P_2]

        permeability_Z = part_parameters[KEY_PERMEABILITY_Z]
        thickness_DM = part_parameters[KEY_THICKNESS_DM]
        porosity_DM = part_parameters[KEY_POROSITY_DM]
        permeability_DM = part_parameters[KEY_PERMEABILITY_DM]
        permeability_DM_Z = part_parameters[KEY_PERMEABILITY_DM_Z]

        #ap = porosity
        #bp = 0.
        #cp = (porosity_2 - porosity) / (p_2^2)
        
        #bvec_p=[porosity; porosity_1; porosity_2]
        #Amat_p=[1 0 0; 1 p_1 p_1^2; 1 p_2 p_2^2]
        #xvec_p=inv(Amat_p)*bvec_p   
        #ap=xvec_p[1]
        #bp=xvec_p[2]
        #cp=xvec_p[3]

        #ap = porosity
        #bp = 0.
        #cp = (porosity_2 - porosity) / (p_2^2)
        #bvec_p=[porosity; porosity_1; porosity_2]
        #Amat_p=[1 0 0; 1 p_1 p_1^2; 1 p_2 p_2^2]
        #xvec_p=inv(Amat_p)*bvec_p  
        #ap=0.5*(ap+xvec_p[1])
        #bp=0.5*(bp+xvec_p[2])
        #cp=0.5*(cp+xvec_p[3])

        p0=0e5
        p1=p_1
        p2=p_2
        eps0=porosity
        eps1=porosity_1
        eps2=porosity_2
        ap=eps0
        bp=(eps1-eps0)/(p1-p0) -eps0/p2+eps1/p2-1/p2*(eps2-eps1)/(p2-p1)*p1
        cp=-1/p2*(eps1-eps0)/(p1-p0)+1/p2*(eps2-eps1)/(p2-p1)

        @assert (ap + bp*(1e5)+ cp * (1e5^2)) <= 0.99 "Part " * name * ": Porosity for 1e5 Pa is greater than threshold =" * string(ap +bp*(1e5) + cp * (1e5^2)) * "; Please modify input values."

        area, volume = __calculate_area_and_volume(cid, thickness, cellgridid, vertices)

        cell_vertices = Vector{Point3{Float64}}(undef, 3)
        cell_vertices .= vertices[cellgridid[cid, :]]
        # create cell struct via basic constructor
        cell = LcmCell(
            cid,
            vertex_ids,
            cell_vertices,
            center, 
            part_id,
            num_neighbours,
            neighbours,
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
            alpha,
            n_CK,
            reference_direction,
            type
        )

        push!(cells, cell)
    end

    theta = Vector{Float64}(undef, N)
    temp_Tmat = Vector{Matrix{Float64}}(undef, N)
    local_verts = Vector{Vector{Point3{Float64}}}(undef, N)
    for cell in cells
        
        _theta, 
        _Tmat, 
        _local_verts = __calculate_local_coordinates(cell)

        theta[cell.id] = _theta
        temp_Tmat[cell.id] = _Tmat
        local_verts[cell.id] = _local_verts
    end

    # loop over cells and add neighbour information
    for cell in cells
        
        for (neighbour_count, neighbour) in enumerate(cells[cell.neighbour_ids])

            face_normal, 
            center_to_center, 
            transformation = __calculate_transformation(
                cell,
                neighbour,
                local_verts[cell.id],
                theta[neighbour.id], 
                temp_Tmat[cell.id],
                vertices
            )

            # calculate face area
            face_area = __calculate_face_area(
                cell,
                neighbour
            )

            # create NeighbouringCell
            neighbouring_cell = NeighbouringCell(
                neighbour.id,
                face_area,
                face_normal,
                center_to_center,
                transformation
            )
            cell.neighbours[neighbour_count] = neighbouring_cell
        end

    end

    # create fields for easy access to inlet/ outlet/ textile cells
    inlet_cells = Vector{Int}(undef, 0)
    outlet_cells = Vector{Int}(undef, 0)
    textile_cells = Vector{Int}(undef, 0)
    named_parts = Vector{NamedPart}(undef, 0)
    for (pid, part) in parts
        type = part[KEY_TYPE]
        part_id = part[KEY_PART_ID]
        name = part[KEY_NAME]

        cell_ids = partID[part_id]
        if type == inlet::CELLTYPE
            append!(inlet_cells, cell_ids)
            push!(named_parts, NamedPart(name, cell_ids))
        elseif type == outlet::CELLTYPE
            append!(outlet_cells, cell_ids)
            push!(named_parts, NamedPart(name, cell_ids))
        else
            append!(textile_cells, cell_ids)
        end
    end

    # calculate permeability ration, which will be added as a general information
    permeability_ratio = __calculate_permeability_ratio(cells, textile_cells)

    mesh = LcmMesh(
        N,
        vertices,
        cells,
        textile_cells,
        inlet_cells,
        outlet_cells,
        named_parts,
        permeability_ratio
    )

    return mesh
end

function __get_neighbours(cellgridid::Matrix{Int})::Tuple{Vector{Vector{Int}}, Vector{CELLTYPE}}
    N = size(cellgridid)[1]
    maxnumberofneighbours = 10

    # Find the set with the IDs of the neighbouring cells
    # and identify wall cells

    # initialize celltype with 1 = for every cell
    celltype = Vector{CELLTYPE}(undef, N)
    celltype .= inner::CELLTYPE

    faces = Array{Int32}(undef, 0, 3)   #three columns: sgrid id1, grid id2, cell id
    i = 1
    for ind = 1:N
        i1 = cellgridid[ind, 1]
        i2 = cellgridid[ind, 2]
        faces = vcat(faces, [min(i1, i2) max(i1, i2) ind])
        i = i + 1
        i1 = cellgridid[ind, 2]
        i2 = cellgridid[ind, 3]
        faces = vcat(faces, [min(i1, i2) max(i1, i2) ind])
        i = i + 1
        i1 = cellgridid[ind, 3]
        i2 = cellgridid[ind, 1]
        faces = vcat(faces, [min(i1, i2) max(i1, i2) ind])
        i = i + 1
    end
    facessorted = sortslices(faces, dims=1)
    vals1 = unique(facessorted[:, 1])

    # this must be generalized, currently only hard-coded number of neighbouring cells of a tria is possible
    # all considered cases had <<10 neighbouring cells 
    cellneighboursarray = zeros(Int, N, maxnumberofneighbours) .- 9


    for i in 1:length(vals1)
        inds2 = findall(isequal(vals1[i]), facessorted[:, 1])
        facesdetail_unsorted = facessorted[inds2, 2:3]
        facesdetail = sortslices(facesdetail_unsorted, dims=1)
        for j = 1:size(facesdetail, 1)
            i1 = facesdetail[j, 2]
            inds3 = findall(isequal(facesdetail[j, 1]), facesdetail[:, 1])
            inds4 = findall(!isequal(j), inds3)
            inds5 = inds3[inds4]
            if isempty(inds5)
                celltype[i1] = wall::CELLTYPE
            else
                if j == 1
                    for k in 1:length(inds5)
                        matrixrow = cellneighboursarray[i1, :]
                        indcolumn = findfirst(isequal(-9), matrixrow)
                        if isnothing(indcolumn)
                            println(matrixrow)
                            error("More than 10 neighbours of one tria is not supported \n")
                        else
                            cellneighboursarray[i1, indcolumn] = facesdetail[inds5[k], 2]
                        end
                    end
                else
                    for k in 1:1
                        matrixrow = cellneighboursarray[i1, :]
                        indcolumn = findfirst(isequal(-9), matrixrow)
                        if isnothing(indcolumn)
                            println(matrixrow)
                            error("More than 10 neighbours of one tria is not supported" * "\n")
                        else
                            cellneighboursarray[i1, indcolumn] = facesdetail[inds5[k], 2]
                        end
                    end
                end
            end
        end
    end

    cellneighbours = []
    for ind in 1:N
        valid_indices = cellneighboursarray[ind, :] .!= -9
        push!(cellneighbours, Int.(cellneighboursarray[ind, valid_indices]))
    end 
    return cellneighbours, celltype
end

"""
    Takes a vector of vertices and cells defined by 3 vertex indices and returns a containing the 
    centers of the cells.

    Vertices and centers are represented as Point3{Float64}.
    Cells are represented as a Matrix, where each row contains 3 vertex indices.
"""
function __calculate_cellcenters(cellgridid::Matrix{Int}, vertices::Vector{Point3{Float64}})::Vector{Point3{Float64}}

    N = size(cellgridid)[1]
    # loop to define cell center coordinates in global CS
    centers = Vector{Point3{Float64}}(undef, N)
    for ind in 1:N
        i1 = cellgridid[ind, 1]
        i2 = cellgridid[ind, 2]
        i3 = cellgridid[ind, 3]
        centers[ind] = (vertices[i1] + vertices[i2] + vertices[i3]) / 3
    end
    return centers
end

"""
    Calculates the permeability of a cell given the permeability and the permeability noise.
"""
function __calculate_permeability(permeability::Float64, permeability_noise::Float64)::Float64
    if permeability_noise == 0.0
        return permeability
    end
    return randn(rng, Float64) * permeability_noise * permeability + permeability
end

"""
    Calculates the porosity of a cell given the porosity and the porosity noise.
"""
function __calculate_porosity(porosity::Float64, porosity_noise::Float64)::Float64
    if porosity_noise == 0.0
        return porosity
    end
    return randn(rng, Float64) * porosity_noise * porosity + porosity
end

"""
    Calculates the local coordinates of a cell given the vertices of the cell.
"""
function __calculate_local_coordinates(cell::LcmCell)
    
    b1 = cell.vertices[2] - cell.vertices[1]
    b1 = b1 / norm(b1)

    a2 = cell.vertices[3] - cell.vertices[1]
    a2 = a2 / norm(a2)

    b2 = a2 - dot(b1, a2) * b1 # b2 is orthogonal to b1
    b2 = b2 / norm(b2)
    b3 = cross(b1, b2)

    # new cell CS is given by the projection of the primary direction in local CS

    Tmat = [b1 b2 b3]
    xvec = cell.reference_direction
    r1 = Tmat \ xvec # ref dir in local CS

    # Calculate the angle by which b1 must be rotated about the b3-axis to
    # match r1 via relation rotation matrix Rz[theta]*[1 0 0]'=r1, i.e.
    # cos(theta)=r1[1] & sin(theta)=r1[2]
    theta = atan(r1[2], r1[1]) # atan with 2 arguments equals atan2 in julia

    # Rotation of theta about nvec=b3 to get c1 & c2
    c = cos(theta)
    s = sin(theta)

    b1 = b3 * dot(b3, b1) + c * cross(cross(b3, b1), b3) + s * cross(b3, b1)
    b2 = b3 * dot(b3, b2) + c * cross(cross(b3, b2), b3) + s * cross(b3, b2)
    b3 = b3 * dot(b3, b3) + c * cross(cross(b3, b3), b3) + s * cross(b3, b3) # axis of rotation, not changed

    # Tmat describes 
    Tmat = [b1 b2 b3]

    #transformation of vertices into local CS
    local_vert1 = Point3{Float64}(Tmat \ (cell.vertices[1] - cell.center))
    local_vert2 = Point3{Float64}(Tmat \ (cell.vertices[2] - cell.center))
    local_vert3 = Point3{Float64}(Tmat \ (cell.vertices[3] - cell.center))

    return theta, Tmat, [local_vert1, local_vert2, local_vert3]
end

"""
    Calculates the area and volume of a cell given its thickness and the vertices of the cell.
"""
function __calculate_area_and_volume(cell_id::Int, thickness::Float64, cellgridid::Matrix{Int}, vertices::Vector{Point3{Float64}})::Tuple{Float64, Float64}
    vid1 = cellgridid[cell_id, 1]
    vid2 = cellgridid[cell_id, 2]
    vid3 = cellgridid[cell_id, 3]
    area = 0.5 * norm(cross(vertices[vid2] - vertices[vid1], vertices[vid3] - vertices[vid1]))
    volume = thickness * area 
    return area, volume
end

"""
    Calculates the area of the face between two cells.
"""
function __calculate_face_area(cell::LcmCell, neighbour::LcmCell)
    common_vertices = intersect(cell.vertices, neighbour.vertices)

    return Float64(0.5) * (cell.thickness + neighbour.thickness) * norm(common_vertices[2] - common_vertices[1])

end

"""
    Calculates the transformation matrix for the local coordinate system of a cell to the local coordinate system of a neighbouring cell.
"""
function __calculate_transformation(
    cell::LcmCell,
    neighbour::LcmCell,
    local_vertices::Vector{Point3{Float64}},
    theta_neighbour::Float64,
    Tmat::Matrix{Float64},
    vertices::Vector{Point3{Float64}}
)

    # for every neighbour find the two indices belonging to the boundary
    # face in between; face direction is from smaller to larger index
    # x0..local coordinates of smaller index
    # r0..vector from smaller to larger index in LCS

    common_1, common_2 = findall(x -> x in neighbour.vertex_ids, cell.vertex_ids) 
    
    x0 = local_vertices[common_1]
    r0 = local_vertices[common_2] - local_vertices[common_1]

    # define xvec as the vector between cell centers ind &
    # neighbouring cell center [A]
    # (in GCS) & transform xvec in local coordinates bvec, this gives
    # A in LCS
    # find normal distance from A in LCS to the cell boundary
    # with that cell center A in flat geometry & face normal vector
    # can be defined
    # Fill the cell arrays

    x = [0., 0., 0.] # column vector, P at origin of local CS
    P = x
    lambda = dot(-x0, r0) / dot(r0, r0)
    Q1 = x0 + lambda * r0
    l1 = norm(Q1)
    nvec = Q1
    nvec = nvec / norm(nvec)
    
    cell_face_normal= Point2{Float64}(nvec[1:2])

    xvec = neighbour.center - cell.center # A in global CS
    x = Tmat \ xvec # A in local CS
    A = x
    lambda = dot(x - x0, r0) / dot(r0, r0)
    Q2 = x0 + lambda * r0
    l2 = norm(A - Q2)

    center_to_center = P[1:2] + (Q1[1:2] - P[1:2]) + (Q2[1:2] - Q1[1:2]) + l2 / l1 * (Q1[1:2] - P[1:2])

    # transfromation matrix for (u,v) of neighrouring cells to local
    # coordinate system()

    #construction of the third one in outside normal direction
    #based on the length of the two non-common edges

    common_vid_1 = cell.vertex_ids[common_1]
    common_vid_2 = cell.vertex_ids[common_2]

    local_vertices_neighbour = []
    push!(
        local_vertices_neighbour,
        (
            common_vid_1, 
            Point3{Float64}(
                [
                    local_vertices[common_1][1],
                    local_vertices[common_1][2],
                    0.
                ]
            )
        )
    )

    push!(
        local_vertices_neighbour,
        (
            common_vid_2, 
            Point3{Float64}(
                [
                    local_vertices[common_2][1],
                    local_vertices[common_2][2],
                    0.
                ]
            )
        )
    )

    # find third vertex of neighbouring cell, that isn't shared with cell ind
    uncommon_vid_neighbour = only(setdiff(neighbour.vertex_ids, [cell.vertex_ids[common_1], cell.vertex_ids[common_2]]))

    # position relative to center in GCS
    xvec = vertices[uncommon_vid_neighbour] - cell.center

    # position in LCS
    x = Tmat \ xvec
    A = x # A in local CS
    lambda = dot(x - x0, r0) / dot(r0, r0)
    Q2 = x0 + lambda * r0
    l2 = norm(A - Q2)
    pos = (P + (Q1 - P) + (Q2 - Q1) + l2 / l1 * (Q1 - P))
    push!(
        local_vertices_neighbour,
        (
            uncommon_vid_neighbour,
            Point3{Float64}(
                [
                    pos[1],
                    pos[2],
                    0.
                ]
            )
        )
    )

    # sort vertices by their global id, from smallest to biggest
    sort!(local_vertices_neighbour, by = x -> x[1])

    f1 = local_vertices_neighbour[2][2] - local_vertices_neighbour[1][2]
    f1 = f1 / norm(f1)
    a2 = local_vertices_neighbour[3][2] - local_vertices_neighbour[1][2]
    a2 = a2 / norm(a2)
    f2 = a2 - dot(f1, a2) / dot(f1, f1) * f1
    f2 = f2 / norm(f2)
    f3 = cross(f1, f2)


    nvec = f3
    xvec = f1

    c1 = nvec * dot(nvec, xvec) + cos(theta_neighbour) * cross(cross(nvec, xvec), nvec) + sin(theta_neighbour) * cross(nvec, xvec)
    xvec = f2
    c2 = nvec * dot(nvec, xvec) + cos(theta_neighbour) * cross(cross(nvec, xvec), nvec) + sin(theta_neighbour) * cross(nvec, xvec)
    xvec = f3
    c3 = nvec * dot(nvec, xvec) + cos(theta_neighbour) * cross(cross(nvec, xvec), nvec) + sin(theta_neighbour) * cross(nvec, xvec)
    f1 = c1
    f2 = c2
    f3 = c3

    Tmat=[f1[1] f2[1] f3[1]; f1[2] f2[2] f3[2]; f1[3] f2[3] f3[3]]


    #(u,v)_e=T*(u,v)_f
    transformation = Tmat[1:2, 1:2]


    return cell_face_normal, center_to_center, transformation
end

"""
    Calculates the ratio between the maximum and the minimum permeability of all textile 
    cells (exluding inlets and outlets). The permeabilites in principal celldirection as well as 
    as in second principal direction (aka alpha * permeability) are taken into account.
"""
function __calculate_permeability_ratio(cells::Vector{LcmCell}, textileCells::Vector{Int})
    
    min_perm = PERMEABILITY_MAX
    max_perm = PERMEABILITY_MIN
    for ind in textileCells
        cell = cells[ind]
        k = cell.permeability
        α = cell.alpha
        min_perm = min(min_perm, k, α * k)
        max_perm = max(max_perm, k, α * k)
    end

    @assert min_perm != 0.0 "Minimum permeability is 0.0"
    return max_perm / min_perm
end

function parse_partfile(partfilename::String, allowed_part_ids::Vector{Int})
    partfile = CSV.File(partfilename)

    # SOME ASSERTIONS
    column_names = ["name", "type", "part_id", "thickness", "porosity", "porosity_noise", "permeability", "permeability_noise", "alpha", "refdir1", "refdir2", "refdir3", "porosity_1", "p_1","porosity_2", "p_2","n_CK","permeability_Z","thickness_DM","porosity_DM","permeability_DM","permeability_DM_Z"]
    @assert issetequal(String.(partfile.names), column_names) "Invalid column name in part file!"

    # get parts / ids from partfile
    pids_part = Int8.(partfile["part_id"])

    # assert that every part id in the partfile is also in the meshfile
    @assert all(map(x -> x in allowed_part_ids, pids_part)) "Invalid part id in part file!"

    # check validity of part description file (correct number and valid types)
    @assert all([!all(map(!isequal(String(type)), ["base", "inlet", "outlet", "patch"])) for type in unique(partfile["type"])]) "Invalid type entry in part description file!"


    parts = Dict{Int, Dict}()
    base_parameters = nothing
    for row in 1:partfile.rows
        pid = partfile["part_id"][row]
        # get group for this part

        refdir = [partfile["refdir1"][row] partfile["refdir2"][row] partfile["refdir3"][row]]

        # nomalize refdir vector # TODO don't like this kind of calculation here, separate
        refdir = refdir / sqrt(dot(refdir, refdir))
        refdir = Point3{Float64}(refdir[:])

        # TODO not happy with this logic here
        typestring = partfile["type"][row]
        type = inner::CELLTYPE
        if typestring == "inlet"
            type = inlet::CELLTYPE
        elseif typestring == "outlet"
            type = outlet::CELLTYPE
        end

        part_parameters = Dict(
            KEY_PART_ID => pid,
            KEY_REFERENCE_DIRECTION => refdir,
            KEY_ALPHA => Float64(partfile["alpha"][row]),
            KEY_PERMEABILITY => partfile["permeability"][row],
            KEY_PERMEABILITY_NOISE => partfile["permeability_noise"][row],
            KEY_POROSITY => partfile["porosity"][row],
            KEY_POROSITY_NOISE => partfile["porosity_noise"][row],
            KEY_THICKNESS => partfile["thickness"][row],
            KEY_TYPE => type,
            KEY_POROSITY_1 => partfile["porosity_1"][row],
            KEY_P_1 => partfile["p_1"][row],
            KEY_POROSITY_2 => partfile["porosity_2"][row],
            KEY_P_2 => partfile["p_2"][row],
            KEY_PERMEABILITY_Z => partfile["permeability_Z"][row],
            KEY_THICKNESS_DM => partfile["thickness_DM"][row],
            KEY_POROSITY_DM => partfile["porosity_DM"][row],
            KEY_PERMEABILITY_DM => partfile["permeability_DM"][row],
            KEY_PERMEABILITY_DM_Z => partfile["permeability_DM_Z"][row],
            KEY_N_CK => partfile["n_CK"][row],
            KEY_NAME => partfile["name"][row]
        )
        parts[pid] = part_parameters
    end
    return parts
end


