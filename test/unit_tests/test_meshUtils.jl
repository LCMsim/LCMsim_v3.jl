using Test
using HDF5
using GeometryBasics
include("../../src/LCMsim_v2.jl")

eps = 1e-20

@testset "test meshUtils" begin
    # load correct results
    HyperMesh_file = "test/unit_tests/test_inputs/permeameter1_HyperMesh.dat"
    Gmsh_file = "test/unit_tests/test_inputs/permeameter1_Gmsh.bdf"
    results_file = "test/unit_tests/true_results/permeameter1.h5" 
    cellgridid = Matrix{Int}(undef, 0, 3)
    vertices = Vector{Point3{Float64}}(undef, 0)
    centers = Vector{Point3{Float64}}(undef, 0)
    volume = Vector{Float64}(undef, 0)
    thickness = Vector{Float64}(undef, 0)
    center_to_center = Vector{Vector{Point2{Float64}}}(undef, 0)
    face_area = Vector{Vector{Float64}}(undef, 0)
    face_normal = Vector{Vector{Point2{Float64}}}(undef, 0)
    neighbours = Vector{Vector{Int}}(undef, 0)
    Tmat = Vector{Vector{Matrix{Float64}}}(undef, 0)
    N = 0
    h5open(results_file, "r") do fid
        cellgridid = read_dataset(fid["/mesh/cells"], "cellgridid")
        N = read_attribute(fid["/mesh/cells"], "N")

        grid = read_dataset(fid["/mesh/nodes"], "grid")
        cellcenters = read_dataset(fid["/mesh/aux"], "cellcenters")
        cellcenter_to_cellcenter = read_dataset(fid["/mesh/aux"], "cellcentertocellcenter")
        cellfacenormal = read_dataset(fid["/mesh/aux"], "cellfacenormal")
        cellneighboursarray = read_dataset(fid["/mesh/aux"], "cellneighboursarray")
        cellfacearea = read_dataset(fid["/mesh/aux"], "cellfacearea")
        T11 = read_dataset(fid["/mesh/aux"], "T11")
        T12 = read_dataset(fid["/mesh/aux"], "T12")
        T21 = read_dataset(fid["/mesh/aux"], "T21")
        T22 = read_dataset(fid["/mesh/aux"], "T22")

        for i in 1:size(grid)[1]
            push!(vertices, Point3{Float64}(grid[i, :]))
        end

        
        for i in 1:N
            push!(centers, Point3{Float64}(cellcenters[i, :]))

            this_center_to_center = Vector{Point2{Float64}}(undef, 0)
            this_neighbours = Vector{Int}(undef, 0)
            this_face_normal = Vector{Point2{Float64}}(undef, 0)
            this_face_area = Vector{Float64}(undef, 0)
            this_Tmat = Vector{Matrix{Float64}}(undef, 0)
            for j in 1:10
                if cellneighboursarray[i, j] != -9.
                    push!(this_center_to_center, Point2{Float64}(cellcenter_to_cellcenter[i, j, :]))
                    push!(this_neighbours, cellneighboursarray[i, j])
                    push!(this_face_normal, Point2{Float64}(cellfacenormal[i, j, :]))
                    push!(this_face_area, cellfacearea[i, j])
                    push!(this_Tmat,
                        [
                            T11[i, j] T12[i, j];
                            T21[i, j] T22[i, j]
                        ]
                    )
                end
            end
            push!(neighbours, this_neighbours)
            push!(center_to_center, this_center_to_center)
            push!(face_normal, this_face_normal)
            push!(face_area, this_face_area)
            push!(Tmat, this_Tmat)
        end

        volume = read_dataset(fid["/mesh/aux"], "cellvolume")
        thickness = read_dataset(fid["/mesh/parameters"], "cellthickness")
    end
    
    @test begin
        test_centers = LCMsim_v2.__calculate_cellcenters(cellgridid, vertices)
        all(test_centers .== centers)
    end

    @test begin
        pass = true
        max = 0.
        ind = 0
        for i in 1:N
            test_area, test_volume = LCMsim_v2.__calculate_area_and_volume(i, thickness[i], cellgridid, vertices)
            pass = pass && abs(test_volume - volume[i]) < eps
            if max < abs(test_volume - volume[i])
                max = test_volume - volume[i]
                ind = i
            end
        end
        @info "Maximum volume error: " * string((ind, max))
        pass
    end

    @test begin
        pass = true
        max = 0.
        ind = 0
        
        for cid in 1:N

            cell = LCMsim_v2.LcmCell(
                    cid, 
                    (cellgridid[cid, 1] , cellgridid[cid, 2], cellgridid[cid, 3]),
                    vertices[cellgridid[cid, :]],
                    Point3{Float64}([0., 0., 0.]),
                    0,
                    1,
                    [0],
                    thickness[cid],
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    Point3{Float64}([0., 0., 0.])
                )

            for (i, neighbour) in enumerate(neighbours[cid])
                neighbour_cell = LCMsim_v2.LcmCell(
                    neighbour, 
                    (cellgridid[neighbour, 1] , cellgridid[neighbour, 2], cellgridid[neighbour, 3]),
                    vertices[cellgridid[neighbour, :]],
                    Point3{Float64}([0., 0., 0.]),
                    0,
                    1,
                    [0],
                    thickness[neighbour],
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    Point3{Float64}([0., 0., 0.]),
                )
                
                test_face_area = LCMsim_v2.__calculate_face_area(
                    cell, 
                    neighbour_cell
                )

                err = abs(face_area[cid][i] - test_face_area)
                if err > max
                    max = err
                    ind = (cid, i)
                end
            end
        end

        @info "Maximum face_area error: " * string((ind, max))
        pass
    end

    @test begin
        pass = true
        max = 0.
        ind = 0
        
        theta = Vector{Float64}(undef, N)
        temp_Tmat = Vector{Matrix{Float64}}(undef, N)
        local_verts = Vector{Vector{Point3{Float64}}}(undef, N)
        for cid in 1:N
            cell = LCMsim_v2.LcmCell(
                    cid, 
                    (cellgridid[cid, 1] , cellgridid[cid, 2], cellgridid[cid, 3]),
                    vertices[cellgridid[cid, :]],
                    centers[cid],
                    0,
                    1,
                    [0],
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    Point3{Float64}([1., 0., 0.])
                )

            _theta, _Tmat, _local_verts = LCMsim_v2.__calculate_local_coordinates(cell)
            theta[cid] = _theta
            temp_Tmat[cid] = _Tmat
            local_verts[cid] = _local_verts
        end

        for cid in 1:N
            for (i, neighbour) in enumerate(neighbours[cid])
                cell = LCMsim_v2.LcmCell(
                    cid, 
                    (cellgridid[cid, 1] , cellgridid[cid, 2], cellgridid[cid, 3]),
                    vertices[cellgridid[cid, :]],
                    centers[cid],
                    0,
                    1,
                    [0],
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    Point3{Float64}([1., 0., 0.])
                )

                neighbour_cell = LCMsim_v2.LcmCell(
                    neighbour, 
                    (cellgridid[neighbour, 1] , cellgridid[neighbour, 2], cellgridid[neighbour, 3]),
                    vertices[cellgridid[neighbour, :]],
                    centers[neighbour],
                    0,
                    1,
                    [0],
                    thickness[neighbour],
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    0.,
                    Point3{Float64}([1., 0., 0.]),
                )

                (  
                    test_face_normal, 
                    test_center_to_center, 
                    transformation
                ) = LCMsim_v2.__calculate_transformation(
                    cell,
                    neighbour_cell, 
                    local_verts[cid], 
                    theta[neighbour], 
                    temp_Tmat[cid], 
                    vertices
                )
                err = maximum([
                    norm(face_normal[cid][i] - test_face_normal),
                    norm(center_to_center[cid][i] - test_center_to_center),
                    maximum(abs.(Tmat[cid][i] .- transformation))
                ])

                if err > max
                    max = err
                    ind = (cid, neighbour)
                end
            end
        end

        @info "Maximum transformation error: " * string((ind, max))
        pass
    end

end