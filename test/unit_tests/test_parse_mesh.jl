using Test
using HDF5
using GeometryBasics
using LinearAlgebra

include("../../src/LCMsim_v2.jl")


@testset "Test mesh parsing" begin
    # load correct results
    HyperMesh_file = "test/unit_tests/test_inputs/permeameter1_HyperMesh.dat"
    Gmsh_file = "test/unit_tests/test_inputs/permeameter1_Gmsh.bdf"
    inp_file = "test/unit_tests/test_inputs/mesh_permeameter2_byPrePoMax.inp"
    results_file = "test/unit_tests/true_results/permeameter1.h5" 
    results2_file = "test/unit_tests/true_results/permeameter2.h5"
    cellgridid = Matrix{Int}(undef, 0, 3)
    grid = Matrix{Float64}(undef, 0, 3)
    # pids = Vector{Int}(undef, 0) # TODO can't be tested because it's not really supported in original code
    N = 0
    h5open(results_file, "r") do fid
        cellgridid = read_dataset(fid["/mesh/cells"], "cellgridid")
        grid = read_dataset(fid["/mesh/nodes"], "grid")
        # pids = read_dataset(fid["/mesh/aux"], "part_ids")
        N = read_attribute(fid["/mesh/cells"], "N")
    end
    vertices = Vector{Point3{Float64}}(undef, size(grid)[1])
    for i in 1:size(grid)[1]
        vertices[i] = Point3{Float64}(grid[i, :])
    end

    N2 = 0
    cellgridid2 = Matrix{Int}(undef, 0, 3)
    grid2 = Matrix{Float64}(undef, 0, 3)
    h5open(results2_file, "r") do fid
        cellgridid2 = read_dataset(fid["/mesh/cells"], "cellgridid")
        grid2 = read_dataset(fid["/mesh/nodes"], "grid")
        # pids = read_dataset(fid["/mesh/aux"], "part_ids")
        N2 = read_attribute(fid["/mesh/cells"], "N")
    end
    vertices2 = Vector{Point3{Float64}}(undef, size(grid2)[1])
    for i in 1:size(grid2)[1]
        vertices2[i] = Point3{Float64}(grid2[i, :])
    end

    # HyperMesh format
    @test begin
        test_cellgridid, test_vertices, test_pids, test_N = LCMsim_v2.__parse_HyperMesh(HyperMesh_file)

        # test expression
        all(test_cellgridid .== cellgridid) && all(isapprox.(test_vertices, vertices, atol=0.)) && test_N == N
    end

    # Gmsh format
    @test begin
        test_cellgridid, test_vertices, test_pids, test_N = LCMsim_v2.__parse_gmsh(Gmsh_file)
    
        # test expression
        all(test_cellgridid .== cellgridid) && all(isapprox.(test_vertices, vertices, atol=1e-5)) && test_N == N
    end

    # inp format / PrePoMax/ Abaqus
    @test begin
        test_cellgridid, test_vertices, test_pids, test_N = LCMsim_v2.__parse_abaqus(inp_file)

    
        # test expression
        all(test_cellgridid .== cellgridid2) && all(isapprox.(test_vertices, vertices2, atol=0.)) && test_N == N2
    end
end