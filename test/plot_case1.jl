using Gtk4Makie; using Gtk4
using GtkObservables
using HDF5
using JLD2
using GeometryBasics
using NativeFileDialog
using LinearAlgebra

i_var=i_var

guipath = savepath

screen = Gtk4Makie.GTKScreen(title="LCMsim v3")    
filename=joinpath(guipath,"start.jld2")
@load filename xyz N cgammavec minval maxval
p1=display(screen, poly(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(minval,maxval)))
ax1=current_axis()  #perspectiveness=0.5, viewmode=:fitzoom,aspect=(ax,ay,az))
fig1=current_figure()
ax1.show_axis=false  
h=grid(screen)
g=Gtk4.GtkGrid()

polys::Vector{Polygon} = []
gamma::Vector{Float64} = []
h5open(joinpath(savepath,"data.h5"), "r") do meshfile

    grid = read_dataset(meshfile["/mesh"], "vertices")
    cells = read_dataset(meshfile["/mesh"], "cells")
    N = size(cells)[1]
    xvec=Array{Float64}(undef, 3, N)
    yvec=Array{Float64}(undef, 3, N)
    zvec=Array{Float64}(undef, 3, N)
    cgammavec=Array{Float64}(undef, 3, N)
    cgammavec_bw=Array{Float64}(undef, 3, N)


    for cid in 1:N
        points::Vector{Point2f} = []
        for j in 1:3
            gid = cells[cid, j]
            x = grid[gid, 1]
            y = grid[gid, 2]
            z = grid[gid, 3]
            push!(points, Point2f(x, y))            
            xvec[j,cid]=x
            yvec[j,cid]=y
            zvec[j,cid]=z
            cgammavec[j,cid]=1.0
            cgammavec_bw[j,cid]=1.0
        end
        push!(polys, Polygon(points))
    end
    xyz = reshape([xvec[:] yvec[:] zvec[:]]', :)

    states = keys(meshfile["/"])
    states = filter(s -> s[1:3] == "sta", states)
    
    if i_var==0
        gammas = [read_dataset(meshfile["/" * state], "gamma") for state in states]
    elseif i_var==1
        gammas = [read_dataset(meshfile["/" * state], "h_out") for state in states]
    elseif i_var==2
        gammas = [read_dataset(meshfile["/" * state], "porosity_out") for state in states]        
    elseif i_var==3
        gammas = [read_dataset(meshfile["/" * state], "gamma_DM") for state in states]       
    elseif i_var==4
        gammas = [read_dataset(meshfile["/" * state], "gamma_out") for state in states]      
    elseif i_var==5
        gammas = [read_dataset(meshfile["/" * state], "w") for state in states]           
    elseif i_var==6
        gammas = [read_dataset(meshfile["/" * state], "p") for state in states]           
    end
    #gammas = [read_dataset(meshfile["/" * state], "p") for state in states]

    gamma = read_dataset(meshfile["properties"], "thickness")

    cgammasvec=Array{Float64}(undef, 3, N,length(gammas))
    cgammasvec_bw=Array{Float64}(undef, 3, N,length(gammas))
    for tid in 1:length(gammas)
        for cid in 1:N
            for j in 1:3
                cgammasvec[j,cid,tid]=gammas[tid][cid]

                #if cgammasvec[j,cid,tid]>=960;
                #    cgammasvec[j,cid,tid]=960
                #end

                if cgammasvec[j,cid,tid]>=0.8;
                    cgammasvec_bw[j,cid,tid]=1.0
                else
                    cgammasvec_bw[j,cid,tid]=0.0
                end
            end
        end
    end
    cgammavec=cgammasvec[:,:,length(gammas)]
    cgammavec_bw=cgammasvec_bw[:,:,length(gammas)]

    times = [read_attribute(meshfile["/" * state], "t") for state in states]
    
    #bounding box
    deltax=maximum(xvec)-minimum(xvec)
    deltay=maximum(yvec)-minimum(yvec)
    deltaz=maximum(zvec)-minimum(zvec)
    mindelta=min(deltax,deltay,deltaz)
    maxdelta=max(deltax,deltay,deltaz)
    if mindelta<maxdelta*0.001
        eps_delta=maxdelta*0.001
    else
        eps_delta=0
    end 
    ax=(deltax+eps_delta)/(mindelta+eps_delta)
    ay=(deltay+eps_delta)/(mindelta+eps_delta)
    az=(deltaz+eps_delta)/(mindelta+eps_delta)

    resolution_val = 800  #1200



    ax1=current_axis()
    if i_var==0
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(0,1))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.5),colorrange=(0,1));  
    elseif i_var==1
        maxval=maximum(cgammavec[:])
        minval=minimum(cgammavec[:])
        #maxval=0.0055
        #minval=0.003
        @info "minval = $minval"
        @info "maxval = $maxval"
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(minval,maxval))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]       
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.75),colorrange=(minval,maxval));  
    elseif i_var==2
        maxval=maximum(cgammavec[:])
        minval=minimum(cgammavec[:])
        #deltaval=max(maxval-minval,0.1*(abs(maxval)))
        #maxval=maxval+deltaval
        #minval=minval-deltaval
        @info "minval = $minval"
        @info "maxval = $maxval"
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(0,1))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]       
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.75),colorrange=(0,1));  
    elseif i_var==3
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(0,1))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.5),colorrange=(0,1));  
    elseif i_var==4
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(0,1))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.5),colorrange=(0,1));  
    elseif i_var==5
        maxval=maximum(cgammavec[:])
        minval=minimum(cgammavec[:])
        @info "minval = $minval"
        @info "maxval = $maxval"
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(minval,maxval))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]       
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.75),colorrange=(minval,maxval));  
    elseif i_var==6
        maxval=maximum(cgammavec[:])
        minval=minimum(cgammavec[:])
        @info "minval = $minval"
        @info "maxval = $maxval"
        empty!(ax1)
        p1 = poly!(connect(xyz, GeometryBasics.Point{3}), connect(1:3*N, TriangleFace); color=cgammavec[:], strokewidth=1,colorrange=(minval,maxval))
        fig1=current_figure(); [delete!(col1) for col1 in fig1.content if col1 isa Colorbar]       
        col1=Colorbar(fig1, colormap = :viridis,  vertical=true, height=Relative(0.75),colorrange=(minval,maxval));  
    end
end