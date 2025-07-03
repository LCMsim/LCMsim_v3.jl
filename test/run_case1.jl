include("../src/LCMsim_v2.jl")

mypath=pwd()
savepath = joinpath(mypath,"test")

i_case=23  #31
if i_case==21  
    #VARI case 1: Validation case without DM from Q. Govignon, S. Bickerton, P.A. Kelly, 
    #             Simulation of the reinforcement compaction and resin flow during the complete resin infusion process    
    meshfile = joinpath(mypath,"test","mesh_21.dat")
    partfile = joinpath(mypath,"test","part_description_21.csv")
    simfile = joinpath(mypath,"test","simulation_params_21.csv")
    i_model=3
    t_max = 200.
elseif i_case==22  
    #VARI case 2: Validation case with DM from from J. Sirtautas, A.K. Pickett, A. George, 
    #             Materialscharacterisation and analysis for flow simulation of liquid resin infusion   
    meshfile = joinpath(mypath,"test","mesh_22.dat")
    partfile = joinpath(mypath,"test","part_description_22.csv")
    simfile = joinpath(mypath,"test","simulation_params_22.csv")
    i_model=3
    t_max = 23.  #23.
elseif i_case==23
    #VARI case 3: Empirical observation from H. M. Andersson, T. S. Lundström, B. R. Gebart, 
    #             Numerical model for vacuum infusion manufacturing of polymer composites
    meshfile = joinpath(mypath,"test","mesh_23.dat")
    partfile = joinpath(mypath,"test","part_description_23.csv")
    simfile = joinpath(mypath,"test","simulation_params_23.csv")
    i_model=3
    t_max = 23.  #23
elseif i_case==0
    #Radial test case 0 with multiple patches
    meshfile = joinpath(mypath,"test","mesh_0.dat")
    partfile = joinpath(mypath,"test","part_description_0.csv")
    simfile = joinpath(mypath,"test","simulation_params_0.csv")
    i_model=2  #3  #  
    t_max = 300.
elseif i_case==1  
    #Radial test case 1 from https://obertscheiderfhwn.github.io/RTMsim/build/tutorials/
    meshfile = joinpath(mypath,"test","mesh_1.dat")
    partfile = joinpath(mypath,"test","part_description_1.csv")
    simfile = joinpath(mypath,"test","simulation_params_1.csv")
    i_model=2  #3  #  
    t_max = 200.
elseif i_case==2
    #Radial test case 2 from https://obertscheiderfhwn.github.io/RTMsim/build/tutorials/
    meshfile = joinpath(mypath,"test","mesh_2.dat")
    partfile = joinpath(mypath,"test","part_description_2.csv")
    simfile = joinpath(mypath,"test","simulation_params_2.csv")
    i_model=2  #3  #  
    t_max = 200.
elseif i_case==31  #FPCM16 VARI test case for verification with OpenFOAM
    meshfile = joinpath(mypath,"test","mesh_31.dat")
    partfile = joinpath(mypath,"test","part_description_31.csv")
    simfile = joinpath(mypath,"test","simulation_params_31.csv")
    i_model=3
    t_max = 30.
elseif i_case==33  #Radial flow similar to 31 but coarse mesh for quick tests
    meshfile = joinpath(mypath,"test","mesh_33.dat")
    partfile = joinpath(mypath,"test","part_description_33.csv")
    simfile = joinpath(mypath,"test","simulation_params_33.csv")
    i_model=3
    t_max = 160.
elseif i_case==7  #AMPCS case 1 from https://www.tandfonline.com/doi/full/10.1080/20550340.2023.2282310, Validation was performed with exp_val=4, now exp_val=25 is default, resulting is small discrepacies
    meshfile = joinpath(mypath,"test","mesh_7.dat")
    partfile = joinpath(mypath,"test","part_description_7.csv")
    simfile = joinpath(mypath,"test","simulation_params_7.csv")
    i_model=2
    t_max = 120.
elseif i_case==8  #AMPCS case 2 from https://www.tandfonline.com/doi/full/10.1080/20550340.2023.2282310, Validation was performed with exp_val=4, now exp_val=25 is default, resulting is small discrepacies
    meshfile = joinpath(mypath,"test","mesh_8.dat")
    partfile = joinpath(mypath,"test","part_description_8.csv")
    simfile = joinpath(mypath,"test","simulation_params_8.csv")
    i_model=2
    t_max = 80.
elseif i_case==5  #Annulus filler test case from https://obertscheiderfhwn.github.io/RTMsim/build/tutorials/
    meshfile = joinpath(mypath,"test","mesh_5.dat")
    partfile = joinpath(mypath,"test","part_description_5.csv")
    simfile = joinpath(mypath,"test","simulation_params_5.csv")
    i_model=2
    t_max = 250.
end
i_var=4     #0..filling,1..thickness,2..porosity,3..filling DM,4..filling combined DM and FP,5..w

@info "mypath = $mypath"
@info "meshfile = $meshfile"
@info "partfile = $partfile"
@info "simfile = $simfile"

t_step = t_max/16
if i_model == 1
    modeltype = LCMsim_v2.model_1
elseif i_model == 2         
    modeltype = LCMsim_v2.model_2
else
    modeltype = LCMsim_v2.model_3
end  

filename_parts=splitpath(meshfile)
meshfilename_parts=splitpath(meshfile)
meshfilename_parts[end]="_" * meshfilename_parts[end]
writefilename=joinpath(meshfilename_parts)

i_advanced=2

LCMsim_v2.create_and_solve(savepath,meshfile,partfile,simfile,modeltype,i_advanced,t_max,t_step,LCMsim_v2.verbose,true,true)

include("plot_case1.jl");

