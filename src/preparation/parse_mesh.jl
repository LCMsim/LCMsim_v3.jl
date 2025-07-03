using GeometryBasics

function __parse_gmsh(filename::String)
    ind = Int32(1)
    gridind = Int32(1)
    orig_sets = []
    origgridid = []
    gridx::Vector{Float64} = []
    gridy::Vector{Float64} = []
    gridz::Vector{Float64} = []
    celloriggridid = []
    pid::Vector{Int8} = []
    cellgridid = Matrix{Int32}(undef, 0, 3)

    open(filename, "r") do fid
        line = 1
        while !eof(fid)
            thisline = readline(fid)
            if length(thisline) >= 8
                card = thisline[1:8]
                if cmp(card, "GRID    ") == 0
                    gridindstring = thisline[9:16]
                    origgridid = vcat(origgridid, parse(Int32, gridindstring))
                    txt = thisline[25:32]
                    txt = replace(txt, " " => "")
                    txt = replace(txt, "E" => "")
                    txt = replace(txt, "e" => "")
                    txt1 = replace(txt, "-" => "e-")
                    txt1 = replace(txt1, "+" => "e+")
                    if cmp(txt1[1], 'e') == 0
                        txt2 = txt1[2:end]
                    else
                        txt2 = txt1
                    end
                    val = parse(Float64, txt2)
                    val1 = val
                    txt = thisline[33:40]
                    txt = replace(txt, " " => "")
                    txt = replace(txt, "E" => "")
                    txt = replace(txt, "e" => "")
                    txt1 = replace(txt, "-" => "e-")
                    txt1 = replace(txt1, "+" => "e+")
                    if cmp(txt1[1], 'e') == 0
                        txt2 = txt1[2:end]
                    else
                        txt2 = txt1
                    end
                    val = parse(Float64, txt2)
                    val2 = val
                    txt = thisline[41:48]
                    txt = replace(txt, " " => "")
                    txt = replace(txt, "E" => "")
                    txt = replace(txt, "e" => "")
                    txt1 = replace(txt, "-" => "e-")
                    txt1 = replace(txt1, "+" => "e+")
                    if cmp(txt1[1], 'e') == 0
                        txt2 = txt1[2:end]
                    else
                        txt2 = txt1
                    end
                    val = parse(Float64, txt2)
                    val3 = val
                    gridx = vcat(gridx, Float64(val1))
                    gridy = vcat(gridy, Float64(val2))
                    gridz = vcat(gridz, Float64(val3))
                    gridind = gridind + 1
                elseif cmp(card, "CTRIA3  ") == 0
                    celloriggridid = vcat(celloriggridid, parse(Int32, thisline[9:16]))
                    pid = vcat(pid, parse(Int8, thisline[17:24]))
                    i1val = parse(Int32, thisline[25:32])
                    i1 = findfirst(isequal(i1val), origgridid)
                    i2val = parse(Int32, thisline[33:40])
                    i2 = findfirst(isequal(i2val), origgridid)
                    i3val = parse(Int32, thisline[41:end])
                    i3 = findfirst(isequal(i3val), origgridid)
                    ivec = [i1, i2, i3]
                    idel = findall(isequal(min(ivec[1], ivec[2], ivec[3])), ivec)
                    deleteat!(ivec, idel)
                    idel = findall(isequal(max(ivec[1], ivec[2])), ivec)
                    deleteat!(ivec, idel)
                    cellgridid = vcat(cellgridid, [min(i1, i2, i3) ivec[1] max(i1, i2, i3)])
                    ind = ind + 1
                end
            end
            line += 1
        end
    end
    N = ind - 1  #total number of cells

    # convert indices
    sets = []
    for set in (orig_sets)
        new_set = []
        for i in 1:length(set)
            i1 = findfirst(isequal(orig_sets[i]), celloriggridid)
            new_set = vcat(new_set, i1)
        end
        push!(sets, new_set)
    end

    grid = [gridx gridy gridz]

    cellgridid = Int32.(cellgridid)
    N = Int32(N)

    vertices = [Point3{Float64}(grid[i, :]) for i in 1:size(grid)[1]]

    return Int.(cellgridid), vertices::Vector{Point3{Float64}}, Int.(pid), Int(N)
end

function __parse_HyperMesh(filename::String)

    # TODO clean this up, do logic here
    grid_vec, ctria_vec = parse_HyperMeshNastran(filename)

    pids = Vector{Int}(undef, 0)
    cellgridid = Matrix{Int}(undef, 0, 3)
    for i in 1:size(ctria_vec)[1]
        push!(pids, ctria_vec[i][1])
        cellgridid = vcat(cellgridid, transpose(sort(ctria_vec[i][2:4])))
    end
    grid = Vector{Point3{Float64}}(undef, 0)
    for i in 1:size(grid_vec)[1]
        push!(grid, Point3{Float64}(grid_vec[i]))
    end

    return cellgridid, grid, pids, size(cellgridid)[1]
end

function __parse_abaqus(filename::String)
    grid = Dict()
    ctria = Dict()
    sets = []
    temp_set = []
    set_parsing_active = false
	element_parsing_active = false
	node_parsing_active = false
    dummy_part_id = 1
	ind=Int(1);
	gridind = Int(1);
	nodeind = Int(1);
	elementind = Int(1);
	setelementind = Int(1);
	setind = Int(0);
	pids = []
    origgridid=[];
    gridx=[];
    gridy=[];
    gridz=[];
    celloriggridid=[];
    cellgridid=Array{Int64}(undef, 0, 3);

    #COb: Read not used sets from file "_filename.inp"
    notusedsets = Vector{Int}(undef, 0)
    inputfile_parts=splitpath(filename)
    inputfile_parts[end]="_" * inputfile_parts[end]
    _inputfile=joinpath(inputfile_parts)
    if isfile(_inputfile) 
        lines = readlines(_inputfile)
        ids = [parse(Int, id) for id in split(lines[1],",")]
        for id in ids
            push!(notusedsets, id)
        end
    end

    lines = readlines(filename)

    #COb: Modified to assign set numbers in increasing order of reading from mesh file, starting from 2
    i_set=2
    for line in lines
	    thisline=replace(line," "=> "");
		len_line=length(thisline)

		if node_parsing_active
			if isempty(thisline)==1 || cmp(thisline[1:1],"*")==0
				node_parsing_active=false
			else
				txt_vec=split(thisline,",")
				gridindstring=txt_vec[1]
				origgridid=vcat(origgridid,parse(Int64,gridindstring));                        
				txt=txt_vec[2]
				txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "");
				txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+");
				if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end;
				val=parse(Float64,txt2);
				val1=val;
				txt=txt_vec[3]
				txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "");
				txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+");
				if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end;
				val=parse(Float64,txt2);
				val2=val;
				txt=txt_vec[4]
				txt=replace(txt," "=> "");txt=replace(txt,"E" => "");txt=replace(txt,"e" => "");
				txt1=replace(txt,"-" => "e-");txt1=replace(txt1,"+" => "e+");
				if cmp(txt1[1],'e')==0;txt2=txt1[2:end];else;txt2=txt1;end;
				val=parse(Float64,txt2);
				val3=val;
				gridx=vcat(gridx,Float64(val1));
				gridy=vcat(gridy,Float64(val2));
				gridz=vcat(gridz,Float64(val3));
				gridind=gridind+1;				
	
                id=parse(Int64,gridindstring)
				pos = length(grid) + 1
				grid[id] = (pos, [Float64(val1), Float64(val2), Float64(val3)])
			end
		end
				
		if element_parsing_active
			if isempty(thisline) == 1 || cmp(thisline[1:1],"*") == 0
				element_parsing_active=false
			else
				txt_vec=split(thisline,",")
				celloriggridid=vcat(celloriggridid,parse(Int64,txt_vec[1]));
				i1val=parse(Int64,txt_vec[2]);
				i1=findfirst(isequal(i1val),origgridid);
				i2val=parse(Int64,txt_vec[3]);
				i2=findfirst(isequal(i2val),origgridid);
				i3val=parse(Int64,txt_vec[4]);
				i3=findfirst(isequal(i3val),origgridid);
				ivec=[i1,i2,i3];
				idel=findall(isequal(min(ivec[1],ivec[2],ivec[3])),ivec);
				deleteat!(ivec,idel)
				idel=findall(isequal(max(ivec[1],ivec[2])),ivec);
				deleteat!(ivec,idel)                        
				cellgridid=vcat(cellgridid,[min(i1,i2,i3) ivec[1] max(i1,i2,i3)]);
				ind=ind+1;    

                id=parse(Int64,txt_vec[1])
				verts=[min(i1,i2,i3) ivec[1] max(i1,i2,i3)]			
			    ctria[id] = (dummy_part_id, verts)	
			end
		end
		
		if set_parsing_active
			if isempty(thisline)==1 || cmp(thisline[1:1],"*")==0
				part_id = i_set
				if (i_set in notusedsets)
                    part_id = 1
                end
                i_set=i_set+1
				
				push!(sets, (part_id, pids))	

				set_parsing_active=false
			elseif isnothing(tryparse(Int64,thisline[1:1]))    #name instead of numbers
				set_parsing_active=false
				setind=setind-1
			else
				txt1=thisline
				txt1=replace(txt1," "=> "")
				txt2=split(txt1,",")
				for i in 1:length(txt2)
					if !isempty(txt2[i])						
                        append!(pids, parse(Int64,txt2[i]))  #push!(pids, parse(Int64,txt2[i]))
					end
				end			
			end
		end		
		
		
		if len_line>=5 &&  (cmp( thisline[1:5],"*Node")==0 || cmp( thisline[1:5],"*NODE")==0 ) #if the first five characters are *Node: isnodedefinition=1; else: isnodedefinition=0
			node_parsing_active=true  #isnodedefinition=Int64(1);
			nodeind=1
			gridind=1
		end
		if (len_line>=17 &&  (cmp(thisline[1:17],"*Element, TYPE=S3")==0 || cmp(thisline[1:17],"*Element, Type=S3")==0 || cmp(thisline[1:17],"*ELEMENT, TYPE=S3")==0 || cmp(thisline[1:17],"*ELEMENT, Type=S3")==0) ) ||
			len_line>=16 &&  (cmp(thisline[1:16],"*Element,TYPE=S3")==0 || cmp(thisline[1:16],"*Element,Type=S3")==0 || cmp(thisline[1:16],"*ELEMENT,TYPE=S3")==0 || cmp(thisline[1:16],"*ELEMENT,Type=S3")==0)  #if the first 17 characters are *Element, TYPE=S3                    
			element_parsing_active=true  #iselementdefinition=Int64(1);
			ind=1
			elementind=1
		end
		if len_line>=6 &&  (cmp( thisline[1:6],"*Elset")==0 || cmp( thisline[1:6],"*ELSET")==0)  #if the first six characters are *ELSET
			set_parsing_active=true  #issetdefinition=Int64(1);
			setelementind=1
			setind=setind+1
            pids=[]  #empty!(pids)
		end
		
		
		
    end

    #COb: if _pset.csv is available in directory with meshfile
    #     change part_ID to an inlet for all cells which lie inside a sphere with radius
    filename_parts=splitpath(filename)
    _psetfile=joinpath( joinpath(filename_parts[1:end-1]) ,"_pset.csv")
    if isfile(_psetfile)
        #read _psetfile: first line radius, following line inlet points
        lines = readlines(_psetfile);n_lines=length(lines);r_p=parse(Float64,lines[1])
        inletpos_xyz=Array{Float64}[]
        for i_line in 2:n_lines; 
            str=split(lines[i_line],",")
            xpos=parse(Float64,str[1])
            ypos=parse(Float64,str[2])
            zpos=parse(Float64,str[3])
            inletpos_xyz=push!(inletpos_xyz,[xpos ypos zpos])
        end
        part_id = 6
        ids = []
        for i_inlet in 1:length(inletpos_xyz)
            #loop over all cells, check if cell center lies inside a sphere with radius r. if so, add cell to set with part_id=6
            xyz_inlet=inletpos_xyz[i_inlet,:]
            for (i, cell) in ctria
                verts = cell[2]
                xvec=0.
                yvec=0.
                zvec=0.
                for i_vert in 1:3
                    xvec=xvec+0.3333*grid[verts[i_vert]][2][1]
                    yvec=yvec+0.3333*grid[verts[i_vert]][2][2]
                    zvec=zvec+0.3333*grid[verts[i_vert]][2][3]
                end
                vec1=[inletpos_xyz[i_inlet][1]-xvec, inletpos_xyz[i_inlet][2]-yvec, inletpos_xyz[i_inlet][3]-zvec ]
                if sqrt(dot(vec1,vec1))<=r_p;     
                    #@info "xvec,yvec,zvec = $xvec,$yvec,$zvec"          
                    id = i
                    #@info "cid $id in inlet $i_inlet"
                    append!(ids, id)
                end
            end     
        end
        ids=unique(ids)
        push!(sets, (part_id, ids))
    end

    ctria_vec = Vector{Any}(nothing, length(ctria))
    grid_vec = Vector{Any}(nothing, length(grid))
    for (part_id, cells) in sets
        for cid in cells
            ctria[cid] = (part_id, ctria[cid][2])
        end
    end

    # rename verts
    for (i, cell) in enumerate(ctria)
        cid = cell[1]
        pid = cell[2][1]
        verts = cell[2][2]
        verts = [grid[gid][1] for gid in verts]
        ctria_vec[i] = [pid, verts...]   #COb: index bug fix
    end

    for vert in grid
        gid = vert[2][1] # gid after renaming
        grid_vec[gid] = vert[2][2]
    end

    pids = Vector{Int}(undef, 0)
    cellgridid = Matrix{Int}(undef, 0, 3)
    for i in 1:size(ctria_vec)[1]
        push!(pids, ctria_vec[i][1])
        cellgridid = vcat(cellgridid, transpose(sort(ctria_vec[i][2:4])))
    end
    grid = Vector{Point3{Float64}}(undef, 0)
    for i in 1:size(grid_vec)[1]
        push!(grid, Point3{Float64}(grid_vec[i]))
    end

    #@info "cellgridid = " * string(cellgridid)
    #@info "grid = " * string(grid)
    #@info "pids = " * string(pids)
    #@info "N = " * string(size(cellgridid)[1])    
    #error("22")

    return cellgridid, grid, pids, size(cellgridid)[1]
end

function parse_mesh(meshfile::String)

    file_extension = splitext(meshfile)[2]

    cellgridid = Matrix{Int}(undef, 0, 3)
    vertices = Vector{Point3{Float64}}(undef, 0)
    pids = Vector{Int}(undef, 0)
    N = 0

    if file_extension == ".bdf"
        @debug "PARSE GMSH"
        cellgridid, vertices, pids, N = __parse_gmsh(meshfile)
    elseif file_extension == ".dat"
        @debug "PARSE HyperMesh"
        cellgridid, vertices, pids, N = __parse_HyperMesh(meshfile)
    elseif file_extension == ".inp"
        @debug "PARSE Abaqus"
        cellgridid, vertices, pids, N = __parse_abaqus(meshfile)
    else
        @error "File extension not supported."
    end

    return cellgridid, vertices, pids, N
end