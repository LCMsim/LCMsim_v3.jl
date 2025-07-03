using Printf

function tabToSpace(s::String, n=8)
    i = findfirst("\t", s)
    while !isnothing(i)
        i = i[1]
        m = i % n
        if m == 1
            tab = repeat(" ", n)
        elseif m == 0
            tab = " "
        else
            tab = repeat(" ", n - (i % n) + 1)
        end
        s = s[1:i - 1] * tab * s[i + 1:end] 
        i = findfirst("\t", s)
    end
    return s
end

function parse_HyperMeshNastran(inputfile::String)
    grid = Dict()
    ctria = Dict()
    sets = []
    temp_set = []
    set_parsing_active = false
    dummy_part_id = 1

    #COb: Read not used sets from file "_inputfilename.dat"
    notusedsets = Vector{Int}(undef, 0)
    inputfile_parts=splitpath(inputfile)
    inputfile_parts[end]="_" * inputfile_parts[end]
    _inputfile=joinpath(inputfile_parts)
    if isfile(_inputfile) 
        lines = readlines(_inputfile)
        ids = [parse(Int, id) for id in split(lines[1],",")]
        for id in ids
            push!(notusedsets, id)
        end
    end

    lines = readlines(inputfile)

    #COb: Modified to assign set numbers in increasing order of reading from mesh file, starting from 2
    i_set=2
    for line in lines
        if set_parsing_active
            if isdigit(last(strip(line)))
                push!(temp_set, line)
                part_id = i_set  #parse(Int, match(Regex("[0-9]+\\s="), temp_set[1]).match[1:end-1]) + 1  
                if (i_set in notusedsets)
                    part_id = 1
                end
                i_set=i_set+1
                ids = []
           
                for l in temp_set
                    ids_string = match(Regex("([0-9]+,)+[0-9]+"), l).match
                    append!(ids, split(ids_string, ","))
                end
           
                ids = [parse(Int, id) for id in ids]
           
                push!(sets, (part_id, ids))
                temp_set = []
                set_parsing_active = false
            else
                push!(temp_set, line)
            end
        elseif !isnothing(match(Regex("GRID\\s*[0-9]+"), line))
            id = parse(Int, match(Regex("\\s[0-9]+\\s"), line).match)
            coords = line[25:end]

            (x, y, z) = parse_coordinate_line(String(coords))

            pos = length(grid) + 1
            grid[id] = (pos, [x, y, z])

        elseif !isnothing(match(Regex("CTRIA3\\s*[0-9]+"), line))
            id = parse(Int, match(Regex("\\s[0-9]+\\s"), line).match)
            verts = match(Regex("(\\s+[0-9]+){3}\$"), line).match
            verts = split(verts, Regex("\\s+"))
            verts = [parse(Int, x) for x in verts[2:end]]

            ctria[id] = (dummy_part_id, verts)

        elseif !isnothing(match(Regex("SET\\s[0-9]+\\s="), line))
            set_parsing_active = true
            push!(temp_set, line)

            if isdigit(last(strip(line)))
                part_id = i_set  #parse(Int, match(Regex("[0-9]+\\s="), temp_set[1]).match[1:end-1]) + 1 
                if (i_set in notusedsets)
                    part_id = 1
                end
                i_set=i_set+1
                ids = []
            
                for l in temp_set
                    ids_string = match(Regex("([0-9]+,)+[0-9]+"), l).match
                    append!(ids, split(ids_string, ","))
                end
            
                ids = [parse(Int, id) for id in ids]
            
                push!(sets, (part_id, ids))
                temp_set = []
                set_parsing_active = false
            end
        end
    end

    #COb: if _pset.csv is available in directory with meshfile
    #     change part_ID to an inlet for all cells which lie inside a sphere with radius
    filename_parts=splitpath(inputfile)
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
        ctria_vec[i] = [pid, verts...]  #ctria_vec[cid] = [pid, verts...]  #COb: index bug fix
    end


    for vert in grid
        gid = vert[2][1] # gid after renaming
        grid_vec[gid] = vert[2][2]
    end

    return grid_vec, ctria_vec
end

function write_GmshNastran(outputfile::String, grid, ctria)

    open(outputfile, "w") do f

        for (gid, vert) in enumerate(grid)
            x = ""
            y = ""
            if vert[1] == 0.
                x = "0.00E+00"
            else
                if vert[1] < 0.
                    x = @sprintf("%1.5f", vert[1])
                else
                    x = @sprintf("%1.6f", vert[1])
                end
            end
            
            if vert[2] == 0.
                y = "0.00E+00"
            else
                if vert[2] < 0.
                    y = @sprintf("%1.5f", vert[2])
                else
                    y = @sprintf("%1.6f", vert[2])
                end
            end

            l = "GRID    " * string(gid) * "\t0\t" * x * y * "0.00E+00" * "\n"
            write(f, tabToSpace(l))
        end

        for (cid, cell) in enumerate(ctria)
            l = "CTRIA3\t" * string(cid) * "\t" * string(cell[1]) * "\t" * string(cell[2]) * "\t" * string(cell[3]) * "\t" * string(cell[4]) * "\n"
            write(f, tabToSpace(l))
        end

        write(f, "ENDDATA\n")
    end

end


function parse_coordinate_line(line::String)::Tuple{Float64, Float64, Float64}
    txt = line[1:8]
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
    x = parse(Float64, txt2)

    txt = line[9:16]
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
    y = parse(Float64, txt2)

    txt = line[17:end]
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
    z = parse(Float64, txt2)

    return x, y, z
end