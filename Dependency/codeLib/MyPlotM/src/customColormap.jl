using PyCall,EzXML
export make_cmap,get_cmap


function get_cmap(cmName)
    path_currModule="$(dirname(pathof(@__MODULE__)))/.."
    if lowercase(cmName) == "sand"
        file="$path_currModule/xml_cmap/08brown-brye1b.xml"
        cmap_cur=make_cmap(file)
        cmap_cur.set_over("k")
        cmap_cur.set_under("w")
    end
    return cmap_cur
end




"""

sample usage:
file="MyPlotM/xml_cmap/08brown-brye1b.xml"
cMap=make_cmap(file)

"""

function make_cmap(xml)
    colorMap=firstelement(root(readxml(xml)))
    pos=[parse(Float32,x["x"]) for x in findall(".//Point", colorMap)]
    r=[parse(Float32,x["r"]) for x in findall(".//Point", colorMap)]
    g=[parse(Float32,x["g"]) for x in findall(".//Point", colorMap)]
    b=[parse(Float32,x["b"]) for x in findall(".//Point", colorMap)]
    if length(pos) != length(r)
        error("pos length must be the same as colors")
    elseif (pos[1] != 0) || (pos[end] != 1)
        error("pos must start with 0 and end with 1")
    end
    cdict = Dict("red" => [], "green" => [], "blue" => [])
    if pos[1] != 0
        push!(cdict["red"], (0, r[1], r[1]))
        push!(cdict["green"], (0, g[1], g[1]))
        push!(cdict["blue"], (0, b[1], b[1]))
    end
    for (n_pos,pos_curr) in enumerate(pos)
        push!(cdict["red"], (pos_curr, r[n_pos], r[n_pos]))
        push!(cdict["green"], (pos_curr, g[n_pos], g[n_pos]))
        push!(cdict["blue"], (pos_curr, b[n_pos], b[n_pos]))
    end
    if pos[end] != 1
        push!(cdict["red"], (1, r[end], r[end]))
        push!(cdict["green"], (1, g[end], g[end]))
        push!(cdict["blue"], (1, b[end], b[end]))
    end
    mpl=pyimport("matplotlib")
    cmap = mpl.colors.LinearSegmentedColormap("my_colormap", cdict, 256)
    return cmap
end


