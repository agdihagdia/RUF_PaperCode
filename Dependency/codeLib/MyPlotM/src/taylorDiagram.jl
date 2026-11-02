

import PyPlot as plt
using PyCall
export taylorDiagram,getTaylorXY

function getTaylorXY(corr,stdd)
    xx,yy=acos.(corr),stdd
    return xx,yy
end

function taylorDiagram(;fig = plt.gcf(),refstd=1,rlocs = vcat([0.1:0.1:0.9;],0.95,0.99),ticks_std=nothing,
    pos=plt.gca().get_position(), level_rmsd=ticks_std,addLabel_std=true,add_clabel=true)
    if isnothing(ticks_std); ticks_std=[0:0.25:1.5;].*refstd;end
    if isnothing(level_rmsd); level_rmsd=ticks_std;end
    
    FA=pyimport("mpl_toolkits.axisartist.floating_axes")
    GF=pyimport("mpl_toolkits.axisartist.grid_finder")
    proj=pyimport("matplotlib.projections")
    tr = proj.PolarAxes.PolarTransform()
    gl2 = GF.FixedLocator(ticks_std)
    # Correlation labels
    # rlocs = [0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1]
    
    tmax = pi/2

    tlocs = acos.(rlocs)        # Conversion to polar angles
    gl1 = GF.FixedLocator(tlocs)    # Positions
    tf1 = GF.DictFormatter(Dict(tlocs[i] => string(rlocs[i]) for i in eachindex(tlocs)))

    
    # Standard deviation axis extent (in units of reference stddev)
    smin,smax=0.0,maximum(ticks_std)

    ghelper = FA.GridHelperCurveLinear(
        tr,
        extremes=(0, tmax, smin, smax),
        grid_locator1=gl1, tick_formatter1=tf1,
        grid_locator2=gl2)
     
    ax = FA.FloatingSubplot(fig, 111, grid_helper=ghelper)
    fig.add_subplot(ax)

    # Adjust axes
    ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
    ax.axis["top"].toggle(ticklabels=true, label=true)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("R")
    # ax.axis["top"].label.set_text("Correlation Coefficient")

    ax.axis["left"].set_axis_direction("bottom")  # "X axis"
    if addLabel_std ;ax.axis["left"].label.set_text("Standard deviation");end

# ax.text(-0.2,0.5,  "Standard deviation",
# horizontalalignment = "center", verticalalignment = "center",
# rotation =  90,
# fontsize = 10)
# ax.axis["right"].set_label("slakdjf")

    # ax.axis["right"].label.set_text("Standard deviation")

    ax.axis["right"].set_axis_direction("top")    # "Y-axis"
    ax.axis["right"].toggle(ticklabels=true)
    ax.axis["right"].major_ticklabels.set_axis_direction("left")

    ax.axis["bottom"].toggle(ticklabels=false, label=false)

    # ax.set_position([0.1, 0.1, 0.9,0.8])
    ax.set_position(pos)
    _ax = ax                   # Graphical axes
    ax = ax.get_aux_axes(tr)   # Polar coordinates
    # Add reference point and stddev contour
    h_obs = ax.plot([0], refstd, "k*",
                    ls="", ms=10, label="Observation")
    _ax.grid()

    # Add Taylor contours centered on the reference point
    rs= collect(range(smin, smax,100)) * ones(100)'
    ts= ones(100) * collect(range(0, tmax,100))'
    rms = sqrt.(refstd^2 .+ rs.^2 .- 2 .* refstd .* rs .* cos.(ts))
    h_contours = ax.contour(ts, rs, rms, level_rmsd,linestyles="dashed",colors="#00b294")
    if add_clabel; plt.clabel(h_contours, inline=true, fontsize=8); end
    return _ax,ax,h_obs,h_contours

end