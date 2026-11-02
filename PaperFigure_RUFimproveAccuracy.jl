using Revise,MyPlotM,StaDataTool,NaqpmsUtilities,utilities,Dates,NaNStatistics
using GradDataM,TimeZones,Statistics,ZZDBM,PlotFunM
using CSV,DataFrames,Base.Threads
using Base.Iterators,PyCall
using Distributed,Base.Filesystem
include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle
matplotlib = pyimport("matplotlib")
configure_matplotlib!(matplotlib)
ccrs = pyimport("cartopy.crs")
ticker=matplotlib.ticker
shpreader = pyimport("cartopy.io.shapereader")
np = pyimport("numpy")
mcolors = pyimport("matplotlib.colors")
fm=matplotlib.font_manager

function get_staIDs(city)
    staIDs=staInfo.staID[in.(staInfo.city,[[city]])]
    return staIDs
end
#draw
function addGrids(ax,rightlabel)
    gl = ax.gridlines(crs=PlaCar,draw_labels=true, dms=true, x_inline=false, y_inline=false, linewidth=0.5, color=[0.1,0.1,0.1], alpha=0.5,zorder=2,linestyle="--") # 添加经纬度网格，这里设定为黑色虚线半透明 
    gl.top_labels  = false;    gl.right_labels = rightlabel;   gl.rotate_labels = false ; gl.left_labels = false
    gl.xlocator=ticker.FixedLocator([40,60,80,100,120,140,160,180])
    gl.xlabel_style = Dict(:size => TICK_SIZE, :fontname => "Arial")
    gl.ylabel_style = Dict(:size => TICK_SIZE, :fontname => "Arial")
end 
function add_SouthChinaAxes(axIn,sta_lon,sta_lat,rmse,color_map,color_norm)
    fig=axIn.get_figure()
    posIn=axIn.get_position()
    in_xmax,in_ymin,in_width,in_height=posIn.xmax,posIn.ymin,posIn.width,posIn.height
    hight2width_fig=fig.figbbox.height/fig.figbbox.width

    width=in_width/6
    height=width*2/hight2width_fig
    position_southChina=[in_xmax-width,in_ymin-height/8,width,height]

    # println(position_southChina)
    ax_southChina=fig.add_axes([0, 0, 0.1, 0.1], projection=axIn.projection)
    ax_southChina.set_extent([106, 120, 2, 24],crs=ccrs.PlateCarree())
    ax_southChina.set_position(position_southChina)
    cn_boundary=ax_southChina.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    ax_southChina.scatter(sta_lon,sta_lat,c=rmse,s=0.2,transform=PlaCar,
        cmap=color_map,norm=color_norm)
    return ax_southChina
end
function add_SouthChinaAxes(axIn,sta_lon,sta_lat,rmse)
    fig=axIn.get_figure()
    posIn=axIn.get_position()
    in_xmax,in_ymin,in_width,in_height=posIn.xmax,posIn.ymin,posIn.width,posIn.height
    hight2width_fig=fig.figbbox.height/fig.figbbox.width

    width=in_width/6
    height=width*2/hight2width_fig
    position_southChina=[in_xmax-width,in_ymin-height/8,width,height]

    # println(position_southChina)
    ax_southChina=fig.add_axes([0, 0, 0.1, 0.1], projection=axIn.projection)
    ax_southChina.set_extent([106, 120, 2, 24],crs=ccrs.PlateCarree())
    ax_southChina.set_position(position_southChina)
    cn_boundary=ax_southChina.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    ax_southChina.scatter(sta_lon,sta_lat,c=rmse,s=0.2,transform=PlaCar,cmap=cmapM)
    return ax_southChina
end
function add_label_atleftup(ax, label,c="black"; xpos=0.02, ypos=0.95, fontsize=PANEL_TITLE_SIZE)
    ax.text(xpos, ypos, label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            va="top",
            ha="left",
            fontproperties=font_propbd,
            color=c)
end
function get_truncated_cmap(base_cmap_name, min_val=0.0, max_val=1.0, n=256)
    base_cmap = plt.get_cmap(base_cmap_name)
    # 在 min_val 到 max_val 之间重采样
    colors_resampled = base_cmap(np.linspace(min_val, max_val, n))
    # 创建新的颜色条对象
    new_cmap = mcolors.LinearSegmentedColormap.from_list("trunc_afmhot", colors_resampled)
    return new_cmap
end

if false
staObs=staDATA("/sh2/home/3clear/wuhj/paper_IETMsand/ori202504.mat")
staInfo=staObs.staInfo
cities=unique(staObs.staInfo.city)
var2plot="PM10"
PlaCar=ccrs.PlateCarree()

matFun3d_v5(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_AllStaOri/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"

leadtimes=collect(1:72)
tStart_all_obs=ZonedDateTime.(DateTime(2025,4,11):Hour(1):DateTime(2025,4,21),tz"Asia/Shanghai")
tStart_all_fore=ZonedDateTime.(DateTime(2025,4,7,20):Hour(1):DateTime(2025,4,21,19),tz"Asia/Shanghai")
tStart_One = ZonedDateTime.(DateTime(2025,4,11),tz"Asia/Shanghai")

filepath_IETM=matFun3d_v5.(tStart_all_fore)
IETM3D_data=staDATA(filepath_IETM,Times=tStart_all_fore)

# obs=Array{Array{Float64},2}(undef,length(cities),length(leadtimes))
# IETM=Array{Array{Float64},2}(undef,length(cities),length(leadtimes))
# IETM_24h=Array{Array{Float64},2}(undef,length(cities),length(leadtimes))
# for n_c in eachindex(cities)
#     city=cities[n_c]
#     for n_l in eachindex(leadtimes)
#         l=leadtimes[n_l] 
#         obs[n_c,n_l]=vec(nanmean(read(staObs,Times=tStart_all_obs,varNames="PM10",staIDs=get_staIDs(city)),dims=3))
#         IETM[n_c,n_l]=vec(nanmean(read(IETM3D_data,Times=tStart_all_obs .- Hour(l),leadTimes=l,varNames="PM10",staIDs=get_staIDs(city)),dims=3))
#     end
# end
# for n_c in eachindex(cities)
#     city=cities[n_c]
#     for n_l in eachindex(leadtimes)
#         IETM_24h[n_c,n_l]=vec(nanmean(read(IETM3D_data,Times=tStart_all_obs .- Hour(24),leadTimes=24,varNames="PM10",staIDs=get_staIDs(city)),dims=3))
#     end
# end

# RMSE_IETM=Array{Float64,2}(undef,length(cities),length(leadtimes))
# R_IETM=Array{Float64,2}(undef,length(cities),length(leadtimes))
# RMSE_IETM24=Array{Float64,2}(undef,length(cities),length(leadtimes))
# R_IETM24=Array{Float64,2}(undef,length(cities),length(leadtimes))
# for n_c in eachindex(cities)
#     city=cities[n_c]
#     for n_l in eachindex(leadtimes)
#         l=leadtimes[n_l]
#         RMSE_IETM[n_c,n_l]=getStatistics(obs[n_c,n_l], IETM[n_c,n_l])[1]
#         R_IETM[n_c,n_l]=getStatistics(obs[n_c,n_l], IETM[n_c,n_l])[3]
#         RMSE_IETM24[n_c,n_l]=getStatistics(obs[n_c,n_l], IETM_24h[n_c,n_l])[1]
#         R_IETM24[n_c,n_l]=getStatistics(obs[n_c,n_l], IETM_24h[n_c,n_l])[3]
#     end
# end

obs=Array{Array{Float64}}(undef,length(cities))
for n_c in eachindex(cities)
    city=cities[n_c]
    obs[n_c]=vec(nanmean(read(staObs,Times=tStart_all_obs,varNames="PM10",staIDs=get_staIDs(city)),dims=3))
end
#多个单次预报拼接读取
leadtimes2plot=[1,6,12,24]
leadtimes24=leadtimes[1:24]
cities2plot=["太原","北京","南京","海口"]
# IETM_Nlead = Array{Array{Float64},2}(undef,length(cities2plot),length(leadtimes2plot))
IETM_Nlead = Array{Array{Float64},2}(undef,length(cities),length(leadtimes24))
all_start_times = Array{Vector{ZonedDateTime}, 2}(undef, length(cities), length(leadtimes24))
total_len = length(tStart_all_obs)
fixed_start = tStart_all_obs[1] - Hour(1)
@time for n_l in eachindex(leadtimes24)
    leadtime = leadtimes24[n_l]
    n_steps = ceil(Int, total_len / leadtime)
    for n_c in eachindex(cities)
        city = cities[n_c]
        combined_data = Float64[]
        sizehint!(combined_data, total_len) 
        current_city_times = ZonedDateTime[]
        sizehint!(current_city_times, n_steps)
        # current_time = tStart_One - Hour(leadtime)
        current_time = fixed_start
        for n_t in 1:n_steps
            push!(current_city_times, current_time)
            # 读取数据
            data_read = read(IETM3D_data, 
                             Times = current_time, 
                             varNames = "PM10", 
                             staIDs = get_staIDs(city))
            
            IETM_OneForecast = vec(nanmean(data_read, dims=3))
            current_filled = length(combined_data)
            remaining_len = total_len - current_filled
            actual_take = min(leadtime, remaining_len)
            append!(combined_data, IETM_OneForecast[1:actual_take])
            current_time += Hour(leadtime)
        end
        IETM_Nlead[n_c, n_l] = combined_data
        all_start_times[n_c, n_l] = current_city_times
    end
end
all_start_times_format = Array{Vector{DateTime}, 2}(undef, length(cities), length(leadtimes24))
for n_c in eachindex(cities)
    for n_l in eachindex(leadtimes24)
        all_start_times_format[n_c,n_l] = DateTime.(all_start_times[n_c,n_l])
    end
end
RMSE_IETM=Array{Float64,2}(undef,length(cities),length(leadtimes24))
R_IETM=Array{Float64,2}(undef,length(cities),length(leadtimes24))
for n_c in eachindex(cities)
    city=cities[n_c]
    for n_l in eachindex(leadtimes24)
        # l=leadtimes24[n_l]
        # leadtime_index = indexin([l],leadtimes)[1]
        RMSE_IETM[n_c,n_l]=getStatistics(obs[n_c], IETM_Nlead[n_c,n_l])[1]
        R_IETM[n_c,n_l]=getStatistics(obs[n_c], IETM_Nlead[n_c,n_l])[3]
    end
end
RMSE_24diff6=RMSE_IETM[:,indexin([6],leadtimes24)] .- RMSE_IETM[:,indexin([24],leadtimes24)]
RMSE_24diff12=RMSE_IETM[:,indexin([12],leadtimes24)] .- RMSE_IETM[:,indexin([24],leadtimes24)]
length(findall(RMSE_24diff6 .< -100)) / length(RMSE_IETM[:,1])
include("/sh2/home/3clear/wuhj/paper_IETMsand/DFCI.jl")
#dtw
# IETM_dtw=Array{Float64,2}(undef,length(cities),length(leadtimes))
# IETM24_dtw=Array{Float64,2}(undef,length(cities),length(leadtimes))
# for n_c in eachindex(cities)
#     city=cities[n_c]
#     for n_l in eachindex(leadtimes)
#         l=leadtimes[n_l]
#         IETM_dtw[n_c,n_l]=dust_forecast_index(obs[n_c],IETM[n_c,n_l],300.0)[2]
#         IETM24_dtw[n_c,n_l]=dust_forecast_index(obs[n_c],IETM_24h[n_c,n_l],300.0)[2]
#     end
# end
IETM_dtw=Array{Float64,2}(undef,length(cities),length(leadtimes2plot))
for n_c in eachindex(cities)
    city=cities[n_c]
    for n_l in eachindex(leadtimes2plot)
        l=leadtimes2plot[n_l]
        IETM_dtw[n_c,n_l]=dust_forecast_index(obs[n_c],IETM_Nlead[n_c,n_l],300.0)[2]
    end
end
end

include(joinpath(@__DIR__, "PaperFigure_RUFimproveAccuracy_data.jl"))

begin
    t2plot=DateTime.(tStart_all_obs)
    plt=get_plt()
    # Tight bounding boxes add a small amount beyond the nominal canvas; this
    # width keeps the exported PDF within a 7-inch double-column layout.
    fig=plt.figure(figsize=(6.55,5.45),dpi=600);plt.clf()
    font_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_prop = fm.FontProperties(fname=font_path,size=TICK_SIZE)
    font_pathbd = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    font_propbd = fm.FontProperties(fname=font_pathbd,size=PANEL_TITLE_SIZE)
    font_axis = fm.FontProperties(fname=font_path,size=AXIS_LABEL_SIZE)
    font_annotation = fm.FontProperties(fname=font_path,size=ANNOTATION_SIZE)
    font_legend = fm.FontProperties(fname=font_path,size=FIGURE_LEGEND_SIZE)
    cmapM=get_truncated_cmap("seismic",0.15,0.85)
    rmse_boundaries = [-150.0, -100.0, -50.0, 0.0, 50.0, 100.0, 150.0]
    rmse_ncolors = length(rmse_boundaries) - 1
    cmapM_discrete = mcolors.ListedColormap(
        cmapM(np.linspace(0.0, 1.0, rmse_ncolors)),
        name="trunc_seismic_discrete",
    )
    rmse_norm = mcolors.BoundaryNorm(rmse_boundaries, cmapM_discrete.N, clip=false)
    proj=ccrs.LambertConformal(central_longitude=105, central_latitude=83,
    standard_parallels=[25,45],globe=nothing)
    ax2 = fig.add_axes([0.08,0.718,0.26,0.204])
    # ax2 = fig.add_axes([0.4,0.6,0.25,0.4],projection=proj)
    ax1 = fig.add_axes([0.395,0.71,0.26,0.22],projection=proj)
    ax3 = fig.add_axes([0.71,0.71,0.26,0.22],projection=proj)
    ax1.set_extent([72.3, 136.3, 16.9, 54.2],crs=PlaCar)
    # ax2.set_extent([72.3, 136.3, 16.9, 54.2],crs=PlaCar)
    ax3.set_extent([72.3, 136.3, 16.9, 54.2],crs=PlaCar)
    cn_boundary=ax1.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    # cn_boundary=ax2.add_geometries(shpreader.Reader("/sh2/home/3clear/wuhj/software/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    # PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    cn_boundary=ax3.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)

    # color_ax2=plt.get_cmap(cmapM)
    # color_rmse = color_ax2(0.2)
    # color_dtw  = color_ax2(0.8)

    rmse_mean = vec(nanmean(RMSE_IETM[:, 1:24], dims=1))
    r_mean  = vec(nanmean(R_IETM[:, 1:24],  dims=1))
    
    w = 0.35
    ax2.bar(leadtimes[1:24] .- w/2, rmse_mean;
            width = w,
            facecolor = OBS_COLOR,
            label = "RMSE")
    ax2.set_ylabel(raw"μg m$^{-3}$",labelpad=0.8,fontproperties=font_axis)
    ax2.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    ax2.tick_params(axis="x", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    ax2.set_xlim(0,25)
    ax2.set_xticks(0:6:24)  
    ax2.set_ylim(0,160)
    ax2r = ax2.twinx()
    ax2r.bar(leadtimes[1:24] .+ w/2, r_mean;
            width = w,
            facecolor = RUF_COLOR,
            alpha = 0.7,
            label = "R")
    # ax2r.set_ylabel(size=5.5,labelpad=0.5,fontproperties=font_prop)
    ax2r.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    ax2r.set_xlim(0,25)
    ax2r.set_ylim(0,1)
    ax2.set_xlabel("Update frequency (h)",labelpad=1,fontproperties=font_axis)
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2r.get_legend_handles_labels()
    # ax2.legend(vcat(h1, h2), vcat(l1, l2), loc="upper center",prop=font_prop)
    leg = ax2r.legend(vcat(h1, h2), vcat(l1, l2), 
    loc="lower center",          
    bbox_to_anchor=(0.5, 0.05),   
    ncol=2, 
    prop=font_legend,
    frameon=true                
    )
    leg.set_zorder(10)

    cities_lon=[nanmean(Float64.(staInfo.lon[staInfo.city .== city])) for city in cities]
    cities_lat=[nanmean(Float64.(staInfo.lat[staInfo.city .== city])) for city in cities]


    RMSE_24diff6=RMSE_IETM[:,indexin([6],leadtimes24)] .- RMSE_IETM[:,indexin([24],leadtimes24)]
    RMSE_24diff12=RMSE_IETM[:,indexin([12],leadtimes24)] .- RMSE_IETM[:,indexin([24],leadtimes24)]
    # R_24diff6=R_IETM[:,4] .- R_IETM[:,2]
    # dtw_diff=IETM_dtw[:,4] .- IETM_dtw[:,2]

    # RMSE_diff2plot=nanmean(RMSE_diff,dims=2)
    # R_diff2plot=nanmean(R_diff,dims=2)
    # dtw_diff2plot=nanmean(dtw_diff,dims=2)
    RMSE_diff2plot_24diff6=RMSE_24diff6
    RMSE_diff2plot_24diff12=RMSE_24diff12

    # mean_RMSE_IETM=nanmean(RMSE_IETM[:,1:24],dims=2)
    # mean_R_IETM=nanmean(R_IETM[:,1:24],dims=2)
    # mean_dtw_IETM=nanmean(IETM_dtw[:,1:24],dims=2)
    # mean_RMSE_IETM24=nanmean(RMSE_IETM24[:,1:24],dims=2)
    # mean_R_IETM24=nanmean(R_IETM24[:,1:24],dims=2)
    # mean_dtw_IETM24=nanmean(IETM24_dtw[:,1:24],dims=2)
    # RMSE_diff2plot=mean_RMSE_IETM .- mean_RMSE_IETM24
    # R_diff2plot=mean_R_IETM .- mean_R_IETM24
    # dtw_diff2plot=mean_dtw_IETM .- mean_dtw_IETM24


    addGrids(ax1,false)
    add_SouthChinaAxes(ax1,cities_lon,cities_lat,RMSE_diff2plot_24diff6,cmapM_discrete,rmse_norm)
    # addGrids(ax2,false)
    # add_SouthChinaAxes(ax2,cities_lon,cities_lat,R_diff2plot,-1,1)
    addGrids(ax3,true)
    add_SouthChinaAxes(ax3,cities_lon,cities_lat,RMSE_diff2plot_24diff12,cmapM_discrete,rmse_norm)

    sc_ax1=ax1.scatter(cities_lon,cities_lat,c=RMSE_diff2plot_24diff6,1,
        cmap=cmapM_discrete,norm=rmse_norm,transform=PlaCar)
    # sc_ax2=ax2.scatter(cities_lon,cities_lat,c=R_diff2plot,1,cmap=cmapM,transform=PlaCar,vmin=-1,vmax=1)
    sc_ax3=ax3.scatter(cities_lon,cities_lat,c=RMSE_diff2plot_24diff12,1,
        cmap=cmapM_discrete,norm=rmse_norm,transform=PlaCar)

    cax_ax1 = fig.add_axes([0.395,0.65,0.26,0.012])
    # cax_ax2 = fig.add_axes([0.4,0.65,0.25,0.01])
    cax_ax3 = fig.add_axes([0.71,0.65,0.26,0.012])


    cbar_ax1=fig.colorbar(sc_ax1,cax_ax1,orientation="horizontal",extend="both",
        boundaries=rmse_boundaries,ticks=rmse_boundaries,spacing="uniform",drawedges=true)
    # cbar_ax2=fig.colorbar(sc_ax2,cax_ax2,orientation="horizontal",extend="both")
    cbar_ax3=fig.colorbar(sc_ax3,cax_ax3,orientation="horizontal",extend="both",
        boundaries=rmse_boundaries,ticks=rmse_boundaries,spacing="uniform",drawedges=true)


    cbar_ax1.ax.tick_params(axis="both",which="both",length=2.8,labelsize=TICK_SIZE,pad=2,width=AXIS_LINEWIDTH)
    # cbar_ax2.ax.tick_params(axis="both",which="both",length=2,labelsize=5.5,pad=2)
    cbar_ax3.ax.tick_params(axis="both",which="both",length=2.8,labelsize=TICK_SIZE,pad=2,width=AXIS_LINEWIDTH)
    ax1.tick_params(axis="x",pad=2)
    ax1.tick_params(axis="y",pad=2)
    # ax2.tick_params(axis="x",pad=2)
    # ax2.tick_params(axis="y",pad=2)

    add_label_atleftup(ax1,"(b) ΔRMSE (U6h-U24h)")
    add_label_atleftup(ax2,"(a) Accuracy Changes")
    add_label_atleftup(ax3,"(c) ΔRMSE (U12h-U24h)")

    ax4 = fig.add_axes([0.08, 0.42, 0.42, 0.16])  # 第一行
    ax5 = fig.add_axes([0.55, 0.42, 0.42, 0.16])  # 第二行
    ax6 = fig.add_axes([0.08, 0.23, 0.42, 0.16])  # 第三行
    ax7 = fig.add_axes([0.55, 0.23, 0.42, 0.16])  # 第四行
    
    
    n = length(leadtimes2plot)
    # cmapM2getcolor = plt.get_cmap(cmapM)
    cmapM2getcolor = get_truncated_cmap("afmhot", 0.25, 0.65)
    colors = n == 1 ? [cmapM2getcolor(0.5)] :
        [cmapM2getcolor((i - 1) / (n - 1) * 0.9) for i in 1:n]
    line_styles = fill("-", length(leadtimes2plot))
    line_markers = fill(nothing, length(leadtimes2plot))
    function draw_text(ax,city)
        metric_x = [0.53, 0.66, 0.79, 0.92]
        ax.text(0.34, 0.98, "RMSE:", transform=ax.transAxes, va="top",
            size=ANNOTATION_SIZE, fontproperties=font_annotation)
        for n_l in eachindex(leadtimes2plot)
            l=leadtimes2plot[n_l]
            ax.text(metric_x[n_l], 0.98,
                "$(round(RMSE_IETM[indexin([city],cities),indexin([l],leadtimes24)][1], digits=1))",
                transform=ax.transAxes, va="top", ha="center", c=colors[n_l],
                size=ANNOTATION_SIZE, fontproperties=font_annotation)
        end
        ax.text(0.34, 0.83, "R:", transform=ax.transAxes, va="top",
            size=ANNOTATION_SIZE, fontproperties=font_annotation)
        for n_l in eachindex(leadtimes2plot)
            l=leadtimes2plot[n_l]
            ax.text(metric_x[n_l], 0.83,
                "$(round(R_IETM[indexin([city],cities),indexin([l],leadtimes24)][1], digits=2))",
                transform=ax.transAxes, va="top", ha="center", c=colors[n_l],
                size=ANNOTATION_SIZE, fontproperties=font_annotation)
        end
        # ax.text(tText[1],ythirdline,"DTW:",size=4.5,fontproperties=font_prop)
        # for n_l in eachindex(leadtimes2plot)
        #     l=leadtimes2plot[n_l]
        #     ax.text(tText[1+n_l],ythirdline,"$(round(IETM_dtw[indexin([city],cities),indexin([l],leadtimes24)][1],
        #     digits=2))",c=colors[n_l],size=4.5,fontproperties=font_prop)
        # end
    end
    function draw_lineFigures_new(ax,t2plot,obs,IETM,city)
        t2plot_label=Dates.format.(tStart_all_obs,"mm-dd")
        city_idx = indexin([city], cities)[1]
        ax.bar(t2plot,obs[city_idx],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
        max_con=Float64[]
        for n_l in eachindex(leadtimes2plot)
            # l=leadtimes2plot[n_l]
            # l_idx = indexin([l], leadtimes24)[1]
            # if l == 12 || l == 24
            #     ax.scatter(all_start_times_format[city_idx,l_idx],IETM_Nlead[city_idx,l_idx][indexin(all_start_times_format[city_idx,l_idx] .+ Hour(l),t2plot)]
            #     ,3,color=colors[n_l],zorder=5)
            # end
            # ax.plot(t2plot,IETM[indexin([city],cities),indexin([l],leadtimes24)][1],color=colors[n_l],linewidth=1,label="RUF_$(l)h")
        
            l = leadtimes2plot[n_l]
            l_idx = indexin([l], leadtimes24)[1]
            first_forecast_times = all_start_times_format[city_idx, l_idx] .+ Hour(1)
            idxs_in_t = indexin(first_forecast_times, t2plot)
            valid_mask = .!isnothing.(idxs_in_t)
            
            if (l == 12) && any(valid_mask)
                plot_x = all_start_times_format[city_idx, l_idx][valid_mask]
                clean_idxs = Int.(idxs_in_t[valid_mask])
                # plot_y = zeros(length(plot_x))
                plot_y = IETM[city_idx, l_idx][clean_idxs]
                ax.scatter(plot_x, plot_y, s=2, color=colors[n_l], 
                           zorder=28, linewidths=0.5, label="Start_$(l)h")
            elseif (l == 24) && any(valid_mask)
                plot_x = all_start_times_format[city_idx, l_idx][valid_mask]
                clean_idxs = Int.(idxs_in_t[valid_mask])
                # plot_y = zeros(length(plot_x))
                plot_y = IETM[city_idx, l_idx][clean_idxs]
                ax.scatter(plot_x, plot_y, s=8, color=colors[n_l], 
                           zorder=16, linewidths=0.5, label="Start_$(l)h")
            end
            Mlinewidth = l == 1 || l == 6 ? 0.6 : DATA_LINEWIDTH
            Mzorder = l == 1 || l == 6 ? l : 40-l
            ax.plot(t2plot, IETM[city_idx, l_idx], color=colors[n_l], 
                    linewidth=Mlinewidth, linestyle=line_styles[n_l], marker=line_markers[n_l],
                    markevery=24, markersize=2.6, label="RUF_U$(l)h", alpha=1,zorder=Mzorder)
            
            push!(max_con,nanmaximum(IETM[city_idx, l_idx]))
        end
        # step = 72
        # xticks = t2plot[1:step:end]
        # xticklabels = t2plot_label[1:step:end]
        daily_step = 24
        all_daily_ticks = t2plot[1:daily_step:end]
        xticklabels_new = String[]
        for (i, t) in enumerate(all_daily_ticks)
            if (i - 1) % 3 == 0
                push!(xticklabels_new, Dates.format(t, "mm-dd"))
            else
                push!(xticklabels_new, "")
            end
        end

        ax.set_xticks(all_daily_ticks)
        ax.set_xticklabels(xticklabels_new)
        # ax.set_xticks(xticks)
        # ax.set_xticklabels(xticklabels)
        ax.set_xlim(DateTime(2025, 4, 11, 0, 0), DateTime(2025, 4, 21))
        ax.set_ylim(0,maximum(max_con) * 1.3)
        ax.tick_params(axis="x", pad=2)
        ax.tick_params(axis="y", pad=1.5)
        ax.tick_params(axis="x", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
        ax.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
        ax.tick_params(axis="x", zorder=100)
        ax.tick_params(axis="y", zorder=100)
    end 
    
    axesM=[ax4,ax5,ax6,ax7]
    for n_c in eachindex(cities2plot)
        city=cities2plot[n_c]
        draw_lineFigures_new(axesM[n_c],t2plot,obs,IETM_Nlead,city)
        draw_text(axesM[n_c],city)
    end
    add_label_atleftup(ax4,"(d) Taiyuan"; ypos=0.99)
    add_label_atleftup(ax5,"(e) Beijing"; ypos=0.99)
    add_label_atleftup(ax6,"(f) Nanjing"; ypos=0.99)
    add_label_atleftup(ax7,"(g) Haikou"; ypos=0.99)
    ax4.set_ylabel(raw"$\mathrm{PM}_{10}$"*raw" (μg m$^{-3}$)",labelpad=0.8,fontproperties=font_axis)
    ax6.set_ylabel(raw"$\mathrm{PM}_{10}$"*raw" (μg m$^{-3}$)",labelpad=0.8,fontproperties=font_axis)
    handles1, labels1 = ax4.get_legend_handles_labels()

    # 2. 定义分类逻辑
    # 我们把图例分为三类：柱状图(Bar/Obs)、折线(RUF)、起报点(Start)
    bar_idx = Int[]
    line_idx = Int[]
    point_idx = Int[]
    other_idx = Int[]

    for (i, label) in enumerate(labels1)
        if startswith(label, "RUF")
            push!(line_idx, i)
        elseif startswith(label, "Start")
            push!(point_idx, i)
        else
            push!(bar_idx, i) 
        end
    end

    get_num(s) = (m = match(r"\d+", s); m === nothing ? 0 : parse(Int, m.match))

    sort!(line_idx, by = i -> get_num(labels1[i]))
    sort!(point_idx, by = i -> get_num(labels1[i]))
    final_indices = [bar_idx; line_idx; point_idx]

    ordered_handles = handles1[final_indices]
    ordered_labels = labels1[final_indices]

    fig.legend(ordered_handles, ordered_labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.115),
        prop=font_legend,
        ncol=4,
        frameon=true,
        handletextpad=0.4,
        columnspacing=0.8
    )
    plt.gcf()
end
output_dir = "/home/zsm/IETM/figure/PapreFigure/20260903"
mkpath(output_dir)
plt.savefig(joinpath(output_dir, "PaperFigure_RUFimproveAccuracy.png"), bbox_inches="tight", pad_inches=0.01,dpi=600)
plt.savefig(joinpath(output_dir, "PaperFigure_RUFimproveAccuracy.pdf"), bbox_inches="tight", pad_inches=0.01)
plt.close("all")
exit()
#备选 
#景德镇 萍乡




maximum(IETM_Nlead[1, l]) * 1.1

IETM_Nlead[indexin(["北京"],cities)
,indexin([6],leadtimes24)][1][indexin(all_start_times_format[indexin(["北京"],cities)[1]
,indexin([6],leadtimes24)[1]] .+ Hour(6),t2plot)]

RMSE_diff=RMSE_IETM[:,1:24] .- RMSE_IETM24[:,1:24]
R_diff=R_IETM[:,1:24] .- R_IETM24[:,1:24]
dtw_diff=IETM_dtw[:,1:24] .- IETM24_dtw[:,1:24]
RMSE_diff2plot=nanmean(RMSE_diff,dims=2)
R_diff2plot=nanmean(R_diff,dims=2)
dtw_diff2plot=nanmean(dtw_diff,dims=2)

# 统计dtw差异小于0并且rmse差异大于 0 的元素
indices = findall((RMSE_diff2plot .> 0) .& (dtw_diff2plot .< 0) .& (R_diff2plot .< 0))
cities_spatial = cities[indices]
# 计算个数
count_negative = length(indices)

# 输出结果
println("Number of elements less than 0: ", count_negative)
println("Indices of elements less than 0: ", indices)
index=indices[13][1]
city_curr=cities[index]
plt=get_plt()
fig=plt.figure(figsize=(6,5),dpi=600);plt.clf()
plt.bar(t2plot,obs[index,15],facecolor="deepskyblue",width=0.02,label="obs")
plt.plot(t2plot,vec(nanmean(read(IETM3D_data,Times=tStart_all_obs,leadTimes=0,varNames="PM10",staIDs=get_staIDs(city_curr)),dims=3)),
label="RUF0")
plt.plot(t2plot,IETM[index,1],label="RUF1")
plt.plot(t2plot,IETM[index,6],label="RUF6")
plt.plot(t2plot,IETM[index,18],label="RUF18")
plt.plot(t2plot,IETM[index,20],label="RUF20")
plt.plot(t2plot,IETM[index,24],label="RUF24")


plt.legend()
plt.gcf()
println("RMSE: ", RMSE_IETM[index,15])
println("R: ", R_IETM[index,15])
println("RMSE_diff: ", RMSE_diff2plot[index])
println("R_diff: ", R_diff2plot[index])
plt.close()
