using Revise,MyPlotM,StaDataM,StaDataTool,utilities,TimeZones
using PlotFunM,Dates,NaNStatistics
# using DBInterface,DataFrames
using PyCall 
using LinearAlgebra
using Statistics
using PyPlot
using DataFrames
include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle
patches = pyimport("matplotlib.patches")
matplotlib = pyimport("matplotlib")
configure_matplotlib!(matplotlib)
plt = matplotlib.pyplot
ticker=matplotlib.ticker
font_manager=matplotlib.font_manager
cm=matplotlib.cm
ccrs = pyimport("cartopy.crs")
shpreader = pyimport("cartopy.io.shapereader")
fm=matplotlib.font_manager
mcolors = pyimport("matplotlib.colors")
np = pyimport("numpy")
include(joinpath(@__DIR__, "PaperFigure_Earlywarningperformance_data.jl"))


function readData(data,Time,leadtime,varname,city)
    result=vec(nanmean(read(data,Times=Time.-Hour(leadtime),leadTimes=leadtime,varNames=varname,staIDs=staInfo.staID[staInfo.city .== city]),dims=3))
    return result
end
function readData(data,Time,varname,city)
    mask = figure3_city.(data.staInfo.city) .== city
    result=vec(nanmean(read(data,Times=Time,varNames=varname,staIDs=data.staInfo.staID[mask]),dims=3))
    return result
end
function add_label_atleftup(ax, label; xpos=0.02, ypos=0.95, fontsize=PANEL_TITLE_SIZE)
    ax.text(xpos, ypos, label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            va="top",
            ha="left",
            fontproperties=font_propbd)
end
function addGrids(ax,bool_right,bool_bottom=false)
    gl = ax.gridlines(crs=PlaCar,draw_labels=true, dms=true, x_inline=false, y_inline=false, linewidth=GRID_LINEWIDTH, color=[0.1,0.1,0.1], alpha=0.5,zorder=2,linestyle="--") # 添加经纬度网格，这里设定为黑色虚线半透明
    gl.top_labels  = false;    gl.right_labels = bool_right;   gl.rotate_labels = false ;gl.left_labels  = false;
    gl.bottom_labels = bool_bottom
    
    gl.xlocator=ticker.FixedLocator([40,60,80,100,120,140,160,180])
    gl.ylocator=ticker.FixedLocator([20,30,40,50,60,70])

    gl.xlabel_style = Dict(:size => TICK_SIZE, :fontname => "Arial")
    gl.ylabel_style = Dict(:size => TICK_SIZE, :fontname => "Arial")
    gl.xpadding = 1   # 经度标签离图框的距离
end 

function add_SouthChinaAxes_zsm(axIn)
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
    # ax_southChina.scatter(cities_lon,cities_lat,color="r",marker="o",1,transform=PlaCar)
    ax_southChina.scatter(sta_lon_nonzero,sta_lat_nonzero,c=advancetimes_nonzero,1,cmap=cmap_scatter,norm=normM,transform=PlaCar)
    return ax_southChina
end
function connect_map_to_ax!(fig, ax_center, lon, lat, ax_target;
    crs, at=(0.5,0.5), color="grey", lw=0.6, style="-|>",ls="--", ms=10, clip_on=false
)
    transA = crs[:_as_mpl_transform](ax_center)   # 允许 (lon,lat) 作为输入
    con = patches.ConnectionPatch(
        (lon, lat), at,
        transA, "axes fraction",
        axesA=ax_center, axesB=ax_target,
        arrowstyle=style, linestyle=ls,
        mutation_scale=ms,
        lw=lw, color=color, shrinkB=2,  # shrinkB 让箭头稍离子图边框
        clip_on=clip_on
    )
    fig.add_artist(con)
    return con
end

# 便捷：根据目标子图在“上/下/左/右”决定终点 at
function edge_anchor(side::Symbol)
    side === :top    && return (0.5, 0.0)  # 子图下边中点
    side === :bottom && return (0.5, 1.0)  # 子图上边中点
    side === :left   && return (1.0, 0.5)  # 子图右边中点
    side === :right  && return (0.0, 0.5)  # 子图左边中点
    return (0.5, 0.5)
end

function plot1ax(ax,n_ax,obs,IETM,leadtimes,city)
    n = length(leadtimes)
    cmapM = get_truncated_cmap("afmhot", 0.25, 0.7)
    colors = n == 1 ? [cmapM(0.5)] : [cmapM((i - 1) / (n - 1)) for i in 1:n]
    line_styles = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
    markers = [nothing, "o", "s", "^", "D"]
    ax.bar(t2plot,readData(obs,tStart_all_obs,var2plot,city),color=NEUTRAL_COLOR,0.02,label="obs")
    for i in 1:length(leadtimes)
        ax.plot(t2plot,readData(IETM,tStart_all_obs,leadtimes[i],var2plot,city),
            color=colors[i], linewidth=DATA_LINEWIDTH, linestyle=line_styles[i],
            marker=markers[i], markevery=12, markersize=2.6, label="RUF_L$(leadtimes[i])h")
    end
    ax.set_xlim(DateTime(2025, 4, 11, 0, 0), DateTime(2025, 4, 14))
    ax.set_xticks(t2plot[1:72:end]) 
    ax.set_xticklabels(t2plot_label[1:72:end], fontsize=TICK_SIZE,fontproperties=font_prop)
    xticks = ax.get_xticklabels()
    # 最左边的
    xticks[1].set_ha("left")    # 左对齐 → 整体往右靠
    # 最右边的
    xticks[end].set_ha("right")  # 右对齐 → 整体往左靠
    ax.tick_params(axis="x", labelsize=TICK_SIZE,pad=2,width=AXIS_LINEWIDTH)
    ax.tick_params(axis="y", labelsize=TICK_SIZE,pad=2,width=AXIS_LINEWIDTH)
    idx = findfirst(==(city), cities)  
    city_name = cities_Eng[idx]
    letter = Char('a' + n_ax) 
    full_label = "($letter) $city_name"
    add_label_atleftup(ax,full_label)
end
# The legacy multi-file preparation is retained for reference but replaced by
# figure3_v2_newdata_prep.jl above.
if false
#突变点识别
include("/sh2/home/3clear/wuhj/paper_IETMsand/change_point_detect.jl")
#初始化
leadtimes=collect(1:72)
cities2detect = unique(staObs.staInfo.city)
derivative_threshold = 0.1
concentration_threshold = 200.0
obs_point=Array{Array{Float64}}(undef,length(cities2detect))
ietm_point=Array{Array{Float64},2}(undef,length(cities2detect),length(leadtimes))
obs_normal_der=Array{Array{Float64}}(undef,length(cities2detect))
ietm_normal_der=Array{Array{Float64},2}(undef,length(cities2detect),length(leadtimes))
bool_ietm=Array{Int64}(undef,length(cities2detect),length(leadtimes))
#读取模拟及观测数据
obs_array=Array{Array{Float64}}(undef,length(cities2detect))
for n_c in eachindex(cities2detect)
    city=cities2detect[n_c]
    # mask=isequal.(staObs.staInfo.city,city) .& isequal.(staObs.staInfo.type,"B")
    mask=isequal.(staObs.staInfo.city,city)
    staids=staObs.staInfo[mask,:staID]
    obs_array[n_c]=vec(nanmean(read(staObs,Times=tStart_all_obs,varNames="PM10",staIDs=staids),dims=3))
end
# data=read(ietm_dataAll,Times=tStart_all_obs,varNames=varName,staIDs=IDs)
ietm3d_array=Array{Array{Float64},2}(undef,length(cities2detect),length(leadtimes))
for n_c in eachindex(cities2detect)
    for n_l in eachindex(leadtimes)
        l=leadtimes[n_l]
        city=cities2detect[n_c]
        # mask=isequal.(staObs.staInfo.city,city) .& isequal.(staObs.staInfo.type,"B")
        mask=isequal.(staObs.staInfo.city,city)
        staids=staObs.staInfo[mask,:staID]
        ietm3d_array[n_c,n_l]=vec(nanmean(read(IETMdata,Times=tStart_all_obs .- Hour(l),leadTimes=l,varNames="PM10",staIDs=staids),dims=3))
    end
end

for n_c in eachindex(cities2detect)
    obs_point[n_c],obs_normal_der[n_c]=detect_starting_points(obs_array[n_c],derivative_threshold,concentration_threshold);
    for n_l in eachindex(leadtimes)
        ietm_point[n_c,n_l],ietm_normal_der[n_c,n_l]=detect_starting_points(ietm3d_array[n_c,n_l],derivative_threshold,concentration_threshold);
    end
end
for n_c in eachindex(cities2detect)
    for n_l in eachindex(leadtimes)

        bool_ietm[n_c,n_l]=judge_forcast(obs_array[n_c],ietm3d_array[n_c,n_l],Int.(obs_point[n_c]),Int.(ietm_point[n_c,n_l]),6);
    end
end
advance_times=bool2leadtime(bool_ietm)
cities=[
    
    "呼和浩特",
    "北京",
    "太原",
    # "长治",
    # "晋城",
    "洛阳",
    # "南阳",
    "襄阳",
    "宜昌",
    "常德",
    # "清远",
    "桂林",
    "茂名",
    "海口"
    # "永州",
    # "邵阳"
]
index_evalution=indexin(cities,cities2detect)
first_point_ietm=Array{Float64,2}(undef,length(cities),12)
first_point_obs=Array{Float64}(undef,length(cities))
for n_c in eachindex(cities)
    first_point_obs[n_c]=obs_point[index_evalution][n_c][1]
    for n_l in 1:12
        if length(ietm_point[index_evalution,12:24][n_c,n_l]) == 0
            first_point_ietm[n_c,n_l] = NaN
        else
            first_point_ietm[n_c,n_l]=ietm_point[index_evalution,12:24][n_c,n_l][1]
        end
    end
end
obs_point[index_evalution]
nanmean(abs.(first_point_ietm .- first_point_obs))
a=nanmean(abs.(first_point_ietm .- first_point_obs),dims=2)
minimum(a)
end

cities=[
    "呼和浩特",
    "北京",
    "太原",
    "洛阳",
    "襄阳",
    "宜昌",
    "常德",
    "桂林",
    "茂名",
    "海口"
]

# index_evalution=indexin(cities,cities2detect)
# first_point_ietm=Array{Float64,2}(undef,length(cities),12)
# first_point_obs=Array{Float64}(undef,length(cities))
# for n_c in eachindex(cities)
#     first_point_obs[n_c]=obs_point[n_c][1]
#     for n_l in 1:12
#         if length(ietm_point[:,12:24][n_c,n_l]) == 0
#             first_point_ietm[n_c,n_l] = NaN
#         else
#             first_point_ietm[n_c,n_l]=ietm_point[:,12:24][n_c,n_l][1]
#         end
#     end
# end
# nanmean(abs.(first_point_ietm .- first_point_obs))
# maximum(abs.(filter(!isnan, first_point_ietm .- first_point_obs)))
# minimum(abs.(filter(!isnan, first_point_ietm .- first_point_obs)))
cities_Eng=[
    
    "Hohhot",
    "Beijing",
    "Taiyuan",
    # "Changzhi",
    # "Jincheng",
    "Luoyang",
    # "Nanyang",
    "Xiangyang",
    "Yichang",
    "Changde",
    # "Qingyuan",#yichang
    "Guilin",
    "Maoming",
    "Haikou"
    # "Yongzhou",
    # "Shaoyang"
    
    

]
#绘图
function get_truncated_cmap(base_cmap_name, min_val=0.0, max_val=1.0, n=256)
    base_cmap = plt.get_cmap(base_cmap_name)
    # 在 min_val 到 max_val 之间重采样
    colors_resampled = base_cmap(np.linspace(min_val, max_val, n))
    # 创建新的颜色条对象
    new_cmap = mcolors.LinearSegmentedColormap.from_list("trunc_afmhot", colors_resampled)
    return new_cmap
end
begin
    
    t2plot_label=Dates.format.(tStart_all_obs,"mm-dd")
    t2plot=DateTime.(tStart_all_obs)
    staInfo_cities=staInfo[indexin(cities,city_labels),:]
    cities_lon=staInfo_cities.lon
    cities_lat=staInfo_cities.lat
    leadtimes2plot=[1,3,6,12,24]
    plt=get_plt()
    proj=ccrs.LambertConformal(central_longitude=105, central_latitude=83,standard_parallels=[25,45],globe=nothing)
    font_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_prop = fm.FontProperties(fname=font_path,size=TICK_SIZE)
    font_pathbd = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    font_propbd = fm.FontProperties(fname=font_pathbd,size=PANEL_TITLE_SIZE)
    font_axis = fm.FontProperties(fname=font_path,size=AXIS_LABEL_SIZE)
    font_legend = fm.FontProperties(fname=font_path,size=FIGURE_LEGEND_SIZE)
    fig=plt.figure(figsize=(6.30,3.45),dpi=600);plt.clf()
    #行_列
    ax1_1 = fig.add_axes([0.045,0.82,0.17,0.135])
    # ax1_2 = fig.add_axes([0.25,0.8,0.15,0.15])
    # ax1_3 = fig.add_axes([0.44,0.8,0.15,0.15])
    # ax1_4 = fig.add_axes([0.63,0.8,0.15,0.15])
    ax1_5 = fig.add_axes([0.785,0.82,0.17,0.135])

    ax2_1 = fig.add_axes([0.045,0.635,0.17,0.135])
    ax2_5 = fig.add_axes([0.785,0.635,0.17,0.135])
    ax3_1 = fig.add_axes([0.045,0.45,0.17,0.135])
    ax3_5 = fig.add_axes([0.785,0.45,0.17,0.135])
    ax4_1 = fig.add_axes([0.045,0.265,0.17,0.135])
    ax4_5 = fig.add_axes([0.785,0.265,0.17,0.135])


    ax5_1 = fig.add_axes([0.045,0.08,0.17,0.135])
    # ax5_2 = fig.add_axes([0.25,0.02,0.15,0.15])
    # ax5_3 = fig.add_axes([0.44,0.02,0.15,0.15])
    # ax5_4 = fig.add_axes([0.63,0.02,0.15,0.15])
    ax5_5 = fig.add_axes([0.785,0.08,0.17,0.135])

    # Keep the central map clear of the two city-series columns.
    ax_center = fig.add_axes([0.225,0.280,0.55,0.7], projection=proj)
    ax_center.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
        PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    ax_center.set_extent([72.3, 136.3, 16.9, 54.2],crs=PlaCar)

    bounds = [0, 3, 6, 9, 12, 24, 36, 48, 60, 72]
    # cmap_scatter = get_truncated_cmap("viridis", 0.1, 0.9)
    cmap_scatter = plt.get_cmap("viridis_r")
    normM = mcolors.BoundaryNorm(bounds, cmap_scatter.N) 
    addGrids(ax_center,false,false)
    sta_lon=combine(groupby(staInfo, :city), :lon => mean => :avg_lon).avg_lon
    sta_lat=combine(groupby(staInfo, :city), :lat => mean => :avg_lat).avg_lat
    mask_zero = advance_times .== 0 # 创建布尔掩码标记0值的位置
    sta_lon_nonzero=sta_lon[.!mask_zero]
    sta_lat_nonzero=sta_lat[.!mask_zero]
    # unique_cities_nonzero=unique_cities[.!mask_zero]
    advancetimes_nonzero = advance_times[.!mask_zero]
    axx=ax_center.scatter(sta_lon_nonzero,sta_lat_nonzero,c=advancetimes_nonzero,2,cmap=cmap_scatter,norm=normM,transform=PlaCar)
    add_SouthChinaAxes_zsm(ax_center)
    add_label_atleftup(ax_center,"(a) Early warning performance of RUF")


    # axs=[ax1_1,ax1_2,ax1_3,ax1_4,ax1_5,ax2_5,ax3_5,ax4_5,ax5_5,ax5_4,ax5_3,ax5_2,ax5_1,
    # ax4_1,ax3_1,ax2_1]
    axs=[ax1_1,ax2_1,ax3_1,ax4_1,ax5_1,ax1_5,ax2_5,ax3_5,ax4_5,ax5_5]
    # axs=reverse(axs)
    n_total = length(axs)
    half_point = ceil(Int, n_total / 2)
    for n_ax in eachindex(axs)
        ax=axs[n_ax]
        city=cities[n_ax]
        plot1ax(ax,n_ax,staObs,IETMdata,leadtimes2plot,city)
        # if n_ax <= 5
        #     connect_map_to_ax!(fig, ax_center, cities_lon[n_ax], cities_lat[n_ax], ax;
        #     crs=PlaCar, at=edge_anchor(:top),
        #     color="k", lw=0.4, style="-", ms=12)
        # elseif n_ax >5 && n_ax <=9 
        #     connect_map_to_ax!(fig, ax_center, cities_lon[n_ax], cities_lat[n_ax], ax;
        #     crs=PlaCar, at=edge_anchor(:right),
        #     color="k", lw=0.4, style="-", ms=12)
        # elseif n_ax >9 && n_ax <=13
        #     connect_map_to_ax!(fig, ax_center, cities_lon[n_ax], cities_lat[n_ax], ax;
        #     crs=PlaCar, at=edge_anchor(:bottom),
        #     color="k", lw=0.4, style="-", ms=12)
        # else
        #     connect_map_to_ax!(fig, ax_center, cities_lon[n_ax], cities_lat[n_ax], ax;
        #     crs=PlaCar, at=edge_anchor(:left),
        #     color="k", lw=0.4, style="-", ms=12)
        # end
        if n_ax <= half_point
            # --- 左侧子图：连接到子图的右中心 (1.0, 0.5) ---
            connect_map_to_ax!(fig, ax_center, cities_lon[n_ax], cities_lat[n_ax], ax;
                crs=PlaCar, 
                at=(1.0, 0.5), # 锚点在子图右侧
                color="grey", lw=0.6, style="-", ms=12)
        else
            # --- 右侧子图：连接到子图的左中心 (0.0, 0.5) ---
            connect_map_to_ax!(fig, ax_center, cities_lon[n_ax], cities_lat[n_ax], ax;
                crs=PlaCar, 
                at=(0.0, 0.5), # 锚点在子图左侧
                color="grey", lw=0.6, style="-", ms=12)
            ax.yaxis.tick_right()
            ax.tick_params(axis="y", labelleft=false, labelright=true)
        end
    end
    fig.text(-0.022, 0.52, raw"$\mathrm{PM}_{10}$"*raw" (μg m$^{-3}$)", rotation=90,
        va="center", ha="center", fontproperties=font_axis)
    
    # 取 ax1_1 的 legend 条目
    handles1, labels1 = ax1_1.get_legend_handles_labels()

    # 取 ax_center 的 legend 条目
    handles2, labels2 = ax_center.get_legend_handles_labels()

    # 拼接
    handles = vcat(handles1, handles2)
    labels  = vcat(labels1,  labels2)

    axouts=fig.add_axes([0.31, 0.225, 0.4, 0.02])
    cbar=fig.colorbar(
        axx,
        cax=axouts,
        orientation="horizontal",
        extend="max",
        ticks=bounds,
    )
    cbar.ax.set_title("Lead time (h)",fontsize=AXIS_LABEL_SIZE,pad=2,fontproperties=font_axis)
    cbar.ax.tick_params(labelsize=TICK_SIZE,length=2.8,pad=1,width=AXIS_LINEWIDTH)
    cbar.ax.set_zorder(10)  # 数字越大越靠上
    # # 在整个 figure 上画图例，位置随意指定
    fig.legend(handles, labels,
               loc="lower center",       # 位置
               bbox_to_anchor=(0.5, 0.055),  # 相对 figure 坐标
               prop=font_legend, ncol=3,
               frameon=true,
               )           # 是否带边框
    # scaleM = 1.5  # 放大 1.5 倍
    # sz = fig.get_size_inches()
    # fig.set_size_inches(sz[1]*scaleM, sz[2]*scaleM)
    plt.gcf()
end
output_dir = "/home/zsm/IETM/figure/PapreFigure/20260826"
mkpath(output_dir)
plt.savefig(joinpath(output_dir, "PaperFigure_Earlywarningperformance.png"),dpi=600,bbox_inches="tight")
plt.savefig(joinpath(output_dir, "PaperFigure_Earlywarningperformance.pdf"),bbox_inches="tight")
plt.close()
println("Saved PaperFigure_Earlywarningperformance.png and .pdf")
