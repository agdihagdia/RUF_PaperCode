using Dates,TimeZones,StaDataM,utilities,NaNStatistics,Statistics,DataFrames,PlotFunM
using JLD2,GradDataM
using Statistics
using PyCall
include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle
include(joinpath(@__DIR__, "PaperFigure_CompareNationAndAllsta_data.jl"))
if false
#初始化
staObs=staDATA("/sh2/home/3clear/wuhj/paper_IETMsand/ori202504.mat")
staInfo=staObs.staInfo
# mask = (staInfo.city = "北京市") .& (staInfo.type .<= "B")
# staInfo = filter(row -> row.city == "北京市" && row.type == "B" && row.province == "北京", staInfo)

include("/sh2/home/3clear/wuhj/paper_IETMsand/DFCI.jl")
lonM,latM=load("/sh2/home/3clear/wuhj/paper_IETMsand/datagrid.jld2","lonM","latM")
sta_lon = staInfo.lon
sta_lat = staInfo.lat

tStart_all=ZonedDateTime.(DateTime(2025,4,8):Hour(1):DateTime(2025,4,13),tz"Asia/Shanghai")
tStart_obs=ZonedDateTime.(DateTime(2025,4,11):Hour(1):DateTime(2025,4,13),tz"Asia/Shanghai")
matFun3d_v5(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_NationalVerified/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"
# matFun2d(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_2d_noAER/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"
matFun3d_v5_allsta(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_AllStaOri/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"

filepath_3d=matFun3d_v5.(tStart_all)
ietm3d_data=staDATA(filepath_3d,Times=tStart_all)
# filepath_2d=matFun2d.(tStart_all)
# ietm2d_data=staDATA(filepath_2d,Times=tStart_all)
filepath_3d_allsta=matFun3d_v5_allsta.(tStart_all)
ietm3d_allstadata=staDATA(filepath_3d_allsta,Times=tStart_all)
leadtimes=collect(1:72)
cities=["延安","太原","大同","忻州"]
function getstaID2plot(city)
    staID2plot=staInfo[isequal.(staInfo.city,city),"staID"]
    return staID2plot
end
#读数据
ietm3d=Array{Array{Float64},2}(undef,length(cities),length(leadtimes))
ietm3d_allsta=Array{Array{Float64},2}(undef,length(cities),length(leadtimes))

obs=Array{Array{Float64},1}(undef,length(cities))
for n_c in eachindex(cities)
    city=cities[n_c]
    obs[n_c]=vec(nanmean(read(staObs,Times=tStart_obs,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
    for n_l in eachindex(leadtimes)
        l=leadtimes[n_l]
        ietm3d[n_c,n_l]=vec(nanmean(read(ietm3d_data,Times=tStart_obs .- Hour(l),leadTimes=l,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
        ietm3d_allsta[n_c,n_l]=vec(nanmean(read(ietm3d_allstadata,Times=tStart_obs .- Hour(l),leadTimes=l,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))

    end
end
All_diff=Array{Float64}(undef,length(cities))
Nation_diff=Array{Float64}(undef,length(cities))
for n_c in eachindex(cities)
    All_max=indexin(maximum(ietm3d_allsta[n_c,2]),ietm3d_allsta[n_c,2])[1]
    Nation_max=indexin(maximum(ietm3d[n_c,2]),ietm3d[n_c,2])[1]
    obs_max=indexin(maximum(obs[n_c]),obs[n_c])[1]
    All_diff[n_c]=abs(All_max - obs_max)
    Nation_diff[n_c]=abs(Nation_max - obs_max)
end

#面图数据
Nation_data2plot=Array{Array{Float32,2},1}(undef,4)
for i in range(0,3)
    Time_grd=ZonedDateTime.(DateTime(2025,4,11,10),tz"Asia/Shanghai")
    Time_grd2Read=Time_grd+Hour(i*2)
    Time_grd_index=indexin([Time_grd2Read],tStart_obs)[1]
    grd_Nation=gradData("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_NationalVerified/DATA/RUF/$(Dates.format(astimezone(Time_grd2Read,tz"Asia/Shanghai"),"yyyymmddHH"))/ruf.ctl")
    grd_Nation_setTime=setTime!(grd_Nation, Time_grd2Read)
    Nation_data2plot[i+1] = read(grd_Nation_setTime, "PM10", 1)
end
Allsta_data2plot=Array{Array{Float32,2},1}(undef,4)
obs_data2plot=Array{Array{Float64},1}(undef,4)

for i in range(0,3)
    Time_grd=ZonedDateTime.(DateTime(2025,4,11,10),tz"Asia/Shanghai")
    Time_grd2Read=Time_grd+Hour(i*2)
    Time_grd_index=indexin([Time_grd2Read],tStart_obs)[1]
    grd_Allsta=gradData("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_AllStaOri/DATA/RUF/$(Dates.format(astimezone(Time_grd2Read,tz"Asia/Shanghai"),"yyyymmddHH"))/ruf.ctl")
    grd_Allsta_setTime=setTime!(grd_Allsta, Time_grd2Read)
    Allsta_data2plot[i+1] = read(grd_Allsta_setTime, "PM10", 1)
    obs_data2plot[i+1] = vec(read(staObs, Times=Time_grd2Read,varNames="PM10",
     staIDs=staInfo.staID))
end

cities_all=unique(staInfo.city)
Time_grd=ZonedDateTime.(DateTime(2025,4,11,10),tz"Asia/Shanghai")
Nation_data2plot_sta=Array{Float64}(undef,length(cities_all),4)
Allsta_data2plot_sta=Array{Float64}(undef,length(cities_all),4)
obs_all=Array{Float64,2}(undef,length(cities_all),4)
for n_c in eachindex(cities_all)
    city=cities_all[n_c]
    for t in 0:3
        n_t=t+1
        obs_all[n_c,n_t]=vec(nanmean(read(staObs,Times=Time_grd+Hour(t*2),varNames="PM10",staIDs=getstaID2plot(city)),dims=3))[1]
        Nation_data2plot_sta[n_c,n_t]=vec(nanmean(read(ietm3d_data,Times=Time_grd+Hour(t*2),leadTimes=0,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))[1]
        Allsta_data2plot_sta[n_c,n_t]=vec(nanmean(read(ietm3d_allstadata,Times=Time_grd+Hour(t*2),leadTimes=0,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))[1]

    end
end


Nation_TS = Array{Float64}(undef,4)
All_TS = Array{Float64}(undef,4)

for t in 0:3
    n_t=t+1
    obs_happen = findall(obs_all[:,n_t] .> 200)
    Nation_happen = findall(Nation_data2plot_sta[:,n_t] .> 200)
    All_happen = findall(Allsta_data2plot_sta[:,n_t] .> 200)

    # Nation_H = length(intersect(obs_happen, Nation_happen))
    # All_H = length(intersect(obs_happen, All_happen))
    # O = length(obs_happen)                 # H + M
    Nation_H = length(intersect(obs_happen, Nation_happen))
    Nation_M = length(setdiff(obs_happen, Nation_happen))
    Nation_F = length(setdiff(Nation_happen, obs_happen))
    Nation_den = Nation_H + Nation_M + Nation_F

    All_H = length(intersect(obs_happen, All_happen))
    All_M = length(setdiff(obs_happen, All_happen))
    All_F = length(setdiff(All_happen, obs_happen))
    All_den = All_H + All_M + All_F

    Nation_TS[n_t] = Nation_den == 0 ? NaN : Nation_H / Nation_den
    All_TS[n_t] = All_den == 0 ? NaN : All_H / All_den
end

(Nation_TS .- All_TS) ./ Nation_TS
end
# plt=get_plt()
# plt.clf()
# plt.plot(Nation_TS,color="r")
# plt.plot(All_TS,color="blue")
# plt.gcf()


#绘图
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

function add_SouthChinaAxes(axIn,data)
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
    ax_southChina.pcolormesh(lonM,latM,data,transform=PlaCar,cmap="afmhot_r",zorder=0)

    return ax_southChina
end


using PyCall
using PlotFunM
matplotlib = pyimport("matplotlib")
ticker=matplotlib.ticker
ccrs = pyimport("cartopy.crs")
shpreader = pyimport("cartopy.io.shapereader")
np=pyimport("numpy")
fm=matplotlib.font_manager
# includet("$(pathof(DASystemM))/../../system/"*"rapidUpdatedForecast_3d.jl")
PlaCar=ccrs.PlateCarree()
t2plot=DateTime.(tStart_obs)
t2plot_label=Dates.format.(tStart_obs,"mm-dd")
#获取四城市坐标
staInfo_4cities=staInfo[indexin(cities,city_labels),:]
cities_lon=staInfo_4cities.lon
cities_lat=staInfo_4cities.lat

function draw_scatter(ax,cities_lon,cities_lat)
    ax.scatter(cities_lon[1],cities_lat[1],marker="o",2,"r",transform=PlaCar,label="Yan'an",zorder=4)
    ax.scatter(cities_lon[2],cities_lat[2],marker="o",2,"b",transform=PlaCar,label="Taiyuan",zorder=4)
    ax.scatter(cities_lon[3],cities_lat[3],marker="o",2,"green",transform=PlaCar,label="Datong",zorder=4)
    ax.scatter(cities_lon[4],cities_lat[4],marker="o",2,"purple",transform=PlaCar,label="Xinzhou",zorder=4)

    # ax.text(cities_lon[1],cities_lat[1],"Yan'an",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax.text(cities_lon[2],cities_lat[2],"Taiyuan",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax.text(cities_lon[3],cities_lat[3]+0.3,"Datong",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax.text(cities_lon[4],cities_lat[4]-0.1,"Xinzhou",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax.text(cities_lon[1], cities_lat[1]-0.3, "Yan'an",
    # fontsize=4.5,
    # ha="center",             # 水平方向居中在点上
    # va="top",                # 文字顶部对齐点 → 显示在点下方
    # transform=PlaCar,
    # fontproperties=font_prop)
    # ax.text(cities_lon[2]+1, cities_lat[2]-0.3, "Taiyuan",
    # fontsize=4.5,
    # ha="left",             # 水平方向居中在点上
    # va="top",                # 文字顶部对齐点 → 显示在点下方
    # transform=PlaCar,
    # fontproperties=font_prop)
end
function ax_baseSet(ax)
    ax.set_aspect("auto")
    ax.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax.set_extent([90, 130, 30, 45],crs=PlaCar)
end
function add_grid_text(ax, x, y, content; 
                       rows=2, cols=2, 
                       col_space=0.2, 
                       row_space=0.1, 
                       transform=false, 
                       kwargs...)
    
    # 设置坐标系：ax.transAxes 表示相对子图位置，nothing 表示使用数据坐标
    tf = transform ? ax.transAxes : ax.transData
    
    # 遍历行与列
    for i in 0:(rows - 1)
        for j in 0:(cols - 1)
            # 计算当前索引
            idx = i * cols + j + 1
            
            # 防止 content 长度不足导致越界
            if idx <= length(content)
                # 计算当前格子的坐标
                # x 轴向右增加，y 轴向下减少（因为 y 通常代表高度）
                curr_x = x + j * col_space
                curr_y = y - i * row_space
                
                ax.text(curr_x, curr_y, content[idx]; 
                        transform=tf,
                        va="center", ha="left",
                        kwargs...)
            end
        end
    end
end
begin
    

    matplotlib = pyimport("matplotlib")
    configure_matplotlib!(matplotlib)
    plt=get_plt()
    # Tight bounding boxes add a small amount beyond the nominal canvas; this
    # width keeps the exported PDF within a 7-inch double-column layout.
    # Reduce only the original map-row heights by 15%, retaining automatic aspect.
    # Keep the lower panels and legends at their existing physical sizes/positions.
    figure_width = 6.55
    original_height = 5.4
    map_height_scale = 0.85
    map_block_height = original_height * (0.93 - 0.31) * map_height_scale
    figure_height = original_height * (1 - (0.93 - 0.31)) + map_block_height
    yscale = original_height / figure_height
    fig=plt.figure(figsize=(figure_width,figure_height),dpi=600);plt.clf()
    proj=ccrs.LambertConformal(central_longitude=105, central_latitude=83,standard_parallels=[25,45],globe=nothing)
    font_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_prop = fm.FontProperties(fname=font_path,size=TICK_SIZE)
    font_pathbd = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    font_propbd = fm.FontProperties(fname=font_pathbd,size=PANEL_TITLE_SIZE)
    font_axis = fm.FontProperties(fname=font_path,size=AXIS_LABEL_SIZE)
    font_annotation = fm.FontProperties(fname=font_path,size=ANNOTATION_SIZE)
    font_legend = fm.FontProperties(fname=font_path,size=FIGURE_LEGEND_SIZE)
    all_stations_color = "#2ca02c"
    title_content = ["RMSE","R"]

    function add_metric_table(axis, nation_rmse, nation_r, all_rmse, all_r)
        box = Dict("facecolor" => "white", "edgecolor" => "none", "alpha" => 0.72, "pad" => 0.25)
        for (x, text_value) in ((0.66, "RMSE"), (0.91, "R"))
            axis.text(x, 0.73, text_value, transform=axis.transAxes, va="top", ha="center",
                color="black", fontproperties=font_annotation, zorder=20, bbox=box)
        end
        for (y, rmse_value, r_value, color_value) in
                ((0.60, nation_rmse, nation_r, "#d62728"),
                 (0.47, all_rmse, all_r, all_stations_color))
            axis.text(0.66, y, string(rmse_value), transform=axis.transAxes, va="top", ha="center",
                color=color_value, fontproperties=font_annotation, zorder=20, bbox=box)
            axis.text(0.91, y, string(r_value), transform=axis.transAxes, va="top", ha="center",
                color=color_value, fontproperties=font_annotation, zorder=20, bbox=box)
        end
    end

    gs = fig.add_gridspec(nrows=3, ncols=4, left=0.075,right=0.98,
        bottom=0.31 * yscale,top=(0.31 * original_height + map_block_height) / figure_height,
        wspace=0.10,hspace=0.0)
    gs4 = fig.add_gridspec(nrows=1, ncols=4, left=0.075,right=0.98,
        bottom=0.045 * yscale,top=0.19 * yscale,wspace=0.28,hspace=0.0)
    
    ax1_1 = fig.add_subplot(gs[1], projection=PlaCar)
    ax1_2 = fig.add_subplot(gs[2], projection=PlaCar)
    ax1_3 = fig.add_subplot(gs[3], projection=PlaCar)
    ax1_4 = fig.add_subplot(gs[4], projection=PlaCar)
    ax2_1 = fig.add_subplot(gs[5], projection=PlaCar)
    ax2_2 = fig.add_subplot(gs[6], projection=PlaCar)
    ax2_3 = fig.add_subplot(gs[7], projection=PlaCar)
    ax2_4 = fig.add_subplot(gs[8], projection=PlaCar)
    ax3_1 = fig.add_subplot(gs[9], projection=PlaCar)
    ax3_2 = fig.add_subplot(gs[10], projection=PlaCar)
    ax3_3 = fig.add_subplot(gs[11], projection=PlaCar)
    ax3_4 = fig.add_subplot(gs[12], projection=PlaCar)
    ax3 = fig.add_subplot(gs4[1])
    ax4 = fig.add_subplot(gs4[2])
    ax5 = fig.add_subplot(gs4[3])
    ax6 = fig.add_subplot(gs4[4])


   
    # ax1_1=fig.add_axes([0.05, 0.7, 0.2, 0.22], projection=proj)  
    # # ax1_colorbar=fig.add_axes([0.07, 0.35, 0.35, 0.02])  
    # ax1_2=fig.add_axes([0.28, 0.7, 0.2, 0.22], projection=proj)  
    # # ax1_3=fig.add_axes([0.02, 0.45, 0.2, 0.22], projection=proj)  
    # # ax1_4=fig.add_axes([0.27, 0.45, 0.2, 0.22], projection=proj)  
    # ax1_3=fig.add_axes([0.51, 0.7, 0.2, 0.22], projection=proj)  
    # ax1_4=fig.add_axes([0.74, 0.7, 0.2, 0.22], projection=proj)  
    
    # ax2_1=fig.add_axes([0.05, 0.45, 0.2, 0.22], projection=proj)  
    # ax2_2=fig.add_axes([0.28, 0.45, 0.2, 0.22], projection=proj)  
    # # ax2_1=fig.add_axes([0.52, 0.7, 0.2, 0.22], projection=proj)  
    # # ax2_2=fig.add_axes([0.77, 0.7, 0.2, 0.22], projection=proj)  
    # ax2_3=fig.add_axes([0.51, 0.45, 0.2, 0.22], projection=proj)  
    # ax2_4=fig.add_axes([0.74, 0.45, 0.2, 0.22], projection=proj)  

    ax2_colorbar=fig.add_axes([0.60, 0.235 * yscale, 0.37, 0.012 * yscale])

    # ax3=fig.add_axes([0.04, 0.13, 0.2, 0.2])  
    # ax4=fig.add_axes([0.28, 0.13, 0.2, 0.2])  
    # ax5=fig.add_axes([0.52, 0.13, 0.2, 0.2])  
    # ax6=fig.add_axes([0.76, 0.13, 0.2, 0.2])  
    axex_3col = [ax3_1,ax3_2,ax3_3,ax3_4]
    for n_ax in eachindex(axex_3col) 
        ax = axex_3col[n_ax]
        ax_baseSet(ax)
        ax.scatter(sta_lon, sta_lat;
        c=obs_data2plot[n_ax],
        s=0.2,
        transform=PlaCar,
        cmap="afmhot_r",
        zorder=2,  
        vmin=0, vmax=1000
        )
        draw_scatter(ax,cities_lon,cities_lat)

    end
    ax1_1.set_aspect("auto")
    ax1_2.set_aspect("auto")
    ax1_3.set_aspect("auto")
    ax1_4.set_aspect("auto")
    ax2_1.set_aspect("auto")
    ax2_2.set_aspect("auto")
    ax2_3.set_aspect("auto")
    ax2_4.set_aspect("auto")
    ax3_1.set_aspect("auto")
    ax3_2.set_aspect("auto")
    ax3_3.set_aspect("auto")
    ax3_4.set_aspect("auto")


    ax1_1.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax1_2.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax1_3.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax1_4.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax2_1.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax2_2.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax2_3.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)
    ax2_4.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 3,linewidth=0.5)


    ax1_1.set_extent([90, 130, 30, 45],crs=PlaCar)
    ax1_2.set_extent([90, 130, 30, 45],crs=PlaCar)
    ax1_3.set_extent([90, 130, 30, 45],crs=PlaCar)
    ax1_4.set_extent([90, 130, 30, 45],crs=PlaCar)
    

    ax2_1.set_extent([90, 130, 30, 45],crs=PlaCar)
    ax2_2.set_extent([90, 130, 30, 45],crs=PlaCar)
    ax2_3.set_extent([90, 130, 30, 45],crs=PlaCar)
    ax2_4.set_extent([90, 130, 30, 45],crs=PlaCar)


    ax1_1.pcolormesh(lonM,latM,Nation_data2plot[1],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    ax1_2.pcolormesh(lonM,latM,Nation_data2plot[2],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    ax1_3.pcolormesh(lonM,latM,Nation_data2plot[3],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    ax1_4.pcolormesh(lonM,latM,Nation_data2plot[4],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    # axes_list = [ax1_1,ax1_2,ax1_3,ax1_4,ax2_1,ax2_2,ax2_3,ax2_4]
    # for i in 1:8
    #     ax  = axes_list[i]
    #     idx = i > 4 ? i - 4 : i
    #     ax.scatter(sta_lon, sta_lat;
    #         c=obs_data2plot[idx],
    #         s=0.2,
    #         transform=PlaCar,
    #         cmap="viridis",
    #         zorder=2,  
    #         vmin=0, vmax=2000
    #     )
        
    # end
    sc_ax2=ax2_1.pcolormesh(lonM,latM,Allsta_data2plot[1],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    cbar=fig.colorbar(sc_ax2,cax=ax2_colorbar,orientation="horizontal",extend="max")
    cbar.ax.set_title(raw"PM$_{10}$ (μg m$^{-3}$)",fontsize=AXIS_LABEL_SIZE,pad=2,loc="center",fontproperties=font_axis)
    cbar.ax.tick_params(axis="both",which="both",length=2.8,labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)

    ax2_2.pcolormesh(lonM,latM,Allsta_data2plot[2],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    ax2_3.pcolormesh(lonM,latM,Allsta_data2plot[3],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)
    ax2_4.pcolormesh(lonM,latM,Allsta_data2plot[4],transform=PlaCar,cmap="afmhot_r",zorder=1,vmin=0,vmax=1000)

    draw_scatter(ax1_1,cities_lon,cities_lat)
    draw_scatter(ax1_2,cities_lon,cities_lat)
    draw_scatter(ax1_3,cities_lon,cities_lat)
    draw_scatter(ax1_4,cities_lon,cities_lat)
    draw_scatter(ax2_1,cities_lon,cities_lat)
    draw_scatter(ax2_2,cities_lon,cities_lat)
    draw_scatter(ax2_3,cities_lon,cities_lat)
    draw_scatter(ax2_4,cities_lon,cities_lat)
    
    ax1_1_TS_content = ["TS:",round(Nation_TS[1],digits=2)]
    add_grid_text(ax1_1,0.7,0.9,ax1_1_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax1_2_TS_content = ["TS:",round(Nation_TS[2],digits=2)]
    add_grid_text(ax1_2,0.7,0.9,ax1_2_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax1_3_TS_content = ["TS:",round(Nation_TS[3],digits=2)]
    add_grid_text(ax1_3,0.7,0.9,ax1_3_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax1_4_TS_content = ["TS:",round(Nation_TS[4],digits=2)]
    add_grid_text(ax1_4,0.7,0.9,ax1_4_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax2_1_TS_content = ["TS:",round(All_TS[1],digits=2)]
    add_grid_text(ax2_1,0.7,0.9,ax2_1_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax2_2_TS_content = ["TS:",round(All_TS[2],digits=2)]
    add_grid_text(ax2_2,0.7,0.9,ax2_2_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax2_3_TS_content = ["TS:",round(All_TS[3],digits=2)]
    add_grid_text(ax2_3,0.7,0.9,ax2_3_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)
    ax2_4_TS_content = ["TS:",round(All_TS[4],digits=2)]
    add_grid_text(ax2_4,0.7,0.9,ax2_4_TS_content,rows=1,cols=2,col_space=0.1, row_space=0.3, color="navy", zorder=10,
    transform=true, fontsize=ANNOTATION_SIZE)

    # ax1.scatter(cities_lon,cities_lat,marker="*",3,"r",transform=PlaCar)
    # ax2.scatter(cities_lon,cities_lat,marker="*",3,"r",transform=PlaCar)
    # ax1.text(cities_lon[1],cities_lat[1],"Yan'an",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax2.text(cities_lon[1],cities_lat[1],"Yan'an",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax1.text(cities_lon[2],cities_lat[2],"Taiyuan",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax2.text(cities_lon[2],cities_lat[2],"Taiyuan",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax1.text(cities_lon[3],cities_lat[3],"Datong",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax2.text(cities_lon[3],cities_lat[3],"Datong",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax1.text(cities_lon[4],cities_lat[4],"Xinzhou",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    # ax2.text(cities_lon[4],cities_lat[4],"Xinzhou",fontsize=4.5,horizontalalignment="left",transform=PlaCar,fontproperties=font_prop)
    addGrids(ax1_1,false)
    addGrids(ax1_2,false)
    addGrids(ax1_3,false)
    addGrids(ax1_4,true)
    addGrids(ax2_1,false)
    addGrids(ax2_2,false)
    addGrids(ax2_3,false)
    addGrids(ax2_4,true)
    addGrids(ax3_1,false,true)
    addGrids(ax3_2,false,true)
    addGrids(ax3_3,false,true)
    addGrids(ax3_4,true,true)

    handles, labels = ax1_1.get_legend_handles_labels()

    # # 在整个 figure 上画图例，位置随意指定
    fig.legend(handles, labels,
               loc="lower left",       # 位置
               bbox_to_anchor=(0.08, 0.225 * yscale),  # 相对 figure 坐标
               prop=font_legend, ncol=4,
               frameon=true, handletextpad=0.3, columnspacing=0.6,
               borderpad=0.25)           # 是否带边框

    ax1_1.set_title("2025-04-11 10:00",fontproperties=font_propbd)
    ax1_2.set_title("2025-04-11 12:00",fontproperties=font_propbd)
    ax1_3.set_title("2025-04-11 14:00",fontproperties=font_propbd)
    ax1_4.set_title("2025-04-11 16:00",fontproperties=font_propbd)
    add_label_atleftup(ax1_1,"(a)")
    add_label_atleftup(ax1_2,"(b)")
    add_label_atleftup(ax1_3,"(c)")
    add_label_atleftup(ax1_4,"(d)")
    add_label_atleftup(ax2_1,"(e)")
    add_label_atleftup(ax2_2,"(f)")
    add_label_atleftup(ax2_3,"(g)")
    add_label_atleftup(ax2_4,"(h)")
    add_label_atleftup(ax3_1,"(i)")
    add_label_atleftup(ax3_2,"(j)")
    add_label_atleftup(ax3_3,"(k)")
    add_label_atleftup(ax3_4,"(l)")


    ax1_1.set_ylabel("Nation", fontproperties=font_propbd, labelpad=3)
    ax1_1.text(-0.1, 0.5, "Nation";
        transform=ax1_1.transAxes,
        ha="right", va="center",
        rotation=90,fontproperties=font_propbd)
    

    ax2_1.text(-0.1, 0.5, "All Stations";
    transform=ax2_1.transAxes,
    ha="right", va="center",
    rotation=90,fontproperties=font_propbd)
    # add_label_atleftup(ax1_1,"(a) Nation 4-11-18")
    ax3_1.text(-0.1, 0.5, "Observation";
    transform=ax3_1.transAxes,
    ha="right", va="center",
    rotation=90,fontproperties=font_propbd)
    # add_label_atleftup(ax1_1,"(a) Nation 4-11-18")
    # # add_label_atleftup(ax1_2,"(b) Nation 4-11-19")
    # # add_label_atleftup(ax1_3,"(c) Nation 4-11-20")
    # # add_label_atleftup(ax1_4,"(d) Nation 4-11-21")

    # # add_label_atleftup(ax2_1,"(e) AllStation 4-11-18")
    # # add_label_atleftup(ax2_2,"(f) AllStation 4-11-19")
    # # add_label_atleftup(ax2_3,"(g) AllStation 4-11-20")
    # # add_label_atleftup(ax2_4,"(h) AllStation 4-11-21")



    ax3.tick_params(axis="x", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    ax3.tick_params(axis="y", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    add_label_atleftup(ax3,"(m) Yan'an")
    ax3.bar(t2plot,obs[1],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    ax3.plot(t2plot,ietm3d[1,2],linewidth=DATA_LINEWIDTH,color=RUF_COLOR,linestyle="-",label="Nation")
    ax3.plot(t2plot,ietm3d_allsta[1,2],linewidth=DATA_LINEWIDTH,color=all_stations_color,linestyle="-",label="All Stations")
    ax3_nation_content = [round(Nation_city_rmse[1], digits=1), round(Nation_city_r[1], digits=2)]
    ax3_all_content = [round(All_city_rmse[1], digits=1), round(All_city_r[1], digits=2)]
    add_metric_table(ax3, ax3_nation_content[1], ax3_nation_content[2], ax3_all_content[1], ax3_all_content[2])
    
    ax3.set_xlim(DateTime(2025, 4, 11, 0, 0), DateTime(2025, 4, 13))
    # ax3.set_ylim(0,4000)
    ax3.set_xticks(t2plot[1:24:end]) 
    ax3.set_xticklabels(t2plot_label[1:24:end], fontsize=TICK_SIZE)
    xticks = ax3.get_xticklabels()
    # 最左边的
    xticks[1].set_ha("left")    # 左对齐 → 整体往右靠
    # 最右边的
    xticks[3].set_ha("right")  # 右对齐 → 整体往左靠
    ax3.set_ylabel(raw"$\mathrm{PM}_{10}$"*raw" (μg m$^{-3}$)",labelpad=0.8,fontproperties=font_axis)
    ax3.tick_params(axis="y", labelsize=TICK_SIZE,pad=1)

    ax4.tick_params(axis="x", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    ax4.tick_params(axis="y", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    add_label_atleftup(ax4,"(n) Taiyuan")
    ax4.bar(t2plot,obs[2],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    ax4.plot(t2plot,ietm3d[2,2],linewidth=DATA_LINEWIDTH,color=RUF_COLOR,linestyle="-",label="Nation")
    ax4.plot(t2plot,ietm3d_allsta[2,2],linewidth=DATA_LINEWIDTH,color=all_stations_color,linestyle="-",label="All Stations")
    
    ax4_nation_content = [round(Nation_city_rmse[2], digits=1), round(Nation_city_r[2], digits=2)]
    ax4_all_content = [round(All_city_rmse[2], digits=1), round(All_city_r[2], digits=2)]
    add_metric_table(ax4, ax4_nation_content[1], ax4_nation_content[2], ax4_all_content[1], ax4_all_content[2])

    ax4.set_xlim(DateTime(2025, 4, 11, 0, 0), DateTime(2025, 4, 13))
    # ax4.set_ylim(0,4000)
    ax4.set_xticks(t2plot[1:24:end]) 
    ax4.set_xticklabels(t2plot_label[1:24:end], fontsize=TICK_SIZE)
    xticks = ax4.get_xticklabels()
    # 最左边的
    xticks[1].set_ha("left")    # 左对齐 → 整体往右靠
    # 最右边的
    xticks[3].set_ha("right")  # 右对齐 → 整体往左靠
    ax4.tick_params(axis="y", labelsize=TICK_SIZE,pad=1)

    # # ax4.set_ylabel(raw"$\mathrm{PM}_{10}$"*raw"(μg m$^{-3}$)",size=5.5,labelpad=0.5,fontproperties=font_prop)
    # # ax4.legend(fontsize=4.5)



    ax5.tick_params(axis="x", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    ax5.tick_params(axis="y", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    add_label_atleftup(ax5,"(o) Datong")
    ax5.bar(t2plot,obs[3],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    ax5.plot(t2plot,ietm3d[3,2],linewidth=DATA_LINEWIDTH,color=RUF_COLOR,linestyle="-",label="Nation")
    ax5.plot(t2plot,ietm3d_allsta[3,2],linewidth=DATA_LINEWIDTH,color=all_stations_color,linestyle="-",label="All Stations")
    
    ax5_nation_content = [round(Nation_city_rmse[3], digits=1), round(Nation_city_r[3], digits=2)]
    ax5_all_content = [round(All_city_rmse[3], digits=1), round(All_city_r[3], digits=2)]
    add_metric_table(ax5, ax5_nation_content[1], ax5_nation_content[2], ax5_all_content[1], ax5_all_content[2])
    
    ax5.set_xlim(DateTime(2025, 4, 11, 0, 0), DateTime(2025, 4, 13))
    # ax5.set_ylim(0,4000)
    ax5.set_xticks(t2plot[1:24:end]) 
    ax5.set_xticklabels(t2plot_label[1:24:end], fontsize=TICK_SIZE)
    xticks = ax5.get_xticklabels()
    # 最左边的
    xticks[1].set_ha("left")    # 左对齐 → 整体往右靠
    # 最右边的
    xticks[3].set_ha("right")  # 右对齐 → 整体往左靠
    ax5.tick_params(axis="y", labelsize=TICK_SIZE,pad=1)

    # ax5.set_ylabel(raw"$\mathrm{PM}_{10}$"*raw"(μg m$^{-3}$)",size=5.5,labelpad=0.5,fontproperties=font_prop)
    # ax5.legend(fontsize=4.5)


    ax6.tick_params(axis="x", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    ax6.tick_params(axis="y", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH)
    add_label_atleftup(ax6,"(p) Xinzhou")
    ax6.bar(t2plot,obs[4],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    ax6.plot(t2plot,ietm3d[4,2],linewidth=DATA_LINEWIDTH,color=RUF_COLOR,linestyle="-",label="Nation")
    ax6.plot(t2plot,ietm3d_allsta[4,2],linewidth=DATA_LINEWIDTH,color=all_stations_color,linestyle="-",label="All Stations")
    
    ax6_nation_content = [round(Nation_city_rmse[4], digits=1), round(Nation_city_r[4], digits=2)]
    ax6_all_content = [round(All_city_rmse[4], digits=1), round(All_city_r[4], digits=2)]
    add_metric_table(ax6, ax6_nation_content[1], ax6_nation_content[2], ax6_all_content[1], ax6_all_content[2])

    ax6.set_xlim(DateTime(2025, 4, 11, 0, 0), DateTime(2025, 4, 13))
    # ax6.set_ylim(0,4000)

    ax6.set_xticks(t2plot[1:24:end]) 
    ax6.set_xticklabels(t2plot_label[1:24:end], fontsize=TICK_SIZE)
    xticks = ax6.get_xticklabels()
    # 最左边的
    xticks[1].set_ha("left")    # 左对齐 → 整体往右靠
    # 最右边的
    xticks[3].set_ha("right")  # 右对齐 → 整体往左靠
    ax6.tick_params(axis="y", labelsize=TICK_SIZE,pad=1)

    line_handles, line_labels = ax3.get_legend_handles_labels()
    line_legend = fig.legend(line_handles, line_labels, loc="lower center", bbox_to_anchor=(0.5, -0.050 * yscale),
        ncol=3, frameon=false, prop=font_legend, handlelength=2.0, columnspacing=1.0)
    for legend_text in line_legend.get_texts()
        legend_text.get_text() == "All Stations" && legend_text.set_color("black")
    end

    # ax6.set_ylabel(raw"$\mathrm{PM}_{10}$"*raw"(μg m$^{-3}$)",size=5.5,labelpad=0.5,fontproperties=font_prop)
    # ax6.legend(fontsize=4.5)

    plt.gcf()

end
output_dir = "/home/zsm/IETM/figure/PapreFigure/20260908"
mkpath(output_dir)
plt.savefig(joinpath(output_dir, "PaperFigure_CompareNationAndAllsta.png"),dpi=600,bbox_inches="tight")
plt.savefig(joinpath(output_dir, "PaperFigure_CompareNationAndAllsta.pdf"),bbox_inches="tight")
plt.close("all")
println("Saved PaperFigure_CompareNationAndAllsta.png and .pdf")
exit()

peak_con_index = indexin.(maximum.(obs),obs)

function calculate_half_statistic(sprit_index,obs,forecast)
    obsLeftHalf = obs[1:sprit_index]
    foreLeftHalf = forecast[1:sprit_index]
    obsRightHalf = obs[sprit_index:end]
    foreRighHalf = forecast[sprit_index:end]
    left_rmse = getStatistics(obsLeftHalf,foreLeftHalf)[1]
    left_r = getStatistics(obsLeftHalf,foreLeftHalf)[3]
    right_rmse = getStatistics(obsRightHalf,foreRighHalf)[1]
    right_r = getStatistics(obsRightHalf,foreRighHalf)[3]
    return left_rmse, left_r, right_rmse, right_r

end
Nation_only = calculate_half_statistic(peak_con_index[2][1],obs[2],ietm3d[2])
All_station = calculate_half_statistic(peak_con_index[2][1],obs[2],ietm3d_allsta[2])
for i in 1:4
    println("rmse：",getStatistics(obs[i],ietm3d_allsta[i])[1])
    println("相关系数：",getStatistics(obs[i],ietm3d_allsta[i])[3])
end
