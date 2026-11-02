using MyPlotM,Dates,TimeZones,NaNStatistics
using Statistics,PyCall
# using .DustForecastMetrics
include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle
matplotlib = pyimport("matplotlib")
configure_matplotlib!(matplotlib)
ticker=matplotlib.ticker
ccrs = pyimport("cartopy.crs")
shpreader = pyimport("cartopy.io.shapereader")
np=pyimport("numpy")
fm=matplotlib.font_manager
mcolors = pyimport("matplotlib.colors")
offsetbox = pyimport("matplotlib.offsetbox")
# MaxNLocator=pyimport("matplotlib.ticker")
# Legacy forecast-system and DFCI includes are not required for this figure.


using PyCall
offsetbox = pyimport("matplotlib.offsetbox")

function add_triple_stats(ax, x, y, m1_stats, m2_stats, obs_stats; 
                          c1=AI_COLOR, c2=RUF_COLOR, c_obs=NEUTRAL_COLOR,
                          h_sep=1, v_sep=2, fontsize=ANNOTATION_SIZE)
    # m1_stats/m2_stats/obs_stats = [mean, max]
    
    # 定义通用的创建 TextArea 的辅助闭包，减少重复代码并强制缩小 pad
    mk_text(txt, color, bold=false) = offsetbox.TextArea(txt, 
        textprops=Dict("fontsize"=>fontsize, "color"=>color, 
                       "fontproperties"=>font_prop,
                       "fontweight"=>(bold ? "bold" : "normal"),
                       "va"=>"baseline"))

    mean_row = offsetbox.HPacker(children=[
            mk_text("mean:", "black", true),
            mk_text("$(round(m1_stats[1], digits=1))", c1),
            mk_text("$(round(m2_stats[1], digits=1))", c2),
            mk_text("$(round(obs_stats[1], digits=1))", c_obs),
        ], align="baseline", pad=0, sep=h_sep)
    max_row = offsetbox.HPacker(children=[
            mk_text("max:", "black", true),
            mk_text("$(round(Int, m1_stats[2]))", c1),
            mk_text("$(round(Int, m2_stats[2]))", c2),
            mk_text("$(round(Int, obs_stats[2]))", c_obs),
        ], align="baseline", pad=0, sep=h_sep)
    box = offsetbox.VPacker(children=[mean_row, max_row], align="left", pad=0, sep=v_sep)

    # --- 放置到图上 ---
    anchored_box = offsetbox.AnnotationBbox(box, xy=(x, y),
                                            xycoords="axes fraction",
                                            box_alignment=(0, 1), 
                                            pad=0.2, # 缩小边框与文字的边距
                                            frameon=true, 
                                            bboxprops=Dict("facecolor"=>"white", "alpha"=>0.8, 
                                                           "edgecolor"=>"none", "boxstyle"=>"round,pad=0.1"))
    
    ax.add_artist(anchored_box)
end
function draw(ax,data,height,xtick,xlabel,bool_x_label)
    vmin=0
    vmax=1500
    ticker_spacing_y=2500
    row=np.arange(1,17)
    X = np.tile(row, (13, 1))
    sc=ax.pcolormesh(X,height,data,vmin=vmin,vmax=vmax,cmap="afmhot_r")
    # y_min = min(ax.get_ylim()[1], ax.get_ylim()[1])
    # y_max = max(ax.get_ylim()[2], ax.get_ylim()[2])
    # ax.set_ylim(y_min, y_max)
    # ax.set_xticklabels(xlabel,rotation=20)
    if bool_x_label
        ax.set_xticks(xtick)
        ax.set_xticklabels(xlabel,fontsize=TICK_SIZE,rotation=45,ha="right")
        ax.tick_params(axis="x", labelsize=TICK_SIZE,pad=2,width=AXIS_LINEWIDTH)
    else
        ax.tick_params(axis="x", which="both", bottom=false, top=false, labelbottom=false)
    end
    ax.yaxis.set_major_locator(ticker.MultipleLocator(ticker_spacing_y))
    ax.tick_params(axis="y", labelsize=TICK_SIZE,pad=2,width=AXIS_LINEWIDTH)
    # ax.set_xlim(xlabel_sub[1], xlabel_sub[end])

    # ax.set_title("$(Dates.format(time,"mm月dd日HH时"))")
    # ax.set_title("$(Dates.format(time,"HH,dd,mm,YY"))")
    return sc
end

function getcity_pos(m_city)

    # staInfo_selected=staInfo[indexin(m_city,staInfo.city),:]
    staInfo_selected=staInfo[staInfo.city .== m_city, :]
    lon=nanmean(staInfo_selected.lon)
    lat=nanmean(staInfo_selected.lat)
    pos=[lon,lat]
    return pos
    
end

# indexin(["海口"],staInfo.city)
function setpic_pos(row,column,row_pos,column_pos,rowspan, colspan)
    ax = plt.subplot2grid((row, column),(row_pos-1, column_pos-1), rowspan=rowspan, colspan=colspan)
    return ax
end

function addGrids(ax)
    gl = ax.gridlines(crs=PlaCar,draw_labels=true, dms=true, x_inline=false, y_inline=false, linewidth=0.5, color=[0.1,0.1,0.1], alpha=0.5,zorder=2,linestyle="--") # 添加经纬度网格，这里设定为黑色虚线半透明 
    gl.top_labels  = false;    gl.right_labels = false;   gl.rotate_labels = false
    gl.xlocator=ticker.FixedLocator([40,60,80,100,120,140,160,180])
    gl.xlabel_style = Dict(:size => TICK_SIZE, :fontname => "Arial")
    gl.ylabel_style = Dict(:size => TICK_SIZE, :fontname => "Arial")
end 

function add_SouthChinaAxes(axIn)
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
    ax_southChina.scatter(cities_lon,cities_lat,c=max_obs,3,cmap=ax1_cmap,zorder=10,vmin=0,vmax=maximum(max_obs) * 0.6,transform=PlaCar)
    return ax_southChina
end
if false
lonM,latM=load("/sh2/home/3clear/wuhj/paper_IETMsand/datagrid.jld2","lonM","latM")
# staInfo=readStaInfo("/sh2/home/3clear/wuhj/online_IETMsand/staInfo.csv")
staObs=staDATA("/sh2/home/3clear/wuhj/paper_IETMsand/ori202504.mat")
staInfo=staObs.staInfo
PlaCar=ccrs.PlateCarree()
# get meshDATA
time2d3d=ZonedDateTime.(DateTime(2025,4,11,18),tz"Asia/Shanghai")
# timeObs=ZonedDateTime.(DateTime(2025,4,11,22):Hour(1):DateTime(2025,4,14,19),tz"Asia/Shanghai")
timeObs=ZonedDateTime.(DateTime(2025,4,11):Hour(1):DateTime(2025,4,15),tz"Asia/Shanghai")
#leadtime3h(图上时间)
tStart_all=ZonedDateTime.(DateTime(2025,4,11,18):Hour(4):DateTime(2025,4,12,18),tz"Asia/Shanghai")
#起报时间
tStart_all_lead=ZonedDateTime.(DateTime(2025,4,11,15):Hour(4):DateTime(2025,4,12,15),tz"Asia/Shanghai")
tStart_allfore=ZonedDateTime.(DateTime(2025,4,8):Hour(1):DateTime(2025,4,15),tz"Asia/Shanghai")
# tStart_allfore=ZonedDateTime.(DateTime(2025,4,7):Hour(1):DateTime(2025,4,14),tz"Asia/Shanghai")

# tStart=DateTime(2023,4,11,8)
tStart=DateTime(2025,4,11,18)
# tStart_1=ZonedDateTime(DateTime(2025,4,11),tz"Asia/Shanghai")
#2023
# city=["呼和浩特","北京","廊坊","天津","沧州","德州","济南","泰安","济宁","枣庄","徐州","宿迁","淮安","南京","上海"]
#2025
cities = [
    "呼和浩特",
    "太原",
    "长治",
    "晋城",
    "洛阳",
    "南阳",
    "襄阳",
    "宜昌",
    "常德",
    "邵阳",
    "永州",
    "桂林",
    "梧州",
    "茂名",
    "海口"
]
var2plot="PM10"
# city=["张掖","银川","鄂尔多斯","朔州","大同","北京","天津","承德","朝阳","阜新","沈阳"]
matFun3d_v5(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_AllStaOri/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"
matFun2d(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_2d_noAER_ori/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"
filepath_2d=matFun2d.(tStart_allfore)
filepath_3d=matFun3d_v5.(tStart_allfore)
# ietm2d=staDATA("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_2d_noAER/DATA/RUF/$(Dates.format(tStart,"yyyymmddHH"))/staData.$(Dates.format(tStart,"yyyymmddHH")).mat")
# ietm3d=staDATA("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v5_2025/DATA/RUF/$(Dates.format(tStart,"yyyymmddHH"))/staData.$(Dates.format(tStart,"yyyymmddHH")).mat")
ietm2d=staDATA(filepath_2d,Times=tStart_allfore)
ietm3d=staDATA(filepath_3d,Times=tStart_allfore)
ctlFun(t)="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_AllStaOri/DATA/RUF/$(Dates.format(t,"yyyymmddHH"))/ruf.ctl"
filepath_ctl=ctlFun.(tStart_allfore)

#测试Dustforcastmetrics
# timeObs=ZonedDateTime.(DateTime(2025,4,10,23):Hour(1):DateTime(2025,4,14),tz"Asia/Shanghai")
# staids=staInfo[isequal.(staInfo.type,"B"),"staID"]

# # @time compute_asymmetric_pms(obs,ietm)
# # @time dtwe = compute_conditional_dtwe(obs, ietm)
# taeses_ietm3d=Array{Float64,2}(undef,length(staids),72)
# @time for n_sta in eachindex(staids)
#     obs=vec(read(staObs,Times=timeObs,varNames="PM10",staIDs=staids[n_sta]))
#     for i in range(1,72)
#         ietm=vec(nanmean(read(ietm3d,Times=timeObs-Hour(i),leadTimes=i,varNames="PM10",staIDs=staids[n_sta]),dims=3))
#         # push!(taeses_ietm3d,compute_taes(obs, ietm,tao=6.0,theta=4000.0,peak_threshold=200.0,dtwe_threshold=200.0))
#         taeses_ietm3d[n_sta,i]=compute_taes(obs, ietm,tao=6.0,theta=4000.0,peak_threshold=200.0,dtwe_threshold=200.0)
#     end
# end
# taeses_ietm2d=Array{Float64,2}(undef,length(staids),72)
# @time for n_sta in eachindex(staids)
#     obs=vec(read(staObs,Times=timeObs,varNames="PM10",staIDs=staids[n_sta]))
#     for i in range(1,72)
#         ietm=vec(nanmean(read(ietm2d,Times=timeObs-Hour(i),leadTimes=i,varNames="PM10",staIDs=staids),dims=3))
#         # push!(taeses_ietm2d,compute_taes(obs, ietm,tao=6.0,theta=4000.0,peak_threshold=200.0,dtwe_threshold=200.0))
#         taeses_ietm2d[n_sta,i]=compute_taes(obs, ietm,tao=6.0,theta=4000.0,peak_threshold=200.0,dtwe_threshold=200.0)
#     end
# end
# plt=get_plt()
# plt.gca()
# fig=plt.figure(figsize=(12,6),dpi=300);plt.clf()
# # plt.plot(obs,label="obs")
# # plt.plot(ietm,label="ietm")
# plt.plot(vec(nanmean(taeses_ietm3d,dims=1)),label="3d")
# plt.plot(vec(nanmean(taeses_ietm2d,dims=1)),label="2d")
# plt.legend()
# plt.gcf()
#################################################################################


function get(start_time,time,city)
    grd1=gradData("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v6_2025_AllStaOri/DATA/RUF/$(Dates.format(start_time,"yyyymmddHH"))/ruf.ctl")
    meshData=read(setTime!(grd1,time),"PM10")
    # prepare staInfo to plot
    staInfo_selected=staInfo[indexin(city,staInfo.city),:]
    # extract data to plot
    oMapper=obsMapper(lonM,latM,staInfo_selected.lon,staInfo_selected.lat)
    data2plot=vcat([map2sta(oMapper,meshData[:,:,z])' for z in axes(meshData,3)]...)
    return data2plot
end
function get_truncated_cmap(base_cmap_name, min_val=0.0, max_val=1.0, n=256)
    base_cmap = plt.get_cmap(base_cmap_name)
    # 在 min_val 到 max_val 之间重采样
    colors_resampled = base_cmap(np.linspace(min_val, max_val, n))
    # 创建新的颜色条对象
    new_cmap = mcolors.LinearSegmentedColormap.from_list("trunc_afmhot", colors_resampled)
    return new_cmap
end
#测试
# grd1=gradData("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d_v5/DATA/RUF/$(Dates.format(tStart,"yyyymmddHH"))/ruf.ctl")
# meshData=read(setTime!(grd1,tStart_all[2]),"PM10")
# # prepare staInfo to plot
# staInfo_selected=staInfo[indexin(city,staInfo.city),:]
# oMapper=obsMapper(lonM,latM,staInfo_selected.lon,staInfo_selected.lat)

function getstaID2plot(city)
    staID2plot=staInfo[isequal.(staInfo.city,city),"staID"]
    return staID2plot
end
function get_ietm(ietm,tStart,leadtime,varname,city)
    ietm2plot=vec(nanmean(read(ietm;Times=tStart,leadTimes=leadtime,varNames=varname,staIDs=getstaID2plot(city)),dims=3))
    # ietm2plot=vec(read(ietm;Times=tStart,leadTimes=leadtime,varNames=varname,staIDs=getstaID2plot(city)[1]))
    
    return ietm2plot
end
#error的读取函数
function get(staFore,tStart_time_nh,varName,city)
    IDs=staInfo.staID[in.(staInfo.city,[[city]])]
    data=read(staFore,Times=tStart_time_nh,varNames=varName,staIDs=IDs)
    return data
end
function get_error(staObs,staFore,tStart_all_obs,city,tStart_time_nh,varName,leadtime)
    #特定城市
    obs=vec(read(staObs;Times=tStart_all_obs,varNames=var2plot,staIDs=getstaID2plot(city)))
    obs=Float64.(obs)
    IETM=get(staFore,tStart_time_nh,varName,city)
    IETM=vec(IETM[:,:,:,leadtime])
    IETM=Float64.(IETM)
    error=getStatistics(obs,IETM)
    error_arr=collect(values(error))
    return error_arr

end
# function get_error(staObs,staFore,tStart_all_obs,city,tStart_time_nh,varName,leadtime,staid)
#     #特定城市
#     obs=vec(read(staObs;Times=tStart_all_obs,varNames=var2plot,staIDs=staid))
#     obs=Float64.(obs)
#     IETM=get(staFore,tStart_time_nh,varName,city)
#     IETM=vec(IETM[:,:,:,leadtime])
#     IETM=Float64.(IETM)
#     error=getStatistics(obs,IETM)
#     error_arr=collect(values(error))
#     return error_arr

# end
# IETM=get(ietm3d,timeObs.-Hour(3),"PM10",cities_error[1])
# IETM=vec(IETM[:,:,:,3])


# for i in eachindex(ietm3d.leadTimes)
#     ietm3d2plot_bj[i]=get_ietm(ietm3d,time2d3d,ietm3d.leadTimes[i],var2plot,"北京")[1]
#     ietm2d2plot_bj[i]=get_ietm(ietm2d,time2d3d,ietm3d.leadTimes[i],var2plot,"北京")[1]
#     ietm3d2plot_nj[i]=get_ietm(ietm3d,time2d3d,ietm3d.leadTimes[i],var2plot,"南京")[1]
#     ietm2d2plot_nj[i]=get_ietm(ietm2d,time2d3d,ietm3d.leadTimes[i],var2plot,"南京")[1]
    
    
    
#     # push!(ietm3d2plot_bj,get_ietm(ietm3d,time2d3d,i,var2plot,"北京")[1])
#     # push!(ietm2d2plot_bj,get_ietm(ietm2d,time2d3d,i,var2plot,"北京")[1])
#     # push!(ietm3d2plot_nj,get_ietm(ietm3d,time2d3d,i,var2plot,"南京")[1])
#     # push!(ietm2d2plot_nj,get_ietm(ietm2d,time2d3d,i,var2plot,"南京")[1])
# end
# ietm3d2plot_bj=get_ietm(ietm3d,timeObs.- Hour(6),6,var2plot,"北京")
# ietm2d2plot_bj=get_ietm(ietm2d,timeObs.- Hour(6),6,var2plot,"北京")
# ietm3d2plot_nj=get_ietm(ietm3d,timeObs.- Hour(6),6,var2plot,"南京")
# ietm2d2plot_nj=get_ietm(ietm2d,timeObs.- Hour(6),6,var2plot,"南京")

# cities_error=["太原","茂名"]
cities_error=["太原","洛阳","常德","海口"]
leadtime=collect(1:72)
#15座城市平均RMSE
# error3d=Array{Float64,2}(undef,length(city))
# error2d=Array{Float64,2}(undef,length(city))
# for n_c in eachindex(city)
#     RUF3d = read(ietm3d,Times=t2plot)
#     error3d[n_l,n_c]=get_error(staObs,ietm3d,timeObs,city[n_c],timeObs.- Hour(i),var2plot,i)[1]
#     error2d[n_l,n_c]=get_error(staObs,ietm2d,timeObs,city[n_c],timeObs.- Hour(i),var2plot,i)[1]


# end
# error3d2plot=vec(nanmean(error3d,dims=2))
# error2d2plot=vec(nanmean(error2d,dims=2))
#全国所有国控站平均RMSE
function get_byid(staFore,tStart_time_nh,varName,staid)
    data=read(staFore,Times=tStart_time_nh,varNames=varName,staIDs=staid)
    return data
end
function get_errorbyid(staObs,staFore,tStart_all_obs,staid,tStart_time_nh,varName,leadtime)
    #特定城市
    obs=vec(read(staObs;Times=tStart_all_obs,varNames=varName,staIDs=staid))
    obs=Float64.(obs)
    IETM=get_byid(staFore,tStart_time_nh,varName,staid)
    IETM=vec(IETM[:,:,:,leadtime])
    IETM=Float64.(IETM)
    error=getStatistics(obs,IETM)
    error_arr=collect(values(error))
    return error_arr

end


# staids=staInfo[isequal.(staInfo.type,"B"),"staID"]
# error3d=Array{Float64,2}(undef,length(ietm3d.leadTimes)-1,length(staids))
# error2d=Array{Float64,2}(undef,length(ietm3d.leadTimes)-1,length(staids))
# @time for n_id in eachindex(staids)
#     staid=staids[n_id]
#     for n_l in eachindex(leadtime)
#         i=leadtime[n_l]
#         # println(n_id)
#         # println(n_l)
#         error3d[n_l,n_id]=get_errorbyid(staObs,ietm3d,timeObs,staid,timeObs.- Hour(i),var2plot,i)[1]
#         error2d[n_l,n_id]=get_errorbyid(staObs,ietm2d,timeObs,staid,timeObs.- Hour(i),var2plot,i)[1]
#     end
#     println("进度:"*"$(n_id/length(staids))")
# end


# rmse_3d=Array{Float64}(undef,length(cities))
# rmse_2d=Array{Float64}(undef,length(cities))
# r_3d=Array{Float64}(undef,length(cities))
# r_2d=Array{Float64}(undef,length(cities))

# for n_c in eachindex(cities)
#     city=cities[n_c]
#     obs_curr=vec(nanmean(read(staObs,Times=t_obs,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
#     ruf3d=vec(nanmean(read(ietm3d,Times=tStart_fore,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
#     ruf2d=vec(nanmean(read(ietm2d,Times=tStart_fore,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
#     mean_3d[n_c]=nanmean(ruf3d)
#     mean_2d[n_c]=nanmean(ruf2d)
#     mean_obs[n_c]=nanmean(obs_curr)
    
# end
#新参数计算
tStart_fore=ZonedDateTime(DateTime(2025,4,11,15),tz"Asia/Shanghai")
t_obs=ZonedDateTime.(DateTime(2025,4,11,15):Hour(1):DateTime(2025,4,11,15) + Hour(72),tz"Asia/Shanghai")
ietm3d_Data=Array{Array{Float64}}(undef,length(cities_error))
ietm2d_Data=Array{Array{Float64}}(undef,length(cities_error))
obs=Array{Array{Float64},1}(undef,length(cities_error))
mean_3d=Array{Float64}(undef,length(cities_error))
mean_2d=Array{Float64}(undef,length(cities_error))
mean_obs=Array{Float64}(undef,length(cities_error))
max_3d=Array{Float64}(undef,length(cities_error))
max_2d=Array{Float64}(undef,length(cities_error))
for n_c in eachindex(cities_error)
    city=cities_error[n_c]
    obs[n_c]=vec(nanmean(read(staObs,Times=t_obs,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
    ietm3d_Data[n_c]=vec(nanmean(read(ietm3d,Times=tStart_fore,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
    ietm2d_Data[n_c]=vec(nanmean(read(ietm2d,Times=tStart_fore,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
    mean_3d[n_c]=nanmean(ietm3d_Data[n_c])
    mean_2d[n_c]=nanmean(ietm2d_Data[n_c])
    mean_obs[n_c]=nanmean(obs[n_c])
    max_3d[n_c]=nanmaximum(ietm3d_Data[n_c])
    max_2d[n_c]=nanmaximum(ietm2d_Data[n_c])
end
cities_all=unique(staObs.staInfo.city)
max_obs = Array{Float64}(undef,length(cities_all))
cities_lon=Array{Float64}(undef,length(cities_all))
cities_lat=Array{Float64}(undef,length(cities_all))
for n_c in eachindex(cities_all)
    city=cities_all[n_c]
    obs_curr=vec(nanmean(read(staObs,Times=t_obs,varNames="PM10",staIDs=getstaID2plot(city)),dims=3))
    max_obs[n_c]=maximum(obs_curr)
    cities_pos=getcity_pos(city)
    cities_lon[n_c]=cities_pos[1]
    cities_lat[n_c]=cities_pos[2]
end
# cities_lon=cities_pos[:][1]



# ================================
#多次预报
# ietm3d2plot_bj=get_ietm(ietm3d,timeObs.- Hour(12),12,var2plot,cities_error[1])
# ietm2d2plot_bj=get_ietm(ietm2d,timeObs.- Hour(12),12,var2plot,cities_error[1])
# ietm3d2plot_nj=get_ietm(ietm3d,timeObs.- Hour(12),12,var2plot,cities_error[2])
# ietm2d2plot_nj=get_ietm(ietm2d,timeObs.- Hour(12),12,var2plot,cities_error[2])
# obs2plot_bj=vec(nanmean(read(staObs;Times=timeObs,varNames=var2plot,staIDs=getstaID2plot(cities_error[1])),dims=3))
# obs2plot_nj=vec(nanmean(read(staObs;Times=timeObs,varNames=var2plot,staIDs=getstaID2plot(cities_error[2])),dims=3))


include("/sh2/home/3clear/wuhj/paper_IETMsand/change_point_detect.jl")
leadtimes=collect(1:72)
derivative_threshold = 0.1
concentration_threshold = 200.0
obs_point,obs_normal_der=detect_starting_points(obs2plot_nj,derivative_threshold,concentration_threshold);
ietm_point,ietm_normal_der=detect_starting_points(ietm3d2plot_nj,derivative_threshold,concentration_threshold)
bool_ietm=judge_forcast(obs2plot_bj,ietm2d2plot_bj,Int.(obs_point),Int.(ietm_point),6);

indexin([maximum(ietm3d2plot_nj)],ietm3d2plot_nj)
indexin([maximum(ietm2d2plot_nj)],ietm2d2plot_nj)
ietm3d_max=maximum(ietm3d2plot_nj)
ietm2d_max=maximum(ietm2d2plot_nj)
obs_max=maximum(obs2plot_nj)
# #计算偏差
bias=abs(ietm2d_max - obs_max) / obs_max
#计算DTW最大改善
maximum(abs.(dtw_3d[2,:] .- dtw_2d[2,:]) ./ dtw_2d[2,:])
# ==============================


ax5data=get(tStart_all_lead[2],tStart_all[2],cities)
ax6data=get(tStart_all_lead[3],tStart_all[3],cities)
ax7data=get(tStart_all_lead[4],tStart_all[4],cities)
ax8data=get(tStart_all_lead[5],tStart_all[5],cities)
ax9data=get(tStart_all_lead[6],tStart_all[6],cities)
ax10data=get(tStart_all_lead[7],tStart_all[7],cities)

row=np.arange(1,17)
X1 = np.tile(row, (13, 1))


####提取高度信息#######

grd=gradData("/sh2/home/3clear/post/Mdata/china_3_27km/naqp/grid/grid.d2.ctl")
staInfo_selected=staInfo[indexin(cities,staInfo.city),:]
oMapper=obsMapper(lonM,latM,staInfo_selected.lon,staInfo_selected.lat)
height=Array{Any,2}(undef,13,length(cities)+1)
fill!(height,1000)
# a=read(grd,"heiz",1)
# height1=map2sta(oMapper,a)
for h in range(1,13)
    height_mesh=read(grd,"heiz",h)
    height_sta=map2sta(oMapper,height_mesh)
    height1floor=map2sta(oMapper,read(grd,"heiz",1))
    for i in eachindex(height_sta)
        height[h,i]=height_sta[i] - height1floor[i]
    end
end
height=reverse(height,dims=1)
# row=np.arange(1,17)
# X = np.tile(row, (13, 1))
# ax7.pcolormesh(X,height,ax10data,vmin=0,vmax=1000,cmap="afmhot_r")

#绘图数据准备
stacon=Dict()
# z=collect(1:12)
# x=[1,2,3,4]
# xtick=collect(0.5:2:15.5)
xtick=collect(1.5:2:15.5)
xcity=cities
allpos=getcity_pos.(cities)
allpos = collect.(zip(allpos...))
# city2plot=city=["呼和浩特","北京","沧州","济南","济宁"
# ,"徐州","淮安","上海"]
city2plot=cities[1:2:15]
citypos=getcity_pos.(city2plot)
citypos = collect.(zip(citypos...))
ty_pos=getcity_pos("太原")
mm_pos=getcity_pos("茂名")

t2plot=DateTime.(t_obs)
t2plot_label=Dates.format.(t_obs,"mm-dd")
end
include(joinpath(@__DIR__, "PaperFigure_VerticaltransportImproves_data.jl"))

function add_label_atleftup(ax, label, position; xpos=0.02, ypos=0.95, fontsize=PANEL_TITLE_SIZE)
    ax.text(xpos, ypos, label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            va="top",
            ha=position,
            fontproperties=font_propbd)
end
begin
    #绘图
    #第二版12x6 标题12 刻度8
    plt=get_plt()
    fig=plt.figure(figsize=(7.15,4.20),dpi=600);plt.clf()
    xlabel=["Hohhot", "Taiyuan", "Changzhi", "Jincheng", "Luoyang", "Nanyang", 
    "Xiangyang", "Yichang", "Changde", "Shaoyang", "Yongzhou", 
    "Guilin", "Wuzhou", "Maoming", "Haikou"]
    city_en=xlabel[1:2:15]

    
    xlabel_sub=xlabel[1:2:end]
    proj=ccrs.LambertConformal(central_longitude=105, central_latitude=83,standard_parallels=[25,45],globe=nothing)
    # ax3 = fig.add_axes([0.08, 0.21, 0.16, 0.16])  
    # ax3_1 = fig.add_axes([0.28, 0.21, 0.16, 0.16])
    # ax4 = fig.add_axes([0.08, 0.05, 0.16, 0.16]) 
    # ax4_1 = fig.add_axes([0.28, 0.05, 0.16, 0.16])
    # cbar_ax = fig.add_axes([0.11, 0.41, 0.3, 0.012]) 
    MaxNLocator=ticker.MaxNLocator
   # --- 1. 定义整体偏移量 (Offset) ---
    # 增加整体的左边距和底边距，以容纳标签，同时让图形在视觉上居中
    # dx = 0.03  # 水平向右平移
    dx = 0.0
    dy = 0.0

    # --- 2. 左侧子图 (ax1, ax3, ax4 等) ---
    # 地图主图
    ax1 = fig.add_axes([0.075 + dx, 0.585 + dy, 0.43, 0.39], projection=proj)
    # 地图对应的水平 Colorbar
    cbar_ax = fig.add_axes([0.125 + dx, 0.475 + dy, 0.33, 0.012])

    # 四个城市时间序列子图
    ax3   = fig.add_axes([0.09 + dx, 0.28 + dy, 0.185, 0.16])
    ax3_1 = fig.add_axes([0.305 + dx, 0.28 + dy, 0.185, 0.16])
    ax4   = fig.add_axes([0.09 + dx, 0.10 + dy, 0.185, 0.16])
    ax4_1 = fig.add_axes([0.305 + dx, 0.10 + dy, 0.185, 0.16])

    # --- 3. 右侧垂直剖面子图 (ax5 - ax10) ---
    # 整体向上平移以对齐 ax1 的顶端
    total_top = 0.945
    total_bottom = 0.095
    n_subs = 6
    h = (total_top - total_bottom) / n_subs
    left = 0.56
    width = 0.32

    ax5  = fig.add_axes([left, total_top - 1*h, width, h]) 
    ax6  = fig.add_axes([left, total_top - 2*h, width, h]) 
    ax7  = fig.add_axes([left, total_top - 3*h, width, h]) 
    ax8  = fig.add_axes([left, total_top - 4*h, width, h]) 
    ax9  = fig.add_axes([left, total_top - 5*h, width, h]) 
    ax10 = fig.add_axes([left, total_bottom,    width, h])
    # 右侧剖面图对应的垂直 Colorbar
    cax = fig.add_axes([0.92, 0.115, 0.012, 0.79])

    ax1.set_extent([72.3, 136.3, 16.9, 54.2],crs=PlaCar)
    cn_boundary=ax1.add_geometries(shpreader.Reader("/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    font_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_prop = fm.FontProperties(fname=font_path,size=TICK_SIZE)
    font_prop_ax1 = fm.FontProperties(fname=font_path,weight="heavy",size=ANNOTATION_SIZE)
    font_pathbd = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    font_propbd = fm.FontProperties(fname=font_pathbd,size=PANEL_TITLE_SIZE)
    font_axis = fm.FontProperties(fname=font_path,size=AXIS_LABEL_SIZE)
    font_legend = fm.FontProperties(fname=font_path,size=PANEL_LEGEND_SIZE)
    font_legend_small = fm.FontProperties(fname=font_path,size=PANEL_LEGEND_SIZE)
    fig.text(0.495, (total_top + total_bottom) / 2, "Height (m)", rotation=90,
        va="center", ha="center", fontproperties=font_axis)
    xticks=[70,80,90,100,110,120,130,140]
    ax1_cmap = get_truncated_cmap("afmhot_r",0.3,1)
    addGrids(ax1)
    add_SouthChinaAxes(ax1)


    sc_ax1=ax1.scatter(cities_lon,cities_lat,c=max_obs,3,cmap=ax1_cmap,zorder=4,vmin=0,vmax=nanmaximum(max_obs) * 0.6,transform=PlaCar)
    ax1.scatter(citypos[1],citypos[2],c=RUF_COLOR,6,zorder=200,transform=PlaCar)
    ax1.scatter(ty_pos[1],ty_pos[2],6,RUF_COLOR,zorder=200,transform=PlaCar)
    ax1.scatter(mm_pos[1],mm_pos[2],6,RUF_COLOR,zorder=200,transform=PlaCar)
    ax1.plot(allpos[1],allpos[2],zorder=150,transform=PlaCar)
    add_label_atleftup(ax1,"(a) Cities for Dust Transport Analysis","left")
    cb = fig.colorbar(sc_ax1, cax=cbar_ax, extend="max", orientation="horizontal")
    cb.set_label(raw"Peak PM$_{10}$ (μg m$^{-3}$)",
        size=ANNOTATION_SIZE, labelpad=2, fontproperties=font_prop)
    cb.ax.xaxis.set_label_position("top")
    cb.ax.tick_params(labelsize=TICK_SIZE,width=AXIS_LINEWIDTH,length=2.8)
    white_bbox = Dict(
    "facecolor" => "white", 
    "edgecolor" => "none", 
    "alpha"     => 0.7,    # 0.7 的透明度可以让底下的地图纹理隐约可见，更高级
    "pad"       => 0.5,    # 边框内衬大小
    "boxstyle"  => "round,pad=0.2" # 圆角矩形
    )
    for i in range(1,length(citypos[1]))
        if i%2==0
        ax1.text(citypos[1][i],citypos[2][i],city_en[i],fontsize=ANNOTATION_SIZE,horizontalalignment="left",zorder=100,transform=PlaCar,fontproperties=font_prop_ax1,bbox=white_bbox)
        else
        ax1.text(citypos[1][i],citypos[2][i],city_en[i],fontsize=ANNOTATION_SIZE,horizontalalignment="right",zorder=100,transform=PlaCar,fontproperties=font_prop_ax1,bbox=white_bbox)
        end
    end
    sc=draw(ax5,ax5data,height,xtick,xlabel_sub,false)
    cbar=fig.colorbar(sc,cax,extend="max")
    font_prop.set_size(TICK_SIZE)
    cbar.ax.set_title("   " * raw"μg m$^{-3}$",fontsize=AXIS_LABEL_SIZE,loc="center",fontproperties=font_axis)
    cbar.ax.tick_params(axis="both",which="both",length=2.8,labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)

    ticker_spacing_x=10
    font_prop.set_size(TICK_SIZE)

    plot_max=2500
    ax3.plot(t2plot,ietm3d_Data[1],RUF_COLOR,linewidth=DATA_LINEWIDTH,linestyle="-",label="RUF")
    ax3.plot(t2plot,ietm2d_Data[1],AI_COLOR,linewidth=DATA_LINEWIDTH,linestyle="--",label=raw"RUF$^{2d}$")
    ax3.bar(t2plot,obs[1],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    add_triple_stats(ax3, 0.05, 0.65,
            [round(mean_2d[1],digits=2), round(max_2d[1],digits=2)] , # RUF2d: mean, max
            [round(mean_3d[1],digits=2), round(max_3d[1],digits=2)], # RUF3d: mean, max
            [round(mean_obs[1],digits=2), round(max_obs_selected[1],digits=2)], # Obs: mean, max
            c1=AI_COLOR, c2=RUF_COLOR, c_obs=NEUTRAL_COLOR,
            h_sep=2, v_sep=1, fontsize=ANNOTATION_SIZE)
    ax3.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    add_label_atleftup(ax3,"(b) Taiyuan","left")
    font_prop.set_size(TICK_SIZE)
    ax3.margins(x=0)
    ax3.set_xlim(first(t2plot), last(t2plot))
    ax3.set_ylim(0,plot_max)
    ax3.set_xticks(t2plot[10:24:end]) 
    ax3.tick_params(axis="x", labelbottom=false)


    ax3_1.plot(t2plot,ietm3d_Data[2],RUF_COLOR,linewidth=DATA_LINEWIDTH,linestyle="-",label="RUF")
    ax3_1.plot(t2plot,ietm2d_Data[2],AI_COLOR,linewidth=DATA_LINEWIDTH,linestyle="--",label=raw"RUF$^{2d}$")
    ax3_1.bar(t2plot,obs[2],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    # add_dual_model_stats(ax3_1, 0.22, 0.5, [round(mean_2d[2],digits=2), round(max_2d[2],digits=2)],
    # [round(mean_3d[2],digits=2), round(max_3d[2],digits=2)], c1="#2ca02c", c2="#d62728", fontsize=4.5)
    add_triple_stats(ax3_1, 0.05, 0.65,
            [round(mean_2d[2],digits=2), round(max_2d[2],digits=2)] , # RUF2d: mean, max
            [round(mean_3d[2],digits=2), round(max_3d[2],digits=2)], # RUF3d: mean, max
            [round(mean_obs[2],digits=2), round(max_obs_selected[2],digits=2)], # Obs: mean, max
            c1=AI_COLOR, c2=RUF_COLOR, c_obs=NEUTRAL_COLOR,
            h_sep=2, v_sep=1, fontsize=ANNOTATION_SIZE)
    ax3_1.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    add_label_atleftup(ax3_1,"(c) Luoyang","left")
    ax3_1.margins(x=0)
    ax3_1.set_xlim(first(t2plot), last(t2plot))
    ax3_1.set_ylim(0,plot_max)
    ax3_1.set_xticks(t2plot[10:24:end]) 
    ax3_1.tick_params(axis="x", labelbottom=false)
    ax3_1.tick_params(axis="y", labelleft=false)





    ax4.plot(t2plot,ietm3d_Data[3],RUF_COLOR,linewidth=DATA_LINEWIDTH,linestyle="-",label="RUF")
    ax4.plot(t2plot,ietm2d_Data[3],AI_COLOR,linewidth=DATA_LINEWIDTH,linestyle="--",label=raw"RUF$^{2d}$")
    ax4.bar(t2plot,obs[3],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    add_triple_stats(ax4, 0.05, 0.65,
                [round(mean_2d[3],digits=2), round(max_2d[3],digits=2)], # RUF2d: mean, max
                [round(mean_3d[3],digits=2), round(max_3d[3],digits=2)], # RUF3d: mean, max
                [round(mean_obs[3],digits=2), round(max_obs_selected[3],digits=2)], # Obs: mean, max
                 c1=AI_COLOR, c2=RUF_COLOR, c_obs=NEUTRAL_COLOR,
                 h_sep=2, v_sep=1, fontsize=ANNOTATION_SIZE)
    ax4.margins(x=0)
    ax4.set_xlim(first(t2plot), last(t2plot))
    ax4.set_ylim(0,plot_max)
    ax4.set_xticks(t2plot[10:24:end]) 
    ax4.set_xticklabels(t2plot_label[10:24:end], fontsize=TICK_SIZE)
    ax4.tick_params(axis="x", pad=2)
    ax4.tick_params(axis="x", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    ax4.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    add_label_atleftup(ax4,"(d) Changde","left")


    
    ax4_1.plot(t2plot,ietm3d_Data[4],RUF_COLOR,linewidth=DATA_LINEWIDTH,linestyle="-",label="RUF")
    ax4_1.plot(t2plot,ietm2d_Data[4],AI_COLOR,linewidth=DATA_LINEWIDTH,linestyle="--",label=raw"RUF$^{2d}$")
    ax4_1.bar(t2plot,obs[4],width=0.02,facecolor=NEUTRAL_COLOR,label="obs")
    add_triple_stats(ax4_1, 0.05, 0.65,
            [round(mean_2d[4],digits=2), round(max_2d[4],digits=2)] , # RUF2d: mean, max
            [round(mean_3d[4],digits=2), round(max_3d[4],digits=2)], # RUF3d: mean, max
            [round(mean_obs[4],digits=2), round(max_obs_selected[4],digits=2)], # Obs: mean, max
            c1=AI_COLOR, c2=RUF_COLOR, c_obs=NEUTRAL_COLOR,
            h_sep=2, v_sep=1, fontsize=ANNOTATION_SIZE)
    ax4_1.margins(x=0)
    ax4_1.set_xlim(first(t2plot), last(t2plot))
    ax4_1.set_ylim(0,plot_max)
    ax4_1.set_xticks(t2plot[10:24:end]) 
    ax4_1.set_xticklabels(t2plot_label[10:24:end], fontsize=TICK_SIZE)
    ax4_1.tick_params(axis="x", pad=2)
    ax4_1.tick_params(axis="x", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    ax4_1.tick_params(axis="y", labelsize=TICK_SIZE,width=AXIS_LINEWIDTH)
    ax4_1.tick_params(axis="y", labelleft=false)
    add_label_atleftup(ax4_1,"(e) Haikou","left")

    legend_handles, legend_labels = ax3.get_legend_handles_labels()
    fig.legend(legend_handles, legend_labels, loc="lower center", bbox_to_anchor=(0.29, 0.015),
        ncol=3, frameon=false, prop=font_legend_small, handlelength=0.9,
        handletextpad=0.25, columnspacing=0.45, labelspacing=0.10, borderaxespad=0.0)
    fig.text(0.018, 0.27, raw"$\mathrm{PM}_{10}$"*raw" (μg m$^{-3}$)", rotation=90,
        va="center", ha="center", fontproperties=font_axis)




    add_label_atleftup(ax5,"2025-4-11 22:00","right",xpos=0.98)
    draw(ax6,ax6data,height,xtick,xlabel_sub,false)
    add_label_atleftup(ax6,"2025-4-12 2:00","right",xpos=0.98)
    

    temp_font = font_propbd.copy()
    temp_font.set_size(PANEL_TITLE_SIZE)
    ax5.set_title("(f) Vertical Distribution of Dust", pad=4,loc="left",fontproperties=temp_font)
    draw(ax7,ax7data,height,xtick,xlabel_sub,false)
    add_label_atleftup(ax7,"2025-4-12 6:00","right",xpos=0.98)
    draw(ax8,ax8data,height,xtick,xlabel_sub,false)
    add_label_atleftup(ax8,"2025-4-12 10:00","right",xpos=0.98)
    draw(ax9,ax9data,height,xtick,xlabel_sub,false)
    add_label_atleftup(ax9,"2025-4-12 14:00","right",xpos=0.98)
    draw(ax10,ax10data,height,xtick,xlabel_sub,true)
    add_label_atleftup(ax10,"2025-4-12 18:00","right",xpos=0.98)
    plt.gcf()

end

output_dir = "/home/zsm/IETM/figure/PapreFigure/20260826"
mkpath(output_dir)
plt.savefig(joinpath(output_dir, "PaperFigure_VerticaltransportImproves.png"),dpi=600,bbox_inches="tight",pad_inches=0.01)
plt.savefig(joinpath(output_dir, "PaperFigure_VerticaltransportImproves.pdf"),bbox_inches="tight",pad_inches=0.01)
plt.close("all")
println("Saved PaperFigure_VerticaltransportImproves.png and .pdf")
exit()


# tStart_all=ZonedDateTime.(DateTime(2023,4,11):Hour(1):DateTime(2023,4,13),tz"Asia/Shanghai")
# tStart_RUF=tStart_all[1]

# #new
# animateAFore(tStart_all;rufFolder="/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d/DATA/RUF/$(Dates.format(astimezone(tStart_all[1],tz"Asia/Shanghai"),"yyyymmddHH"))")
# #old
# # animateAFore(tStart_all;rufFolder="/sh2/home/3clear/wuhj/paper_IETMsand/backup/IETM3d_DATA/DATA/RUF/$(Dates.format(astimezone(tStart_all[1],tz"Asia/Shanghai"),"yyyymmddHH"))")
