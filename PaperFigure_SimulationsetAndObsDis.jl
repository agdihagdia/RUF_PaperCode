using Revise,MyPlotM,GradDataM,MeshDataM,StaDataM,StaDataTool,MeshDataTool,NaqpmsUtilities,utilities,Dates,TimeZones,NaNStatistics
using GradDataM,Dates,TimeZones,Statistics,ZZDBM,DASystemM,PlotFunM
using PyCall
using Conda
includet("$(pathof(DASystemM))/../../system/"*"rapidUpdatedForecast_3d.jl")
# Conda.add("owslib.wms")
ccrs = pyimport("cartopy.crs")
matplotlib = pyimport("matplotlib.pyplot")
cfeature = pyimport("cartopy.feature")
shpreader = pyimport("cartopy.io.shapereader")
wms = pyimport("owslib.wms")
np = pyimport("numpy")
# pd = pyimport("pandas")
matplotlib = pyimport("matplotlib")
plt = matplotlib.pyplot
ticker=matplotlib.ticker
WebMapService=wms.WebMapService
fm=matplotlib.font_manager
mpimg = pyimport("matplotlib.image")
np = pyimport("numpy")
PIL = pyimport("PIL.Image")



function read_tiff(filepath)
    img = PIL.open(filepath)
    img_array = np.array(img)
    return img_array
end

function addGrids(ax)
    gl = ax.gridlines(crs=PlaCar,draw_labels=true, dms=true, x_inline=false, y_inline=false, linewidth=0.5, color=[0.1,0.1,0.1], alpha=0.5,zorder=2,linestyle="--") # 添加经纬度网格，这里设定为黑色虚线半透明 
    gl.top_labels  = false;    gl.right_labels = false;   gl.rotate_labels = false
    gl.xlocator=ticker.FixedLocator([40,60,80,100,120,140,160,180])
    gl.xlabel_style = Dict(:size => 5.5, :fontname => "Helvetica")  
    gl.ylabel_style = Dict(:size => 5.5, :fontname => "Helvetica")  
end 


function add_SouthChinaAxes(axIn)
    fig=axIn.get_figure()
    posIn=axIn.get_position()
    in_xmax,in_ymin,in_width,in_height=posIn.xmax,posIn.ymin,posIn.width,posIn.height
    hight2width_fig=fig.figbbox.height/fig.figbbox.width

    width=in_width/7
    height=width*2/hight2width_fig
    position_southChina=[in_xmax-width,in_ymin-height/5.1,width,height]

    # println(position_southChina)
    ax_southChina=fig.add_axes([0, 0, 0.1, 0.1], projection=axIn.projection)
    ax_southChina.set_extent([106, 120, 2, 24],crs=ccrs.PlateCarree())
    ax_southChina.set_position(position_southChina)
    cn_boundary=ax_southChina.add_geometries(shpreader.Reader("/sh2/home/3clear/wuhj/software/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
    PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
    ax_southChina.scatter(Bsta_lon,Bsta_lat,marker=".",label="National Monitoring Station",zorder=3,s=0.1,transform=PlaCar)
    ax_southChina.scatter(Csta_lon,Csta_lat,marker=".",label="County Monitoring Station",zorder=2,s=0.05,transform=PlaCar)
    ax_southChina.scatter(Dsta_lon,Dsta_lat,marker=".",label="Township Monitoring Station",zorder=1,s=0.01,transform=PlaCar)
    return ax_southChina
end


# filepath="/sh2/home/3clear/wuhj/paper_IETMsand/figure_zsm/data/SR_LR.tif"
# img_array = read_tiff(filepath)
PlaCar=ccrs.PlateCarree()
# staInfo=readStaInfo("/sh2/home/3clear/wuhj/online_IETMsand/staInfo.csv")
staInfo=staDATA("/sh2/home/3clear/wuhj/paper_IETMsand/ori202504.mat")
lonM,latM=load("datagrid.jld2","lonM","latM")


cities_second = [
    "呼和浩特",
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
    "清远",
    "茂名"
    
]
cities_first=["太原","北京","南京","海口"]
df_sub = filter(:city => c -> c in cities_first, staInfo.staInfo)
cities_first_mean = combine(groupby(df_sub, :city),
                  :lon => mean => :lon_mean,
                  :lat => mean => :lat_mean)
fcities_lon=cities_first_mean.lon_mean
fcities_lat=cities_first_mean.lat_mean
df_sub = filter(:city => c -> c in cities_second, staInfo.staInfo)
cities_second_mean = combine(groupby(df_sub, :city),
                  :lon => mean => :lon_mean,
                  :lat => mean => :lat_mean)
scities_lon=cities_second_mean.lon_mean
scities_lat=cities_second_mean.lat_mean



# lon=staInfo.lon
# lat=staInfo.lat
# type=staInfo.type

Bsta_lon=staInfo.staInfo.lon[in.(staInfo.staInfo.type,[["B"]])]
Bsta_lat=staInfo.staInfo.lat[in.(staInfo.staInfo.type,[["B"]])]
Csta_lon=staInfo.staInfo.lon[in.(staInfo.staInfo.type,[["C"]])]
Csta_lat=staInfo.staInfo.lat[in.(staInfo.staInfo.type,[["C"]])]
Dsta_lon=staInfo.staInfo.lon[in.(staInfo.staInfo.type,[["D"]])]
Dsta_lat=staInfo.staInfo.lat[in.(staInfo.staInfo.type,[["D"]])]
Csta_lon_sub=Csta_lon[1:2:end]
Csta_lat_sub=Csta_lat[1:2:end]
Dsta_lon_sub=Dsta_lon[1:4:end]
Dsta_lat_sub=Dsta_lat[1:4:end]



vars = matread("/sh2/home/3clear/wuhj/software/codeLib/MyPlotM/china_basic_map/cnShp.mat")
# 提取经纬度
lon = vec(vars["cnLong"])
lat = vec(vars["cnLat"])
# lam_lon,lam_lat=proj.transform_points(ccrs.PlateCarree(), lon, lat)


# boundary_d2_lon = vcat(lonM[:,1],lonM[end,:],reverse(lonM[:,end]),reverse(lonM[1,:]))
# boundary_d2_lat = vcat(latM[:,1],latM[end,:],reverse(latM[:,end]),reverse(latM[1,:]))



#new
grd1=gradData("/sh2/home/3clear/post/Mdata/china_3_27km/naqp/grid/grid.d1.ctl")
grd2=gradData("/sh2/home/3clear/post/Mdata/china_3_27km/naqp/grid/grid.d2.ctl")
lonM1=read(grd1,"XLONG")[:,:,1]
latM1=read(grd1,"XLAT")[:,:,1]
lonM2=read(grd2,"XLONG")[:,:,1]
latM2=read(grd2,"XLAT")[:,:,1]
boundary_d1_lon = vcat(lonM1[:,1],lonM1[end,:],reverse(lonM1[:,end]),reverse(lonM1[1,:]))
boundary_d1_lat = vcat(latM1[:,1],latM1[end,:],reverse(latM1[:,end]),reverse(latM1[1,:]))
boundary_d2_lon = vcat(lonM2[:,1],lonM2[end,:],reverse(lonM2[:,end]),reverse(lonM2[1,:]))
boundary_d2_lat = vcat(latM2[:,1],latM2[end,:],reverse(latM2[:,end]),reverse(latM2[1,:]))



minimum(boundary_d2_lat)
maximum(boundary_d2_lon)


function add_label_atleftup(ax, label; xpos=0.02, ypos=0.95, fontsize=6)
    ax.text(xpos, ypos, label,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="bold",
            va="top",
            ha="left",
            fontproperties=font_propbd)
end
#绘图
begin
plt=get_plt()
fig=plt.figure(figsize=(6,3),dpi=600);plt.clf()
plt.subplots_adjust(left=0.06, bottom=0.03, right=0.989,  top=0.98, wspace=0.3,  hspace=0.45)
matplotlib.rcParams["font.family"] = "Helvetica"
proj=ccrs.LambertConformal(central_longitude=105, central_latitude=83,standard_parallels=[25,45],globe=nothing)
font_path = "/sh2/home/3clear/wuhj/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
font_prop = fm.FontProperties(fname=font_path,size=6.5)
font_pathbd = "/sh2/home/3clear/wuhj/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
font_propbd = fm.FontProperties(fname=font_pathbd,size=5.5)
ax1=fig.add_axes([0.1, 0.3, 0.40, 0.55], projection=proj)
# ax1=fig.add_subplot(121,projection=proj)
# plotChina(ax1)
cn_boundary=ax1.add_geometries(shpreader.Reader("/sh2/home/3clear/wuhj/software/codeLib/MyPlotM/china_basic_map/cntry02.dbf").geometries(),
PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)
# ax1.imshow(img_array, cmap="terrain", transform=ccrs.PlateCarree()) 
# ax1.set_facecolor([0, 170/255 ,255/255])
# ax1.set_extent([63.4, 146.6, 2.2, 63.4],crs=PlaCar)
ax1.set_extent([63.4, 146.6, 0, 65],crs=PlaCar)
# ax1.add_wms(wms="http://portal.cubewerx.com/cubewerx/cubeserv.cgi/simple",layers=["Foundation.GTOPO30"])
ax1.add_wms(wms="https://basemap.nationalmap.gov/arcgis/services/USGSImageryOnly/MapServer/WMSServer",
           layers=["0"])
ax1.plot(boundary_d2_lon,boundary_d2_lat,"r", linewidth=0.8,transform=PlaCar)
ax1.plot(boundary_d1_lon,boundary_d1_lat,"y",linewidth=0.8,transform=PlaCar)
addGrids(ax1)
# ax1.set_title("a) Research Area",fontsize=6,fontweight="heavy",pad=3,fontproperties=font_prop)
# add_label_atleftup(ax1,"a) Research Area")
ax1.text(0.02, 0.95, "(a) Research Area",
transform=ax1.transAxes,
fontsize=6,
fontweight="bold",
va="top",
ha="left",
fontproperties=font_propbd,
color="white")
ax1.set_aspect("auto")
# ax2=fig.add_subplot(122,projection=proj)
ax2=fig.add_axes([0.56, 0.3, 0.40, 0.55], projection=proj)
ax2.set_extent([72.3, 136.3, 16.9, 54.2],crs=PlaCar)

# cn_boundary=ax2.add_geometries(shpreader.Reader("/sh2/home/3clear/wuhj/software/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
cn_boundary=ax2.add_geometries(shpreader.Reader("/sh2/home/3clear/wuhj/software/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf").geometries(),
PlaCar, facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5)

ax2.scatter(Bsta_lon,Bsta_lat,marker=".",label="National Station",zorder=3,s=0.1,transform=PlaCar)
ax2.scatter(Csta_lon,Csta_lat,marker=".",label="County Station",zorder=2,s=0.05,transform=PlaCar)
ax2.scatter(Dsta_lon,Dsta_lat,marker=".",label="Township Station",zorder=1,s=0.01,transform=PlaCar)
ax2.scatter(fcities_lon,fcities_lat,marker="*",label="Selected major cities",color="r",zorder=4,s=1.0,transform=PlaCar)
ax2.scatter(scities_lon,scities_lat,marker="*",label="Selected secondary cities",color="r",zorder=4,s=0.02,transform=PlaCar)

addGrids(ax2)
add_SouthChinaAxes(ax2)
# ax2.set_title("b) Distribution of Sites within the Research Area",fontsize=6,fontweight="heavy",pad=3,fontproperties=font_prop)
add_label_atleftup(ax2,"(b) Distribution of Sites")
ax2.legend(fontsize=4.5,loc="lower left")
ax2.set_aspect("auto")

# ax2.set_title("National-County-Township Monitoring Stations",fontsize=4)


# plt.savefig("/sh2/home/3clear/wuhj/paper_IETMsand/figure_zsm/figure2_8_16_test.png",dpi=600)
# plt.tight_layout()
plt.gcf()
end
plt.savefig("/sh2/home/3clear/wuhj/paper_IETMsand/figure_zsm/figure2_9_5.png",dpi=600)
plt.close()