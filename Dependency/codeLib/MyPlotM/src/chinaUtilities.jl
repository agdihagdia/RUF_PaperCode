
using PyCall,MAT
import PyPlot as plt

export plotCnShp,getProj,getSouthChinaAxes,setAxes_China,plotChina,add_SouthChinaAxes,add_subAxes,addChinaColorbar,getChinaShpDir
function plotCnShp(ax=plt.gca())
    

    ccrs = pyimport("cartopy.crs")
    cfeature=pyimport("cartopy.feature")
    shpreader=pyimport("cartopy.io.shapereader")
    #  gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=true, dms=true, x_inline=false, y_inline=false, linewidth=1, color="black", alpha=0.5, linestyle="--") # 添加经纬度网格，这里设定为黑色虚线半透明 
    # ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=true) # 添加经纬度网格，这里设定为黑色虚线半透明 
    # gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true, x_inline=false, y_inline=false, linewidth=0.33, color='k',alpha=0.5)
    china = shpreader.Reader("$(dirname(pathof(@__MODULE__)))/../china_basic_map/bou2_4l.dbf").geometries() # 读取中国的国境线、省界及九段线等信息，这里读取 dbf 或 shp 文件都可 

    cn_boundary=ax.add_geometries(china,ccrs.PlateCarree(), facecolor="none", edgecolor="gray", zorder = 1,linewidth=0.5) # 添加 shp 信息到经纬度投影
    return cn_boundary
end

function getChinaShpDir()
    china = "$(dirname(pathof(@__MODULE__)))/../china_basic_map/bou2_4l.shp"
    return china
end
function addChinaColorbar(ax=plt.gca();pos=[0.9, 0.02, 0.03, 0.55],extend="both")    
    # h_bar=plt.colorbar(extend="both",location="right",shrink=0.8) 
    axins = ax.inset_axes(pos)
    h_bar = plt.colorbar(cax=axins,extend=extend)
    return h_bar
end

function plotChina(ax=plt.gca();linewidth=1)
    vars=matread("$(dirname(pathof(@__MODULE__)))/../china_basic_map/cnShp.mat")
    hp=plt.plot(vec(vars["cnLong"]),vec(vars["cnLat"]),color="gray",linewidth=linewidth)
    return hp
end

function getSouthChinaAxes(fig=plt.gcf(),axIn=plt.gca();pos="LowerRight")
    # if isempty(posIn)
        posIn=axIn.get_position()
        in_xmax,in_ymin,in_width,in_height=posIn.xmax,posIn.ymin,posIn.width,posIn.height
    # else
    #     in_xmin,in_ymin,in_width,in_height=posIn...
    #     in_xmax=in_xmin+in_width
    # end
    ccrs = pyimport("cartopy.crs")
    hight2width_fig=fig.figbbox.height/fig.figbbox.width
    width=in_width/8
    height=width*2/hight2width_fig
    if pos == "LowerRight"
        position_southChina=[in_xmax-width,in_ymin,width,height]
    elseif  pos == "LowerLeft"
        position_southChina=[posIn.xmin,in_ymin,width,height]
    else
        @error "pos should be LowerRight or LowerLeft "
    end
    println(position_southChina)
    ax_southChina = fig.add_axes(position_southChina, projection=axIn.projection)
    ax_southChina.set_extent([106, 120, 2, 24],crs=ccrs.PlateCarree())
    return ax_southChina
end

function add_SouthChinaAxes(axIn=plt.gca();pos="LowerRight",disp=false)
    fig=axIn.get_figure()
    posIn=axIn.get_position()
    in_xmax,in_ymin,in_width,in_height=posIn.xmax,posIn.ymin,posIn.width,posIn.height
    ccrs = pyimport("cartopy.crs")
    hight2width_fig=fig.figbbox.height/fig.figbbox.width
    if pos == "LowerRight"
        width=in_width/7
        height=width*2/hight2width_fig
        position_southChina=[in_xmax-width,in_ymin-height/8,width,height]
    elseif  pos == "LowerLeft"
        width=in_width/6
        height=width*2/hight2width_fig
        position_southChina=[posIn.xmin,in_ymin-height/7.5,width,height]
    elseif length(pos) == 4
        position_southChina=pos
    else
        @error "pos should be LowerRight or LowerLeft "
    end
    if disp;println(position_southChina);end
    ax_southChina=fig.add_axes([0, 0, 0.1, 0.1], projection=axIn.projection)
    ax_southChina.set_extent([106, 120, 2, 24],crs=ccrs.PlateCarree())
    ax_southChina.set_position(position_southChina)
    return ax_southChina
end

function add_subAxes(ax,pos_relative)
    
    ax_sub = ax.get_figure().add_axes([0, 0, 0.1, 0.1])
    # plt.draw()
    posn = ax.get_position()
    ax_sub.set_position([posn.x0 + posn.width *pos_relative[1],          posn.y0 + posn.height *pos_relative[2], 
                  posn.width *pos_relative[3], posn.height*pos_relative[4]])
    return ax_sub
end


function setAxes_China(ax=plt.gca())
    ccrs = pyimport("cartopy.crs")
    ax.set_extent([80, 130, 15, 55],crs=ccrs.PlateCarree())
end