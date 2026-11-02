module MyPlotM

include("customColormap.jl")
include("dateUtilities.jl")
include("pollutantUtilities.jl")
include("chinaUtilities.jl")
include("taylorDiagram.jl")
include("textUtilities.jl")
export getProj,resize_colobar,get_plt,changeFont,get_font
export cm2inch,inset_axes,get_mticker,quickPlot
using PyCall
import PyPlot as plt
# 2022-12-20 add plotChina to plot China bounding using plot intead of cartopy
# sample usage ==========================
# proj,PlaCar=getProj("Robinson")

function cm2inch(in)
    return in/2.54
end

function get_plt()
    changeFont(plt)
    plt.rc("xtick",direction="in")  
    plt.rc("ytick",direction="in")  
    return plt
end

# ========================= for fonts ====================================
function changeFont(plt,fontName="Microsoft YaHei";fontsize=8)#"DejaVu Sans"
    # when Font family ['SimHei'] not found, 
    # 1. copy the font file ../font/SimHei to  /public/home/cnemc_user1/software/anaconda3/envs/conda_jl/lib/python3.10/site-packages/matplotlib/mpl-data/fonts/ttf/
    # 2. delete catch $ rm ~/.cache/matplotlib/fontlist-v*
    # using PyCall
    # m=pyimport("matplotlib")
    # print(m.matplotlib_fname())
    # /public/home/cnemc_user1/software/anaconda3/envs/conda_jl/lib/python3.10/site-packages/matplotlib/mpl-data/matplotlibrc
    
    plt.rc("font",family=fontName)   # 步骤一（设置字体）
    plt.rc("font",size=fontsize)   # 步骤一（设置字体）
    plt.rc("axes",unicode_minus=false)   # 步骤二（解决坐标轴负数的负号显示问题）
    # plt.ylabel("直接信任度值", fontsize=13, fontdict={"family": "KaiTi", "size": 15, "color": "k"}) 
    return plt
end
function inset_axes(argIn...)
    inset_locator  = pyimport("mpl_toolkits.axes_grid1.inset_locator")
    return inset_locator.inset_axes(argIn...)
end
function get_font(fontName)
    fmanager = pyimport("matplotlib.font_manager")
    if lowercase.(fontName) == "arial"
        my_font = fmanager.FontProperties(fname="$(dirname(pathof(@__MODULE__)))/../font/ARIALN.TTF", size=10)
    elseif lowercase.(fontName) == "simsun"
        my_font = fmanager.FontProperties(fname="$(dirname(pathof(@__MODULE__)))/../font/SimSun.ttf", size=10)
    elseif lowercase.(fontName) == "simsunb"
        my_font = fmanager.FontProperties(fname="$(dirname(pathof(@__MODULE__)))/../font/simsunb.ttf", size=10)
    end

    return my_font
end
function resize_colobar(ax,cbar)
    # Tell matplotlib to re-draw everything, so that we can get
    # the correct location from get_position.
    posn=ax.get_position()
    cbar.ax.set_position([posn.x0 + posn.width + 0.01, posn.y0, 0.04, posn.height])
end
function get_mticker()
    return pyimport("matplotlib.ticker")
end

"""
sample usage: 
proj,PlaCar=getProj()
"""
function getProj(projName="default")
    ccrs = pyimport("cartopy.crs")
    #  proj = ccrs.Robinson() # 
    if projName == "default"
        proj = ccrs.LambertConformal(central_longitude=110, central_latitude=90,standard_parallels=(25, 47))
    elseif projName == "Robinson"
        proj = ccrs.Robinson()
    elseif projName == "Orthographic"
        proj = ccrs.Orthographic()
    elseif projName == "PlateCarree"
        proj = ccrs.PlateCarree()
    end
    PlaCar=ccrs.PlateCarree()
    return proj,PlaCar
end

function quickPlot(x_all...)
    fig=plt.figure(figsize=(4,2 * length(x_all)));
    gs = fig.add_gridspec(nrows=length(x_all), ncols=1, left=0.05,right=0.98,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)
    for (n,x) in enumerate(x_all)
        ax=plt.subplot(gs[n])
        plt.plot(x)
    end
    plt.gcf()
end

end # module
