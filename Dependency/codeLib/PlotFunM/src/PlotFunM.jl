module PlotFunM
using MyPlotM
export plotCnShp,getProj,getSouthChinaAxes,setAxes_China,resize_colobar,get_plt,set_dateAxis,changeFont,get_font,varNameDisp,addPollYLabel
export cm2inch,inset_axes,add_subAxes,get_mticker,add_SouthChinaAxes,plotChina,setDateTick,quickPlot
# plt=get_plt()
import PyPlot as plt
include("plotDataBase.jl")
include("plotStaData.jl")
include("plotQualityControl.jl")
include("plotNaqpms.jl")
include("plotSandStorm.jl")
include("plotObsFromTxt.jl")
include("plotSHAP.jl")
end # module PlotFunM
