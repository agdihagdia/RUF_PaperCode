using MyPlotM,NaNStatistics,PyCall,Dates,utilities

"""

sample usage:
fileName="/public/obs/20230517080000.txt"
pollName="PM10"
plt=get_plt()

fig = plt.figure(figsize=(7,5),dpi=150)
plt.clf()
gs = fig.add_gridspec(1,1,left=0.05,right=0.999,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)
plt.subplot(gs[1])
plotObsScatter(fileName,pollName)
plt.gcf()
"""
function plotObsScatter(fileName::AbstractString,pollName;cLimit=nothing,cmap="jet")
    df_obs=readAObsTxt(fileName)
    ind_targetPoll=indexin([lowercase(pollName)],lowercase.(names(df_obs)))
    lon=df_obs.lon;lat=df_obs.lat;
    data2plot=vec(Matrix(df_obs[!,ind_targetPoll]))
    if isnothing(cLimit);cLimit=[0,nanpctile(data2plot,98)];end
    
    # cmap_cur = plt.cm.get_cmap(cmap,len(boundaries))
    h_s=plt.scatter(lon,lat,c=data2plot,s=3,cmap=cmap)
    plt.clim(cLimit);
    h_china=plotChina()
    h_bar=addChinaColorbar()
    # h_bar=plt.colorbar(extend="both") 
    h_bar.ax.set_title("   " * getUnitStr(pollName))
    plt.axis([72, 137, 17, 54])
    plt.title( varNameDisp(pollName))
    return h_s,h_china,h_bar
end

