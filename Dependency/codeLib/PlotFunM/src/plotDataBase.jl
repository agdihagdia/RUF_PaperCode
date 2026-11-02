using DataBaseM,LocalizerM,NaNStatistics,PyCall,Dates,TimeZones,utilities
export plotTargetAndNearby,plotObsByCity,plotObsScatter,plotObsScatterAnimation

function plotTargetAndNearby(Times,pollName,staID_target;staID_nearby=nothing,tableName="original.pollutant")
    staID_target,Times=in2vec.([staID_target,Times])
    if isnothing(staID_nearby)
        staInfo_valid=getValidStaInfo(Times[1:min(10,length(Times))],pollName;tableName=tableName)
        lonlat=Matrix(staInfo_valid[:,["lon","lat"]])
        locM=localizerM(lonlat,lonlat)
        ind_target=indexin(staID_target,staInfo_valid.staID)
        staInfo_nearby=staInfo_valid[LocalizerM.getNearest(locM,ind_target,1:11),:]
        # staID_ref=staID_valid[LocalizerM.getNearest(locM,1,1:21)]
        # staID_nearby=staInfo[getLocCoef(localizer,ind,20e3),"staID"]
    end
    # select Times and staID and pollName
    obs_target=vec(getObsPoll(Times,pollName,staID_target;tableName=tableName))

    # prepare titleName
    # staInfo_target=DataBaseM.getStaInfo(staID_target)
    titleName=join([staInfo_nearby.province[1],staInfo_nearby.city[1],staInfo_nearby.staName[1],staInfo_nearby.staID[1]],"  ")

    # obs_nearby=getObsPoll(Times,pollName,staID_nearby)
    obs_nearby=dropdims(getObsPoll(Times,pollName,staInfo_nearby.staID;tableName=tableName),dims=2)
    plt.figure(figsize=(12,4))
    h_target=plt.plot(Times,obs_target,linewidth=1,zorder=10,"r",label="target")
    h_nearby=plt.plot(Times,obs_nearby,linewidth=1,zorder=1,c="#95a5a6",label="nearby")
    # plt.legend()
    plt.legend(handles=vcat(h_target,h_nearby[1]))
    plt.xlim(Times[[1,end]])
    plt.title(titleName,family="SimHei")
end

"""

sample usage:
Times=collect(DateTime(2023,3,1):Hour(1):DateTime(2023,4,20))
city2plot="北京市"
plt.figure(figsize=(12,4))
plotObsByCity(Times,"PM10",city2plot;tableName="original.pollutant")
plt.gcf()
"""
function plotObsByCity(Times,pollName,cityName;useNation=true,tableName="verified.pollutant_julia")
    staInfo=DataBaseM.getStaInfo()
    staInfo=staInfo[isequal.(staInfo.city,cityName),:]
    if useNation
        staInfo=staInfo[isequal.(staInfo.type,"Nation"),:]
    end
    obs=dropdims(getObsPoll(Times,pollName,staInfo.staID;tableName=tableName),dims=2)
    plt.plot(Times,obs,linewidth=1,zorder=1,label=staInfo.staName)
    L = plt.legend()
    plt.setp(L.texts, family="SimHei")
    plt.xlim(Times[[1,end]])
    addPollYLabel(pollName)
    plt.title(cityName * " " * varNameDisp(pollName),family="SimHei")
end
plotObsScatter(Times::Union{DateTime,ZonedDateTime},pollName;useNation=true,tableName="verified.pollutant_julia")   = plotObsScatter([Times],pollName;useNation=useNation,tableName=tableName) 
function plotObsScatter(Times::AbstractVector,pollName;useNation=true,tableName="verified.pollutant_julia")
    staInfo=DataBaseM.getStaInfo()
    if useNation
        staInfo=staInfo[isequal.(staInfo.type,"Nation"),:]
    end
    obs=vec(nanmean(getObsPoll(Times,pollName,staInfo.staID;tableName=tableName),dims=1))
    h_s=plt.scatter(staInfo.lon,staInfo.lat,c=obs,cmap="jet")
    h_china=plotChina()
    h_bar=plt.colorbar()
    h_bar.ax.set_title(getUnitStr(pollName))
    plt.axis([72, 137, 17, 54])
    plt.title( varNameDisp(pollName),family="SimHei")
    return h_s,h_china,h_bar
end

function plotObsScatterAnimation(Times,pollName;useNation=false,tableName="verified.pollutant_julia",
    cLimit=nothing,dpi=150,fileName="animation.mp4",axLimit=[72, 137, 17, 54],figsize=(7,5))
    plt=get_plt()
    animation=pyimport("matplotlib.animation")    
    staInfo=DataBaseM.getStaInfo()
    if useNation
        staInfo=staInfo[isequal.(staInfo.type,"Nation"),:]
    end
    obs=dropdims(getObsPoll(Times,pollName,staInfo.staID;tableName=tableName),dims=2)
    
    if isnothing(cLimit)
        cLimit=[0,nanpctile(obs,99)]
    end
    
    fig = plt.figure(figsize=figsize,dpi=dpi)
    gs = fig.add_gridspec(1,1,left=0.05,right=0.999,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)
    function plot_aframe(i)
        n=i+1
        plt.clf()
        data2plot=obs[n,:]
        ax=plt.subplot(gs[1])
        ind=sortperm(data2plot)
        h_s=plt.scatter(staInfo.lon[ind],staInfo.lat[ind],c=data2plot[ind],s=5,cmap="jet")
        h_china=plotChina()
        h_bar=plt.colorbar()
        h_bar.ax.set_title(getUnitStr(pollName))
        plt.axis(axLimit)
        plt.clim(cLimit)
        plt.title( varNameDisp(pollName) * " at " * Dates.format(Times[n],"yyyy-mm-dd HH:00:00"),family="SimHei")
        return [h_s , h_china , h_bar]
    end
    
    ani = animation.FuncAnimation(fig, plot_aframe, frames=length(Times), interval=200)
    ani.save(fileName)
    return ani    
end
