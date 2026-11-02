using StaDataM,Dates,TimeZones,NaNStatistics
export plotASta,animateStaScatter,plotStatistics,plotObsScatter
function plotASta(staObjAll,tagS,staID2plot;Times=staObjAll[1].times,varNames=staObjAll[1].varNames)
    staInfo2plot=getStaInfo(staObjAll[1],staID2plot)
    # prepare titleName ----------
    titleName=""
    for field in [ :province  ,:city,:staName, :staID ]
        if hasproperty(staInfo2plot,field)
            titleName *= staInfo2plot[1,field] * " "
        end
    end
    # prepare cMaps and t2plot
    cMapS=["#ff4757","#ffa502","#1e90ff","#2ed573","#ff7f00","#ffff33","#a65628","#f781bf","#999999"]
    t2plot=DateTime.(Times)
    # ------------------ plot -----------
    fig=plt.figure(figsize=(8,1.8*length(varNames)));
    gs = fig.add_gridspec(nrows=6, ncols=1, left=0.06,right=0.99,bottom=0.02,top=0.96, hspace=0.13, wspace=0.15)
    for (n_var,varName) in enumerate(varNames)
        ax=plt.subplot(gs[n_var])
        for (n_tag,tag) in enumerate(tagS)
            data2plot=vec(read(staObjAll[n_tag];Times=Times,varNames=varName,staIDs=staInfo2plot.staID))
            if n_tag ==1 
                plt.plot(t2plot,data2plot,label=tag,color=cMapS[n_tag],marker="o",linestyle="",markersize=3,zorder=0)
            else
                plt.plot(t2plot,data2plot,label=tag,color=cMapS[n_tag])
            end
            plt.xlim(t2plot[[1,end]])
            addPollYLabel(varName)
            plt.legend()
        end
    end
    ht=plt.figtext(0.5,0.97,titleName,weight="bold",size=13,family="SimHei",ha="center")
end

function animateStaScatter(staObj;Times=staObj.times,varName=staObj.varNames[1],staInfo=staObj.staInfo,
    cLimit=[0,1.5*nanpctile(read(staObj;Times=Times,varNames=varName,staIDs=staInfo.staID),99)], 
     cMap="jet",  dpi=100,axLimit=[72, 137, 17.5, 55],figsize=(6,4.9),markersize=1,
    outFileName="animation of $varName started at $(Dates.format(Times[1],"yyyy-mm-dd HH")).mp4")   
    
    data_sta=dropdims(read(staObj;Times=Times,varNames=varName,staIDs=staInfo.staID),dims=2)

    plt=get_plt()
    animation=pyimport("matplotlib.animation")
    fig = plt.figure(figsize=figsize,dpi=dpi)
    gs = fig.add_gridspec(1,1,left=0.05,right=0.98,bottom=0.04,top=0.93, hspace=0.1, wspace=0.1)
    function plot_aframe(i)
        n=i+1
        plt.clf()
        ax=plt.subplot(gs[1])
        hp=plt.scatter(staInfo.lon,staInfo.lat,s=markersize,c=data_sta[n,:],cmap=cMap)
        plt.clim(cLimit);
        hcn=plotChina()
        plt.axis(axLimit)
        h_bar=addChinaColorbar()
        h_bar.ax.set_title("   " * getUnitStr(varName))

        plt.title("$varName at $(Dates.format(Times[n],"yyyy-mm-dd HH"))")
        hp_all=[hp,hp]
        return hp_all
    end
    
    ani = animation.FuncAnimation(fig, plot_aframe, frames=length(Times), interval=200, blit=true)    
    ani.save(outFileName)
    return ani
end

function plotObsScatter(staObj::StaDataM.staDATAStru;t=staObj.times[1],varName=staObj.varNames[1],staInfo=staObj.staInfo,
    cLimit=[0,nanpctile(staObj.staDataAll,99.9)],  cMap="jet",axLimit=[72, 137, 17.5, 55], newFig=false)

    if newFig
        plt=get_plt()
        fig = plt.figure(figsize=(6,5),dpi=100)
        gs = fig.add_gridspec(nrows=1, ncols=1,left=0.05,right=0.98,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)
        ax=plt.subplot(gs[1])
    else
        ax=plt.gca()
    end
    
    data2plot=vec(read(staObj;Times=t,varNames=varName,staIDs=staInfo.staID))
    hp=plt.scatter(staInfo.lon,staInfo.lat,c=data2plot,s=2,cmap=cMap)
    plt.clim(cLimit);
    h_china=plotChina()
    h_bar=addChinaColorbar()
    h_bar.ax.set_title("   " * getUnitStr(varName))
    plt.axis(axLimit)
    plt.title("$(varNameDisp(varName)) at $(Dates.format(t,"yyyy-mm-dd HH"))")
    
    return hp

end
function plotStatistics(staObs,staSim;TimeS=staObs.times, varName=staObs.varNames,staIDs=staObs.staInfo.staID,    
    gs = plt.gcf().add_gridspec(nrows=2, ncols=2,left=0.07,right=0.99,bottom=0.05,top=0.92, hspace=0.2, wspace=0.1)
    )
    # plot distribution of mean bias rmse and correlation between two staDataObjS
    # staDataObjS contain {ObsStaDATA simStaDATA}. only one varName is supported 

    
    # input check
    staInfo=StaDataM.getStaInfo(staObs,staIDs)
    
    # ************* prepare data to plot ***********
    observation=dropdims(read(staObs;Times=TimeS,staIDs=staIDs,varNames=varName),dims=4)
    simulation=dropdims(read(staSim;Times=TimeS,staIDs=staIDs,varNames=varName),dims=4)
    obsmean=vec(nanmean(observation,dims=1))
    simmean=vec(nanmean(simulation,dims=1))
    bias=vec(nanmean(simulation-observation,dims=1))
    rmse=vec(sqrt.(nanmean((simulation-observation).^2,dims=1)))
    corrLoc=[ nancor(observation[:,1,n],simulation[:,1,n]) for n in axes(observation,3)]
    # collect statistics to staInfo
    staInfo.obsMean=obsmean;staInfo.simMean=simmean;staInfo.bias=bias;staInfo.rmse=rmse;staInfo.corr=corrLoc

    function setAxes(ax)
        # get the axLimit
        xLimit=ax.get_xlim();yLimit=ax.get_ylim();
        plotChina()
        ax.set_xlim(xLimit);ax.set_ylim(yLimit)
        plt.colorbar()
    end
    plt.clf()
    # plot observation mean
    ax=plt.subplot(gs[1])
    plt.scatter(staInfo.lon,staInfo.lat,c=staInfo.obsMean,s=2,cmap="rainbow")
    setAxes(ax)
    plt.title("观测平均")    
    # plot bias
    ax=plt.subplot(gs[2])    
    plt.scatter(staInfo.lon,staInfo.lat,c=staInfo.bias,s=2,cmap="seismic",vmin=-1.1*nanmaximum(abs.(staInfo.bias)),vmax=1.1*nanmaximum(abs.(staInfo.bias)))
    setAxes(ax)
    plt.title("偏差")
    # plot rmse    
    ax=plt.subplot(gs[3])
    plt.scatter(staInfo.lon,staInfo.lat,c=staInfo.rmse,s=2,cmap="Oranges")
    setAxes(ax)
    plt.title("均方根误差")  
    # plot corr    
    ax=plt.subplot(gs[4])
    plt.scatter(staInfo.lon,staInfo.lat,c=staInfo.corr,s=2,cmap="BrBG",vmin=-1,vmax=1)
    setAxes(ax)
    plt.title("相关系数")  

    # plt.gcf()
    return staInfo
end
