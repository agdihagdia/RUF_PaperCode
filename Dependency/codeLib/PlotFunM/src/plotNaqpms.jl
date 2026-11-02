using GradDataM,MyPlotM,PyCall,Dates,TimeZones,NaNStatistics
export animateNaqpms,animateNaqpmsWithObs,animateNaqpmsWithWind
"""
sample usage:

grd=gradData("lonlatM_d2.ctl";disp=false)
lonM=read(grd,"lonM",1)
latM=read(grd,"latM",1)

fig = plt.figure(figsize=(14,4))
animateNaqpms(ctlFileName,lonM,latM,varName)
"""
function animateNaqpms(ctlFileName::AbstractString,Times,lonM,latM,varName;nz=1,varType="direct",addSand=true,cLimit=nothing,TimeZoneDisp=nothing,cMap="jet",dpi=100,
    outFileName="animation of $varName started at $(Dates.format(Times[1],"yyyy-mm-dd HH")).mp4")

    # prepare data2plot ==========================
    data2plot=Array{Float32,3}(undef,(size(lonM)...,length(Times)))
    grd=gradData(ctlFileName;disp=false)

    # prepare readFun
    # local readFun::Function
    function readFun(grd)
        if varType == "combine"
            if transUnit
                data=toMass!(combinePollFromGrd(grd,varName,nz;addSand=addSand))
            else
                data=combinePollFromGrd(grd,varName,nz;addSand=addSand)
            end
        elseif varType == "direct"
            data=read(grd,varName,nz)
        else
            error("varType must be either combine or direct")
        end
        return data
    end
    # read data2plot 
    for (n_t,t) in enumerate(Times)
        data2plot[:,:,n_t]=readFun(setTime!(grd,t))
    end
    # plot the animation ==========================
    if isnothing(cLimit); cLimit=[0,nanpctile(data2plot,99.9)]; end
    ani=animateNaqpms(data2plot::AbstractArray,astimezone.(Times,TimeZoneDisp),lonM,latM,varName; cLimit=cLimit, outFileName=outFileName,cMap=cMap,dpi=dpi)
    return ani
end

function animateNaqpms(data2plot::AbstractArray,Times,lonM,latM,varName;cLimit=[0,nanpctile(data2plot,99.9)],cMap="jet",dpi=100,
    outFileName="animation of $varName started at $(Dates.format(Times[1],"yyyy-mm-dd HH")).mp4")
    
    # plot ======================================
    
    plt=get_plt()
    animation=pyimport("matplotlib.animation")
    # fig = plt.figure(figsize=(6,5))
    # gs = fig.add_gridspec(nrows=1, ncols=1,left=0.05,right=0.98,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)

    fig = plt.figure(figsize=(6,4.9),dpi=dpi)
    gs = fig.add_gridspec(1,1,left=0.05,right=0.98,bottom=0.04,top=0.93, hspace=0.1, wspace=0.1)

    function plot_aframe(i)
        n=i+1
        plt.clf()
        ax=plt.subplot(gs[1])
        hp=plt.pcolormesh(lonM,latM,data2plot[:,:,n],cmap=cMap,shading="gouraud")
        plt.clim(cLimit);
        hcn=plotChina()
        # plt.colorbar(extend="both")
        h_bar=addChinaColorbar()
        h_bar.ax.set_title("   " * getUnitStr(varName))
        plt.axis([72, 137, 17.5, 55])
        plt.title("$(varNameDisp(varName)) at $(Dates.format(Times[n],"yyyy-mm-dd HH"))")
        hp_all=[hp,hp]
        return hp_all
    end
    
    ani = animation.FuncAnimation(fig, plot_aframe, frames=length(Times), interval=200, blit=true)
    
    ani.save(outFileName)

    return ani
end
    
function animateNaqpmsWithWind(dataPoll::AbstractArray,dataU::AbstractArray,dataV::AbstractArray,Times,lonM,latM,varName;cLimit=[0,nanpctile(dataPoll,99.9)],cMap="jet",dpi=100,
    outFileName="animation of $varName started at $(Dates.format(Times[1],"yyyy-mm-dd HH")).mp4",windSpace=10,quiverWidth=0.002)
    
    # plot ======================================
    
    plt=get_plt()
    animation=pyimport("matplotlib.animation")
    # fig = plt.figure(figsize=(6,5))
    # gs = fig.add_gridspec(nrows=1, ncols=1,left=0.05,right=0.98,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)

    fig = plt.figure(figsize=(6,4.9),dpi=dpi)
    gs = fig.add_gridspec(1,1,left=0.05,right=0.98,bottom=0.04,top=0.93, hspace=0.1, wspace=0.1)

    function plot_aframe(i)
        n=i+1
        plt.clf()
        ax=plt.subplot(gs[1])
        hp=plt.pcolormesh(lonM,latM,dataPoll[:,:,n],cmap=cMap,shading="gouraud")
        plt.clim(cLimit);
        hcn=plotChina()
        # plt.colorbar(extend="both")
        h_bar=addChinaColorbar()
        h_bar.ax.set_title("   " * getUnitStr(varName))
        # add wind quiver
        plt.quiver(lonM[1:windSpace:end,1:windSpace:end],latM[1:windSpace:end,1:windSpace:end],dataU[1:windSpace:end,1:windSpace:end,n],dataV[1:windSpace:end,1:windSpace:end,n],
            alpha=0.9,width=quiverWidth)#headwidth=3,headlength=3)#,color=[1,1,1] .* 0.8)
        plt.axis([72, 137, 17.5, 55])
        plt.title("$(varNameDisp(varName)) at $(Dates.format(Times[n],"yyyy-mm-dd HH"))")
        hp_all=[hp,hp]
        return hp_all
    end
    
    ani = animation.FuncAnimation(fig, plot_aframe, frames=length(Times), interval=200, blit=true)
    
    ani.save(outFileName)

    return ani
end
    
function animateNaqpmsWithObs(ctlFileName::AbstractString,staObs,Times,lonM,latM,varName;nz=1,varType="direct",addSand=true,cLimit=nothing,
    outFileName="animation of $varName started at $(Dates.format(Times[1],"yyyy-mm-dd HH")).mp4")

    # prepare data2plot ==========================
    data_mesh2plot=Array{Float32,3}(undef,(size(lonM)...,length(Times)));
    grd=gradData(ctlFileName;disp=false)
    staInfo=staObs.staInfo
    data_staObs=dropdims(read(staObs;Times=Times,varNames=varName,staIDs=staInfo.staID),dims=2)

    # prepare readFun
    # local readFun::Function
    function readFun(grd)
        if varType == "combine"
            if transUnit
                data=toMass!(combinePollFromGrd(grd,varName,nz;addSand=addSand))
            else
                data=combinePollFromGrd(grd,varName,nz;addSand=addSand)
            end
        elseif varType == "direct"
            data=read(grd,varName,nz)
        else
            error("varType must be either combine or direct")
        end
        return data
    end
    # read data2plot 
    for (n_t,t) in enumerate(Times)
        data_mesh2plot[:,:,n_t]=readFun(setTime!(grd,t))
    end
    # plot the animation ==========================
    if isnothing(cLimit); cLimit=[0,nanpctile(data_mesh2plot,99.9)]; end
    plt=get_plt()
    animation=pyimport("matplotlib.animation")
    fig = plt.figure(figsize=(6,5))
    gs = fig.add_gridspec(nrows=1, ncols=1,left=0.05,right=0.98,bottom=0.05,top=0.93, hspace=0.1, wspace=0.1)
    function plot_aframe(i)
        n=i+1
        plt.clf()
        ax=plt.subplot(gs[1])
        hp=plt.pcolormesh(lonM,latM,data_mesh2plot[:,:,n],cmap="jet",shading="gouraud")
        plt.clim(cLimit);
        plt.scatter(staInfo.lon,staInfo.lat,s=1,c=data_staObs[n,:],cmap="jet")
        plt.clim(cLimit);
        hcn=plotChina()
        plt.axis([72, 137, 17.5, 55])
        plt.colorbar()

        plt.title("$varName at $(Dates.format(Times[n],"yyyy-mm-dd HH"))")
        hp_all=[hp,hp]
        return hp_all
    end
    
    ani = animation.FuncAnimation(fig, plot_aframe, frames=length(Times), interval=200, blit=true)
    
    ani.save(outFileName)

    return ani
end
