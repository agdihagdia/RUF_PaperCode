
using PyCall,MyPlotM
export animateSandStorm
function animateSandStorm(data2plot::AbstractArray,Times,lonM,latM;cLimit=[0,2000],cMap="afmhot_r",dpi=100,interval=200,
    outFileName="dust storm started at $(Dates.format(Times[1],"yyyy-mm-dd HH")).mp4")
    
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
        h_bar.ax.set_title("   " * raw"μg m$^{-3}$")
        plt.axis([72, 137, 17.5, 55])
        plt.title("沙尘  $(Dates.format(Times[n],"yyyy-mm-dd HH"))")
        hp_all=[hp,hp]
        return hp_all
    end
    
    ani = animation.FuncAnimation(fig, plot_aframe, frames=length(Times), interval=interval, blit=true)
    
    ani.save(outFileName)

    return ani
end
    