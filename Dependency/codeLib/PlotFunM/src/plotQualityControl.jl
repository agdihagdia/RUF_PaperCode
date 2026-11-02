using StaDataM,MyPlotM,LocalizerM
export plotQC_aSta,getValidStaInfo_nearby
# function plotQC_aSta(Times,obs,q_outlier;estimate=nothing,pdf_misfit=nothing,ax=nothing)    

#     # create figure if needed
#     if isnothing(ax)
#         fig, ax1 = plt.subplots(figsize=(12,4))
#     else
#         ax1=ax
#     end

#     # plot observation 
#     plt.plot(Times,obs,linewidth=3,label="original")
#     # plot estimation
#     if !isnothing(estimate)
#         plt.plot(Times,estimate,linewidth=1,label="estimate")
#     end
#     # plot outlier
#     plt.plot(Times[q_outlier],obs[q_outlier],"r*",label="outlier")

#     plt.legend()
#     plt.xlim(Times[[1,end]])
    
#     # plot pdf
#     if !isnothing(pdf_misfit)
#         ax1.patch.set_visible(false)
#         ax2 = ax1.twinx()
#         ax2.set_zorder(ax1.get_zorder()-1)
#         ax2.plot(Times,pdf_misfit,"g.",label="pdf")
#         ax2.set_ylabel("pdf", color="g")
#         ax2.set_yscale("log")
#         ax2.set_xlim(Times[[1,end]])

#         return ax1,ax2
#     end
#     return ax1
# end

function getValidStaInfo_nearby(staID,staInfo_ref;indNear=1:11)
    lonlat=Matrix(staInfo_ref[:,["lon","lat"]])
    locM=localizerM(lonlat,lonlat)
    ind_target=indexin([staID],staInfo_ref.staID)
    staInfo_ref=staInfo_ref[LocalizerM.getNearest(locM,ind_target[1],indNear),:]
    return staInfo_ref
end
function plotQC_aSta(staOri,analysis,Times,varName,staID;staInfo_ref=staOri.staInfo)  
    
    # --------------------- get data to plot
    staInfo_target=StaDataM.getStaInfo(staOri,staID)
    # get the staInfo of nearby sites
    staInfo_nearby=getValidStaInfo_nearby(staID,staInfo_ref)
    # get the original data
    obsTarget=vec(read(staOri;Times=Times,varNames=varName,staIDs=staID));
    obsRef=read(staOri;Times=Times,varNames=varName,staIDs=staInfo_nearby.staID);
    obsRef=reshape(obsRef,length(Times),length(staInfo_nearby.staID))
    # get the idx for analysis results
    q_t=in.(Times,[staOri.times])
    idx_v= findfirst(isequal.(varName,staOri.varNames))
    idx_s= findfirst(isequal.(staID,staOri.staInfo.staID))

    # --------------------- plot
    t2plot=DateTime.(Times)
    plt.close("all")
    fig=plt.figure(figsize=(10,5),dpi=200)
    gs = fig.add_gridspec(nrows=2, ncols=1,left=0.055,right=0.96,bottom=0.05,top=0.94, hspace=0.1, wspace=0.1)
    # plot observations
    plt.subplot(gs[1])
    labelNear=fill("_附近站点",length(staInfo_nearby.staID));labelNear[1]="附近站点";
    h_near=plt.plot(t2plot,obsRef,c=[1,1,1].* 0.8,label=labelNear)
    h_before=plt.plot(t2plot,obsTarget,"b-",label="质控前")
    q_outlier=analysis["q_outlier"][q_t,idx_v,idx_s]
    h_outlier=plt.plot(t2plot[q_outlier],obsTarget[q_outlier],"r*",label="算法识别异常值")
    plt.xlim(t2plot[1],t2plot[end])
    plt.legend()
    plt.ylabel( "$(varNameDisp(varName)) ($(getUnitStr(varName)))" ,fontsize=10)
    titleName=staInfo_target.province[1] * " " * staInfo_target.city[1] * " " * staInfo_target.staName[1] * " " * staInfo_target.staID[1] * " " * varNameDisp(varName) 
    plt.title(titleName,fontsize=12)
    # plot analysis results
    plt.subplot(gs[2])
    plt.plot(t2plot,analysis["pdf_st"][q_t,idx_v,idx_s],".",label="时空合理性")
    plt.plot(t2plot,analysis["pdf_time"][q_t,idx_v,idx_s],label="时间合理性")
    plt.plot(t2plot,analysis["pdf_space"][q_t,idx_v,idx_s],label="空间合理性")
    plt.plot(t2plot,analysis["pdf_stag"][q_t,idx_v,idx_s],label="时间梯度合理性")
    plt.plot(t2plot,analysis["pdf_cali"][q_t,idx_v,idx_s],label="周期合理性")
    plt.plot(t2plot,analysis["pdf_course"][q_t,idx_v,idx_s],label="前后相邻合理性")
    if isequal(uppercase.(varName),"PM10")
        val_outweight=analysis["q_PM25OutweightPM10"][q_t,idx_s] .* 1e-8
        plt.plot(t2plot,val_outweight,"*",label="倒挂异常")
    end
    plt.xlim(t2plot[1],t2plot[end])
    plt.yscale("log")
    plt.ylabel("合理性指数（越大越合理）",fontsize=10)
    plt.legend()
    return staInfo_nearby
end
