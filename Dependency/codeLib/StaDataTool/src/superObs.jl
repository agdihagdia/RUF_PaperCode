using StaDataM
export mergeStaByGrid
function mergeStaByGrid(staObj,lonM,latM;weightedAverage=false)
    # merge the stations in the same grid
    # if weightedAverage is true, you should provide staObj.staInfo.weight
    staInfo_ori=staObj.staInfo
    inds=getNearest(lonM,latM,staInfo_ori.lon,staInfo_ori.lat;outType="linearind")
    mergedStaID_all=string.("merged_",inds)
    oriStaID=staInfo_ori.staID
    return mergeStation(staObj, oriStaID, mergedStaID_all;weightedAverage=weightedAverage)
end