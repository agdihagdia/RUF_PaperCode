using QualityControlM,StaDataM
export qcStaData

function qcStaData(staObs;analysis=false,ind2ref_all=axes(staObs.staInfo,1))
    staObsAll_ori=staObs.staDataAll
    staInfo=staObs.staInfo
    varNames=staObs.varNames 
    leadTimes=[0]
    if analysis == false
        staObsAll_ori=staObs.staDataAll
        ind2ref_all=findall(in.(staInfo.type,[["B","C"]]))
        staObsAll_qc=QCraw_v2(staObsAll_ori,staInfo,varNames;probabilityLimit=1e-10,analysis=false,disp=false,ind2ref_all=ind2ref_all);
        return staDATA(staObs.times,varNames,staInfo,leadTimes,staObsAll_qc)
    else
        staObsAll_qc,analysis=staObsQualityControl(staObsAll_ori,staInfo,varNames;analysis=analysis,ind2ref_all=ind2ref_all)
        return staDATA(staObs.times,varNames,staInfo,leadTimes,staObsAll_qc),analysis
    end
end