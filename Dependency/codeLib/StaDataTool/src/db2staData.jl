using DataBaseM,StaDataM

export staDATAFromDB
"""
# sample usage:

Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,6)
varNameS=["PM25","PM10","SO2","NO2","CO","O3"]
staInfo=getStaInfo()
staObs=staDATAFromDB(Times,varNameS,staInfo)
write2mat(staObs,"test.mat")
"""
function staDATAFromDB(Times,varNameS,staInfo;tableName="verified.pollutant_julia")
    staDataAll=getObsPoll(Times,varNameS,staInfo.staID;tableName=tableName)   
    return staDATA(Times,varNameS,staInfo,staDataAll)
end
