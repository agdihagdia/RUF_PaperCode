using StaDataM,Dates,TimeZones,DataFrames,DelimitedFiles,CSV,utilities

export txt2staData,readAObsTxt
# updates:
# 2022-12-17 add support for disabling multi-threads in readAObsTxt to enable multi-threads in txt2staData
# Sample Usage:
#     # read a file -----------
#     fileName="20211129180000.txt";
#     staObs=txt2staData(fileName)

#     # read multiple files -----------
#     fileFun_obs(t)="/public/home/tangxiao/data_tangxiao/konglei/20DXL/OI_3km/obsinput/"*Dates.format(astimezone(t,tz"Asia/Shanghai"),"yyyymmddHH0000.txt")
#     t=ZonedDateTime(2022,9,28,tz"Asia/Shanghai")
#     fileName=fileFun_obs(t)
#     readAObsTxt(fileName)
#     staObs=txt2staData(fileFun_obs(t))
#     staD=txt2staData(fileFun_obs.(t .+ Hour.([-2*24:4*24;])),staObs.staInfo;disp=true)

function txt2staData(fileName::String;tzInFile=tz"Asia/Shanghai")
    
    varNames=["pm25","pm10","co","no2","so2","o3"]
    dataInFile=readAObsTxt(fileName)
    staInfo=dataInFile[:,1:3]
    tStr=basename(fileName)[1:10]
    Times=ZonedDateTime(DateTime(tStr,"yyyymmddHH"),tzInFile)
    staDataAll=Matrix(dataInFile[!,varNames])'
    staDataAll=reshape(staDataAll,(1,length(varNames),size(staInfo,1),1))
    staDataAll[staDataAll.<0].=NaN
    leadTimes=[0]
    return staDATA(Times,varNames,staInfo,leadTimes,staDataAll)
end
function txt2staData(fileName::String,staInfo;tzInFile=tz"Asia/Shanghai")    
    if eltype(staInfo.staID) !== String 
        @error "types of staID in staInfo should be abString!!"
        return
    end
    varNames=["pm25","pm10","co","no2","so2","o3"]
    dataInFile=readAObsTxt(fileName)
    tStr=basename(fileName)[1:10]
    Times=ZonedDateTime(DateTime(tStr,"yyyymmddHH"),tzInFile)
    lia,lob=ismember(staInfo.staID,dataInFile.staID)
    staDataAll=fill(NaN,1,length(varNames),length(staInfo.staID),1)
    staDataAll[1,:,lia,1]=Matrix(dataInFile[lob[lia],varNames])'
    staDataAll[staDataAll.<0].=NaN
    leadTimes=[0]
    return staDATA(Times,varNames,staInfo,leadTimes,staDataAll)
end
function txt2staData(fileName_all::Vector{String},staInfo;tzInFile=tz"Asia/Shanghai",disp=false)
    if eltype(staInfo.staID) !== String 
        @error "types of staID in staInfo should be String!!"
        return
    end
    varNames=["pm25","pm10","co","no2","so2","o3"]
    Times=[ZonedDateTime(DateTime(basename(file_curr)[1:10],"yyyymmddHH"),tzInFile) for file_curr in fileName_all]
    staDataAll=fill(NaN,length(Times),length(varNames),size(staInfo,1),1)
    Threads.@threads for n_file in eachindex(fileName_all)
    # for n_file in eachindex(fileName_all)
        file_curr=fileName_all[n_file]
        if !isfile(file_curr)
            if disp;println("File not found: ",file_curr);end
            continue
        end
        if disp;println("Reading ",file_curr);end
        # read a file ------------------------------  
        dataInFile=readAObsTxt(file_curr;par=false)
        # extract wanted data and write to staDataAll -------------
        loc2extract=indexin(staInfo.staID,dataInFile.staID );
        qOut_validSta=.~isnothing.(loc2extract);
        data_aFile=fill(NaN,length(varNames),size(staInfo,1))
        data_aFile[:,qOut_validSta]=Matrix(dataInFile[loc2extract[qOut_validSta],varNames])'
        staDataAll[n_file,:,:,:]=data_aFile
    end
    staDataAll[staDataAll.<0].=NaN
    leadTimes=[0]
    return staDATA(Times,varNames,staInfo,leadTimes,staDataAll)
end
