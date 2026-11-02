module StaDataM
using Base
export staDATAStru,staDATA,write2mat,toSimForm,toForeForm,getValidBySta,extract,getStaInfo,mergeStation,staread
export mergeCity,mergeStaInfo_ByCity
using MAT, TimeZones,Dates,DataFrames,Statistics,NaNStatistics

# 2022-12-20 add a read method to read forecast data by tStart and Times:
    # Base.read(obj::staDATAStru,tStart=obj.times[1],Times=tStart.+Hour.(obj.leadTimes); varNames=obj.varNames,
    #                             staIDs=obj.staInfo.staID)
# 2023-2-16 add support for staID of type SentinelArrays.ChainedVector by replace Vector by AbstractVector in in2vec
# 2023-2-16 add support extract method with staID not found ing obj.staInfo.staID
# 2023-2-16 add read write support staInfo with columns more than staID,lon,lat
# 2023-2-16 force staID to string string in staDATA(fileName_all::Vector{String}
mutable struct staDATAStru
    times::Vector
    varNames::Vector
    staInfo::DataFrame
    leadTimes::Vector
    staDataAll::Array
end
# function toStandardForm_staInfo(staInfo)    
#     staInfo.staID=string.(staInfo.staID)
#     for col in names(staInfo)
#         type=eltype(staInfo[1,col])
#         if type <: AbstractChar; 
#             staInfo[!,col]=string.(strip.(staInfo[!,col]))
#         else
#             if type <: AbstractString; type=String;end
#             staInfo[!,col]=convert(Vector{type},staInfo[!,col])
#         end
#     end
#     return staInfo
# end
function toStandardForm_staInfo(staInfo)    
    # 1. 确保 staID 列是字符串类型
    staInfo.staID = string.(staInfo.staID)
    
    # 2. 遍历每一列
    for col in names(staInfo)
        # 获取该列的元素类型（比取第一行更准确，防止第一行有 missing 或特殊字符）
        col_type = eltype(staInfo[!, col])
        
        if col_type <: AbstractChar || col_type <: AbstractString
            # 如果是字符或字符串，统一去除两端空格并转为标准 String
            # skipmissing 可以防止列里有 missing 值时报错
            staInfo[!, col] = [passmissing(strip)(x) for x in staInfo[!, col]]
            # 或者如果确保没有 missing，直接用：
            # staInfo[!, col] = string.(strip.(staInfo[!, col]))
            
        else
            # 如果是其他类型（如 Int, Float），确保类型转换正确
            staInfo[!, col] = convert(Vector{col_type}, staInfo[!, col])
        end
    end
    
    return staInfo
end
function staDATA(fileName::String)
    if !isfile(fileName);@error "fill not exit: $fileName";end
    vars = matread(fileName)
    staInfo=DataFrame(vars["staInfo"])
    staInfo=toStandardForm_staInfo(staInfo)
    # staInfo.staID=string.(staInfo.staID)
    # for col in names(staInfo)
    #     type=eltype(staInfo[1,col])
    #     if type <: AbstractChar; 
    #         staInfo[!,col]=string.(staInfo[!,col])
    #     else
    #         if type <: AbstractString; type=String;end
    #         staInfo[!,col]=convert(Vector{type},staInfo[!,col])
    #     end
    # end
    # staInfo=toDataFrame_staInfo(staInfo)
    staDataAll=vars["staDataAll"]
    varName=string.(uppercase.(in2vec(vars["varName"])))
    leadTimes=vars["leadTimes"]
    if leadTimes isa Real;leadTimes=[leadTimes];end
    time=ZonedDateTime.(in2vec(vars["time"]), "yyyy-mm-ddTHH:MM:SSzzzz")

    return staDATAStru(time,varName,staInfo,leadTimes,staDataAll)
end
function write2mat(obj::staDATAStru,fileName)
    # file= matopen(fileName, "w")
    matopen(fileName, "w") do file
        # replace missing with NaN or empty string for staInfo
        staInfo2write=copy(obj.staInfo)
        for col in names(staInfo2write)
            if eltype(staInfo2write[!,col]) <: Union{AbstractString,Missing}
                staInfo2write[!,col]=convert.(String,replace!(staInfo2write[!,col],missing=>""))
            elseif eltype(staInfo2write[!,col]) <: Union{AbstractFloat,Missing}
                staInfo2write[!,col]=convert.(Float64,replace!(staInfo2write[!,col],missing=>NaN))
            end
        end
        staInfo2write=Dict(col=>staInfo2write[:,col] for col in names(staInfo2write))   
        write(file, "leadTimes", obj.leadTimes)
        write(file, "staDataAll", obj.staDataAll)
        write(file, "varName", obj.varNames)
        write(file, "staInfo", staInfo2write)
        tString=Dates.format.(obj.times, "yyyy-mm-ddTHH:MM:SSzzzz")
        write(file, "time", tString)
    end
    # close(file)
end
# function toDataFrame_staInfo(staInfo)
#     for key in keys(staInfo)
#         staInfo[key]=vec(staInfo[key])
#     end
#     staInfo["staID"]=string.(staInfo["staID"])
#     staInfo=DataFrame(staInfo)
#     return staInfo
# end
function staread(fileName_all;Times=nothing,varNames=nothing,staInfo=nothing,leadTimes=nothing,disp=false)    
    fileName_all=in2vec(fileName_all)  
    staObj=staDATA(fileName_all;Times=Times,varNames=varNames,staInfo=staInfo,leadTimes=leadTimes,disp=disp) 
    return staObj
end
function staDATA(fileName_all::Vector{String};Times=nothing,varNames=nothing,staInfo=nothing,leadTimes=nothing,disp=false)    
    Times,varNames,leadTimes=in2vec.([Times,varNames,leadTimes])
    local staDataAll
    for n_file=1:length(fileName_all) 
        fileName=fileName_all[n_file]
        if !isfile(fileName)
            if disp;println("File not found: ",fileName);end
            continue
        end
        if disp;print("Reading ",fileName,"  ");end
        vars_afile = matread(fileName)
        # ---------- prepare varName_curr leadTimes_curr time_curr
        varName_curr=string.(uppercase.(in2vec(vars_afile["varName"])))
        leadTimes_curr=in2vec(vars_afile["leadTimes"])
        time_curr=ZonedDateTime.(vars_afile["time"], "yyyy-mm-ddTHH:MM:SSzzzz")
        staInfo_curr=DataFrame(vars_afile["staInfo"])
        staInfo_curr=toStandardForm_staInfo(staInfo_curr)

        # staInfo_curr.staID=string.(staInfo_curr.staID)
        # staInfo.staID=string.(staInfo.staID)
        # staInfo_curr=toDataFrame_staInfo(staInfo_curr)
        
        # ---------- prepare varNames leadTimes time if not specified by inputs
        if ! @isdefined staDataAll
            if isnothing(Times);Times=time_curr;end
            if isnothing(varNames);varNames=varName_curr;end
            if isnothing(staInfo);staInfo=staInfo_curr;end
            if isnothing(leadTimes);leadTimes=leadTimes_curr;end 
            varNames=uppercase.(varNames);
            leadTimes=vec(collect(leadTimes))
            staDataAll=fill(NaN,length(Times),length(varNames),size(staInfo,1),length(leadTimes))
        end
        # ----------- extract data       
        LiaTime,LocTime = ismember(Times,time_curr);
        LiavarName,LocvarName = ismember(varNames,varName_curr);
        LiaStaID,LocStaID = ismember(staInfo.staID,staInfo_curr.staID);
        LiaLead,LocLead = ismember(leadTimes,leadTimes_curr);
        # output
        if disp; println("tvsl",[sum(LiaTime) ,sum(LiavarName), sum(LiaStaID) , sum(LiaLead)]); end
        if sum(LiaTime) >0 && sum(LiavarName) >0 && sum(LiaStaID) >0 && sum(LiaLead) >0 
            staDataAll[LiaTime,LiavarName,LiaStaID,LiaLead] = vars_afile["staDataAll"][
                LocTime[LiaTime],LocvarName[LiavarName],LocStaID[LiaStaID],LocLead[LiaLead]
                ]
        end
    end
    return staDATAStru(Times,varNames,staInfo,leadTimes,staDataAll)
end
function staDATA(Times,varNames,staInfo,staDataAll)
    return staDATA(Times,varNames,staInfo,[0],staDataAll)
end
function staDATA(Times,varNames,staInfo,leadTimes,staDataAll)    
    Times,varNames,leadTimes=in2vec.([Times,varNames,leadTimes])

    return staDATAStru(Times,varNames,staInfo,leadTimes,staDataAll)
end
function ismember(a,b)
    if eltype(a) == eltype(b)
    elseif eltype(a)  <: AbstractString && eltype(b)  <: AbstractString
    elseif eltype(a)  <: Number && eltype(b)  <: Number
    else
        error("element types not comparable ($(a[1]),$(b[1])) ($(eltype(a)), $(eltype(b)))")
    end
    locb=indexin(a,b);
    lia=.~isnothing.(locb);
    return lia,locb
end
function in2vec(in)
    if  in isa StepRange;       in=collect(in);        end
    if !(in isa AbstractVector) && !isnothing(in);      in=[in];          end
    return in
end
function Base.read(obj::staDATAStru,tStart::ZonedDateTime;Times=tStart.+Hour.(obj.leadTimes), varNames=obj.varNames,
                            staIDs=obj.staInfo.staID)
    leads2extract=[Hour(x).value for x in Times-tStart]
    return Base.read(obj;Times=tStart,varNames=varNames,staIDs=staIDs, leadTimes=leads2extract)
end
function Base.read(obj::staDATAStru;Times=obj.times, varNames=obj.varNames,
                            staIDs=obj.staInfo.staID, leadTimes=obj.leadTimes)
    # find position
    Times,varNames,leadTimes,staIDs=in2vec.([Times,varNames,leadTimes,staIDs])
    
    LiaTime,LocTime = ismember(Times,obj.times);
    LiavarName,LocvarName = ismember(uppercase.(varNames),uppercase.(obj.varNames));
    LiaStaID,LocStaID = ismember(staIDs,obj.staInfo.staID);
    LiaLead,LocLead = ismember(leadTimes,obj.leadTimes);
    
    # output
    outData=fill(NaN,(length(Times),length(varNames),length(staIDs),length(leadTimes)))
    outData[LiaTime,LiavarName,LiaStaID,LiaLead]=obj.staDataAll[
        LocTime[LiaTime],LocvarName[LiavarName],LocStaID[LiaStaID],LocLead[LiaLead]
        ]
    return outData
end
function Base.put!(obj::staDATAStru,data2write;Times=obj.times, varNames=obj.varNames,
                            staIDs=obj.staInfo.staID, leadTimes=obj.leadTimes)
    Times,varNames,leadTimes,staIDs=in2vec.([Times,varNames,leadTimes,staIDs])
    
    # find position
    LiaTime,LocTime = ismember(Times,obj.times);
    LiavarName,LocvarName = ismember(uppercase.(varNames),uppercase.(obj.varNames));
    LiaStaID,LocStaID = ismember(staIDs,obj.staInfo.staID);
    LiaLead,LocLead = ismember(leadTimes,obj.leadTimes);
    
    # output
    obj.staDataAll[
        LocTime[LiaTime],LocvarName[LiavarName],LocStaID[LiaStaID],LocLead[LiaLead]
        ]=data2write[LiaTime,LiavarName,LiaStaID,LiaLead]
end
function extract(obj::staDATAStru;Times=obj.times, varNames=obj.varNames,
                            staIDs=obj.staInfo.staID, leadTimes=obj.leadTimes)
    Times,varNames,staIDs=in2vec.([Times,varNames,staIDs])
    staInfo=getStaInfo(obj,staIDs)
    staDataAll=Base.read(obj;Times=Times, varNames=varNames,staIDs=staInfo.staID,leadTimes=leadTimes)   
    return staDATAStru(Times,varNames,staInfo,leadTimes,staDataAll)
end
function getStaInfo(obj::staDATAStru,staIDs)
    staIDs=in2vec(staIDs)
    ind_out=indexin(staIDs,obj.staInfo.staID)
    q_ismissing=isnothing.(ind_out)
    if sum(q_ismissing)>0
        @warn "stations not found, omiting: $(staIDs[q_ismissing])"
    end
    return obj.staInfo[ind_out[.!q_ismissing],:]
end
# convert form foreform(tStart,var,sta,lead) to sim form(times,var,sta), filling the gaps between tStarts with lead time specified
function toSimForm(obj::staDATAStru,startLead::Real;Times=obj.times, varNames=obj.varNames,
    staIDs=obj.staInfo.staID)
    Times,varNames,staIDs=in2vec.([Times,varNames,staIDs])
    staInfo_out=getStaInfo(obj,staIDs)
    leads2extract=startLead.+[1:(Hour(Times[2]-Times[1]).value);] .-1
    Times_out=[tStart+Hour(lead) for tStart in Times for lead in leads2extract]
    data_ori=Base.read(obj;Times=Times, varNames=varNames,staIDs=staIDs,leadTimes=leads2extract)
    staDataAll=reshape(
                    permutedims(data_ori,(4,1,2,3))
                    ,length(Times_out),length(varNames),length(staIDs),1)
    
    return staDATAStru(Times_out,varNames,staInfo_out,[startLead],staDataAll)
end

# convert form foreform(tStart,var,sta,lead) to sim form(times,var,sta,startLead), filling the gaps between tStarts with lead time specified
function toSimForm(obj::staDATAStru,startLead_all::AbstractVector;Times=obj.times, varNames=obj.varNames,
        staIDs=obj.staInfo.staID)
    staInfo=obj.staInfo
    staSim_all=collect(toSimForm(obj,startLead;varNames=varNames,staIDs=staIDs,Times=Times) for startLead in startLead_all)
    Times_out=sort(unique(vcat(collect(s.times for s in staSim_all)...)))
    staDataAll=fill(NaN,length(Times_out),length(varNames),length(staInfo.staID),length(startLead_all))
    for (nStart,startLead) in enumerate(startLead_all)
        staDataAll[:,:,:,nStart] = Base.read(staSim_all[nStart];Times=Times_out)
    end
    staSim=staDATA(Times_out,varNames,staInfo,startLead_all,staDataAll)
    return staSim
end
# convert form foreform(tStart,var,sta,lead) to sim form(times,var,sta) for single tStart
function toSimForm(obj::staDATAStru;Times=obj.times, varNames=obj.varNames,
    staIDs=obj.staInfo.staID, leadTimes=obj.leadTimes)
    if length(Times) > 1
        @error "multiple tStarts found, using other method instead: toSimForm(obj,startLead)"
    end
    Times,varNames,staIDs=in2vec.([Times,varNames,staIDs])
    staInfo_out=getStaInfo(obj,staIDs)
    Times_out=Times[1] .+ Hour.(obj.leadTimes)
    data_ori=Base.read(obj;Times=Times, varNames=varNames,staIDs=staIDs,leadTimes=leadTimes)
    staDataAll= permutedims(data_ori,(4,2,3,1))
    
    return staDATAStru(Times_out,varNames,staInfo_out,[0],staDataAll)
end
function getValidBySta(obj::staDATAStru;Times=obj.times, varNames=obj.varNames,
            staIDs=obj.staInfo.staID, leadTimes=obj.leadTimes)
    data_tmp=Base.read(obj;Times=Times, varNames=varNames,staIDs=staIDs,leadTimes=leadTimes)
    validNum=sum(.!isnan.(data_tmp),dims=(1,2,4))[1,1,:,1]
    return validNum
end
function toForeForm(obj::staDATAStru)
    if length(obj.leadTimes) > 1
        staFore=toForeForm_multiLead(obj::staDATAStru,obj.times,obj.leadTimes)
    else
        staFore=toForeForm_singleLead(obj::staDATAStru,obj.times[1],[0:length(obj.times)-1;])
    end
    return staFore
end

function toForeForm(obj::staDATAStru,tStart_all,leads_all)
    if length(obj.leadTimes) > 1
        staFore=toForeForm_multiLead(obj::staDATAStru,tStart_all,leads_all)
    else
        staFore=toForeForm_singleLead(obj::staDATAStru,tStart_all,leads_all)
    end
    return staFore
end
# ------------------ convert SimForm to ForeForm for multiple lead time in SimFore
function toForeForm_multiLead(obj::staDATAStru,tStart_all,leads_all)
    tStart_all,leads_all=in2vec.([tStart_all,leads_all])
    staDataAll=fill(NaN,length(tStart_all),length(obj.varNames),size(obj.staInfo,1),length(leads_all))
    for (n_start,tStart) in enumerate(tStart_all)
        for (n_lead_inFore,lead) in enumerate(leads_all)
            t_curr=tStart .+ Hour.(lead)
            idx_t=indexin([t_curr], obj.times)[1]
            n_lead_inSim=indexin([lead],obj.leadTimes)[1]
            if !isnothing(idx_t) && !isnothing(n_lead_inSim)
                staDataAll[n_start,:,:,n_lead_inFore] =  obj.staDataAll[idx_t,:,:,n_lead_inSim]
            end
        end
    end
    return staDATAStru(tStart_all,obj.varNames,obj.staInfo,leads_all,staDataAll)
end

# ------------------ convert SimForm to ForeForm for single lead time in SimFore
function toForeForm_singleLead(obj::staDATAStru,tStart_all,leads_all)
    tStart_all,leads_all=in2vec.([tStart_all,leads_all])
    staDataAll=fill(NaN,length(tStart_all),length(obj.varNames),size(obj.staInfo,1),length(leads_all))
    for (n_start,tStart) in enumerate(tStart_all)
        Times_tmp=tStart .+ Hour.(leads_all)
        LiaTime,LocTime = ismember(Times_tmp,obj.times);
        staDataAll[n_start,:,:,LiaTime]=permutedims(
                                    obj.staDataAll[LocTime[LiaTime],:,:,:],
                                    (4,2,3,1)
                                    )
    end
    return staDATAStru(tStart_all,obj.varNames,obj.staInfo,leads_all,staDataAll)
end
# toForeForm(obj::staDATAStru)=toForeForm_singleLead(obj,obj.times[1],[0:length(obj.times)-1;])

function mergeStation(obj::staDATAStru, oriStaID::AbstractVector, mergedStaID_all::AbstractVector;weightedAverage=weightedAverage)
    # mergge stations by staIDs
    # the length of oriStaID and mergedStaID_all should be the same, same staID in mergedStaID_all will be merged
    # if weightedAverage is true, you should provide obj.staInfo.weight

    oriStaInfo = obj.staInfo
    if weightedAverage
        if !("weight" in names(oriStaInfo))
            @error "When weightedAverage is true, you should provide staObj.staInfo.weight"
        else
            oriStaWeight=oriStaInfo.weight
        end
    else
        oriStaWeight=ones(length(oriStaInfo.staID))
    end
    # exclude stations not in oriStaInfo
    q_validOriSta = in.(oriStaID, [oriStaInfo.staID])
    oriStaID = oriStaID[q_validOriSta]
    mergedStaID_all = mergedStaID_all[q_validOriSta]
    # initialize mergedStaInfo
    mergedStaID_unique = unique(mergedStaID_all)
    mergedStaInfo = DataFrame(staID = mergedStaID_unique, lon = fill(NaN, length(mergedStaID_unique)), lat = fill(NaN, length(mergedStaID_unique)), 
                        oriStaIDs = fill([], length(mergedStaID_unique)),weight= fill(NaN, length(mergedStaID_unique)))
    # initialize mergedStaDataAll
    sizeMegStaDataAll = collect(size(obj.staDataAll))
    sizeMegStaDataAll[3] = length(mergedStaInfo.staID)
    mergedStaDataAll = fill(NaN, sizeMegStaDataAll...)
    # construct mergedStaInfo and mergedStaDataAll
    for (n_sta,staID_curr) in enumerate(mergedStaInfo.staID)
        q2merge_aSta = isequal.(staID_curr, mergedStaID_all)
        weights_curr=oriStaWeight[q2merge_aSta] ./ sum(oriStaWeight[q2merge_aSta])
        mergedStaInfo.lon[n_sta] = sum(oriStaInfo.lon[q2merge_aSta] .* weights_curr)
        mergedStaInfo.lat[n_sta] = sum(oriStaInfo.lat[q2merge_aSta] .* weights_curr)
        mergedStaDataAll[:, :, n_sta, :] = sum(obj.staDataAll[:, :, q2merge_aSta, :] .* reshape(weights_curr,1,1,:,1), dims=3)        
        mergedStaInfo.weight[n_sta] = sum(oriStaWeight[q2merge_aSta])
        mergedStaInfo.oriStaIDs[n_sta] = oriStaInfo.staID[q2merge_aSta]
    end
    return staDATAStru(obj.times,obj.varNames,mergedStaInfo,obj.leadTimes,mergedStaDataAll)
end

function mergeCity(staObj::staDATAStru;Times=staObj.times,varNames=staObj.varNames,staInfo=staObj.staInfo,city_all=unique(staInfo.city),leads_all=staObj.leadTimes)
    if !(varNames isa AbstractVector);varNames=[varNames];end
    cityInfo=mergeStaInfo_ByCity(staInfo;city_all=city_all)
    # exclude 对照点
    # q_valid=.!occursin.("(对照点)",staInfo.staName)
    # staInfo=staInfo[q_valid,:]
    staDataAll=fill(NaN,length(Times),length(varNames),length(city_all),length(leads_all))
    for (n_city,city) in enumerate(city_all)
        staIDs_aCity=staInfo.staID[isequal.(staInfo.city,city)]
        staDataAll[:, :, n_city:n_city,:]=nanmean(read(staObj;Times=Times,varNames=varNames,staIDs=staIDs_aCity,leadTimes=leads_all); dims=3)
    end
    cityData=staDATA(Times,varNames,cityInfo,leads_all,staDataAll)
    return cityData
end
function mergeStaInfo_ByCity(staInfo;city_all=unique(staInfo.city))
    q_invalid=isempty.(city_all)
    city_all=string.(city_all[.!q_invalid])
    cityInfo=DataFrame(lon=NaN,lat=NaN,staID=city_all)
    for (n_city,city) in enumerate(city_all)
        q_currCity=isequal.(staInfo.city,city)
        cityInfo.lon[n_city]=nanmean(staInfo.lon[q_currCity])
        cityInfo.lat[n_city]=nanmean(staInfo.lat[q_currCity])
    end
    # cityInfo.staID=unifyCityName(cityInfo.staID)
    return cityInfo
end

end# module
