using StaDataM

export staDATAFromDB,makeDF2insert
"""
# sample usage:

Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,6)
varNameS=["PM25","PM10","SO2","NO2","CO","O3"]
staInfo=getStaInfo()
staInfo=staInfo[in.(staInfo.type, [["B","C","D"]]),:];
staObs=staDATAFromDB(Times,varNameS,staInfo)
write2mat(staObs,"test.mat")
"""
function staDATAFromDB(Times,varNameS,staInfo;tableName="tc_prod_zhfx.station_qx_data_original",connInfo=getConnInfo())
    if eltype(Times) == DateTime # use Asia/Shanghai as the default timezone
        Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
    else # convert to timezone of Asia/Shanghai
        Times=astimezone.(Times,tz"Asia/Shanghai")
    end
    
    if tableName in ["tc_prod_zhfx.station_qx_data_original","tc_prod_zhfx.station_qx_data_original_mark"]
        staDataAll=getObsPoll(Times,varNameS,staInfo.staID;tableName=tableName,connInfo=connInfo) 
    elseif tableName == "tc_prod_zhfx.t_station_data_qc"
        staDataAll=getQCData(Times,varNameS,staInfo.staID;tableName=tableName,connInfo=connInfo) 
    else
        error("tableName not supported")
    end
    return staDATA(Times,varNameS,staInfo,staDataAll)
end


function makeDF2insert(Times,varNameS,staIDs,obsOri,analysis)
    function setNothing(x)
        x=replace(x,NaN=>nothing)
        x=replace(x,missing=>nothing)
        return x
    end
    obsAll=setNothing(obsOri)
    # convert nan to missing for pdf_time, pdf_space, pdf_stag, pdf_cali, pdf_course, pdf_union, pdf_st in analysis
    for key in ["pdf_time","pdf_space","pdf_stag","pdf_cali","pdf_course","pdf_st"]
        analysis[key]=setNothing(analysis[key])
    end

    df_all=Array{Any,2}(undef, length(Times),length(varNameS))

    pdf_range=fill(1.0,size(obsOri))
    pdf_range[analysis["q_outRange"]] .= 0.0
    for (n_t,t_curr) in enumerate(Times)
        # prepare df to write in a time
        # =========== prepare sql statements =================
        for n_poll in eachindex(varNameS)
            q_validSta=.!isnothing.(obsAll[n_t,n_poll,:])
            if sum(q_validSta) > 0
                varName= varNameS[n_poll]
                validNum=sum(q_validSta)
                df_all[n_t,n_poll]=DataFrame(data_time=fill(t_curr,validNum),station_code=staIDs[q_validSta],element_type=fill(varName,validNum),
                    value=obsAll[n_t,n_poll,q_validSta], pdf_time=analysis["pdf_time"][n_t,n_poll,q_validSta] , pdf_space = analysis["pdf_space"][n_t,n_poll,q_validSta], 
                    pdf_stag= analysis["pdf_stag"][n_t,n_poll,q_validSta],  pdf_cali= analysis["pdf_cali"][n_t,n_poll,q_validSta], 
                    pdf_course= analysis["pdf_course"][n_t,n_poll,q_validSta],     pdf_union= analysis["pdf_st"][n_t,n_poll,q_validSta],
                    isOutRange= analysis["q_outRange"][n_t,n_poll,q_validSta], isValid= .!analysis["q_outlier"][n_t,n_poll,q_validSta])
                # add isPM25OutweightPM10 if "pm10" in lowercase.(varNameS)
                if "pm10" in lowercase.(varNameS) 
                    if lowercase.(varName) == "pm10"
                        isPM25OutweightPM10=analysis["q_PM25OutweightPM10"][n_t,q_validSta]
                    else
                        isPM25OutweightPM10=fill(nothing,sum(q_validSta))
                    end
                    df_all[n_t,n_poll][!,:isPM25OutweightPM10] = isPM25OutweightPM10
                end
            end
        end
    end
    return vcat(df_all...)
end
