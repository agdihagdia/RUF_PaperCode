module ZZDBM

using MySQL,DataFrames,Dates,DBInterface,TimeZones
export getObsPoll,getStaInfo,getValidStaID,getValidStaInfo,query2df,insertDF,getConn,getQCData,getQCAnalysis,getObsMet,getConnInfo
import utilities:ismember,in2vec

include("additional.jl")
"""

updates:
2023-04-24: use ZonedDateTime in getObsPoll

"""
mutable struct connInfoStr
    host::String
    user::String
    passwd::String
    port::Int
end
function getConnInfo(configName::AbstractString)
    if configName=="local"
        return connInfoStr(ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],39030)
    elseif configName=="direct"
        return connInfoStr(ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],9030)
    elseif configName=="n1"
        return connInfoStr(ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],9030)
    elseif configName=="n2"
        return connInfoStr(ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],9030)
    elseif configName=="n3"
        return connInfoStr(ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],9030)
    else
        error("configName should be either local or direct")
    end
end
# function getConnInfo(;host=ENV["RUF_DB_HOST"],user=ENV["RUF_DB_USER"],passwd=ENV["RUF_DB_PASSWORD"],port=39030)
function getConnInfo(;host=ENV["RUF_DB_HOST"],user=ENV["RUF_DB_USER"],passwd=ENV["RUF_DB_PASSWORD"],port=39030)
    return connInfoStr(host,user,passwd,port)
end

function getConn(conInfo::connInfoStr)
    return DBInterface.connect(MySQL.Connection,conInfo.host,conInfo.user,conInfo.passwd,port=conInfo.port,write_timeout=3600,read_timeout=3600,connect_timeout=5)
end


"""

"""
function getConn(source="local")
    # conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],port=9030)
    if source == "direct"
        conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],port=9030,write_timeout=3600,read_timeout=3600,connect_timeout=5)
    elseif source == "local"
        conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],port=39030,write_timeout=3600,read_timeout=3600,connect_timeout=5)
    else
        error("source should be either direct or local")
    end
    return conn
end


"""
Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,5,2)
varNameS=["PM25","PM10","SO2","NO2","CO","O3"]
staInfo=getValidStaInfo(Times)
@time getObsPoll(Times,varNameS,staInfo.staID)
"""
function getObsPoll(Times,varNameS,staIDs;tableName="tc_prod_zhfx.station_qx_data_original",disp=false,connInfo = getConnInfo())
    conn=getConn(connInfo)
    Times,varNameS,staIDs=in2vec.([Times,varNameS,staIDs])
    if eltype(Times) == DateTime # use Asia/Shanghai as the default timezone
        Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
    else # convert to timezone of Asia/Shanghai
        Times=astimezone.(Times,tz"Asia/Shanghai")
    end

    obsData=fill(NaN,length(Times),length(varNameS),length(staIDs))
    
    if length(staIDs) < 2000
        # read the specified stations
        sqlquery=join(["SELECT data_time as obsTime,station_code as staID,",join(varNameS,","),
            " FROM ",tableName,
            " WHERE duration_code = 6 and (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )",
            " AND  (station_code in ('",join(staIDs,"','"),"') ) "
            ]);
    else
        # read all stations, becase ready by station will be too slow if there are too many stations
        sqlquery=join(["SELECT data_time as obsTime,station_code as staID,",join(varNameS,","),
            " FROM ",tableName,
            " WHERE duration_code = 6 and (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )"
            ]);
    end
    if disp;println("SQL Query: ", sqlquery);end
    response=DBInterface.execute(conn, sqlquery)
    close(conn)
    data_all=DataFrame(response)
    data_all.obsTime=ZonedDateTime.(data_all.obsTime,tz"Asia/Shanghai")

    # extract the wanted data and put into obsData
    # for n_t in eachindex(Times)
    Threads.@threads for n_t in eachindex(Times)
        # if disp;println("processing station: ",n_sta," ",staIDs[n_sta]);end
        ind_all_currT=findall(data_all.obsTime.==Times[n_t])
        if length(ind_all_currT) >0
            q_validInd,indOut=ismember(data_all.staID[ind_all_currT],staIDs)
            data2write=Float64.(replace!(Matrix(data_all[ind_all_currT[q_validInd],varNameS]),missing=>NaN)')
            data2write[data2write .< 0 ] .= NaN
            obsData[n_t,:,indOut[q_validInd]]=data2write
        end
    end
    # for n_t in eachindex(Times)
    #     # Threads.@threads for n_sta in eachindex(staIDs)
    #     # if disp;println("processing station: ",n_sta," ",staIDs[n_sta]);end
    #     t_curr=Times[n_t]
    #     data_aT=
    #     indSQL_all_currSta=findall(data_all.staID.==staIDs[n_sta])
    #     obsData_aSta=fill(NaN,length(Times),length(varNameS))
    #     if length(indSQL_all_currSta) >0
    #         q_validInd,indOut=ismember(data_all.obsTime[indSQL_all_currSta],Times)
    #         obsData_aSta[indOut[q_validInd],:]=dataMatrix_all[indSQL_all_currSta[q_validInd],:]
    #     end
    #     obsData[:,:,n_sta]=obsData_aSta
    # end
    return obsData
end

"""

Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,5,2)
varNameS=["pressure","temperature","humidity","winSpeed","winDir"]
staInfo=getValidStaInfo(Times;tableName="tc_prod_zhfx.weather_data_original")
# staInfo=staInfo[[ occursin("ee", x) for x in staInfo.type],:]
@time getObsMet(Times,varNameS,staInfo.staID)
"""
function getObsMet(Times,varNameS,staIDs;tableName="tc_prod_zhfx.weather_data_original",disp=false,connInfo = getConnInfo())
    conn=getConn(connInfo)
    Times,varNameS,staIDs=in2vec.([Times,varNameS,staIDs])
    if eltype(Times) == DateTime
        Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
    end
    # dic_inside2pub=Dict("PRS"=>"pressure","TEM"=>"temperature","RHU"=>"humidity","WIN_S_INST"=>"winSpeed","WIN_D_INST"=>"winDir")
    dic_pub2inside=Dict("pressure"=>"PRS","temperature"=>"TEM","humidity"=>"RHU","winSpeed"=>"WIN_S_AVG_10mi","winDir"=>"WIN_D_AVG_10mi")
    varNameS_inside=map(x->dic_pub2inside[x],varNameS)

    # SELECT station_code as staID,data_time as obsTime,prs ,tem ,rhu as humidity,win_s_inst as winSpeed,win_d_inst as winDir FROM tc_prod_zhfx.weather_data_original where station_code= 50136;
    obsData=fill(NaN,length(Times),length(varNameS),length(staIDs))
    
    if length(staIDs) < 1000
        # read the specified stations
        sqlquery=join(["SELECT station_code as staID,data_time as obsTime,",join(varNameS_inside,","),
            " FROM ",tableName,
            " WHERE duration_code = 6 and (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )",
            " AND  (station_code in ('",join(staIDs,"','"),"') ) "
            ]);
    else
        # read all stations, becase ready by station will be too slow if there are too many stations
        sqlquery=join(["SELECT station_code as staID,data_time as obsTime,",join(varNameS_inside,","),
            " FROM ",tableName,
            " WHERE duration_code = 6 and (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )"
            ]);
    end
    if disp;println("SQL Query: ", sqlquery);end
    response=DBInterface.execute(conn, sqlquery)
    close(conn)
    data_all=DataFrame(response)
    data_all.obsTime=ZonedDateTime.(data_all.obsTime,tz"Asia/Shanghai")

    # extract the wanted data and put into obsData
    # for n_t in eachindex(Times)
    Threads.@threads for n_t in eachindex(Times)
        # if disp;println("processing station: ",n_sta," ",staIDs[n_sta]);end
        ind_all_currT=findall(data_all.obsTime.==Times[n_t])
        if length(ind_all_currT) >0
            q_validInd,indOut=ismember(data_all.staID[ind_all_currT],staIDs)
            data2write=Float64.(replace!(Matrix(data_all[ind_all_currT[q_validInd],varNameS_inside]),missing=>NaN)')
            # data2write[data2write .< 0 ] .= NaN
            obsData[n_t,:,indOut[q_validInd]]=data2write
        end
    end
    return obsData
end



"""

Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,5,2)
varNameS=["PM25","PM10","SO2","NO2","CO","O3"]
staInfo=getValidStaInfo(Times)
@time getQCData(Times,varNameS,staInfo.staID)
"""
function getQCData(Times,varNameS,staIDs;tableName="tc_prod_zhfx.t_station_data_qc",disp=false,connInfo = getConnInfo())
    conn=getConn(connInfo)
    Times,varNameS,staIDs=in2vec.([Times,varNameS,staIDs])
    if eltype(Times) == DateTime
        Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
    end

    obsData=fill(NaN,length(Times),length(varNameS),length(staIDs))
    for n_var in eachindex(varNameS)
        # conn = getConn()
        if length(staIDs) < 2000
            # read the specified stations
            sqlquery=join(["SELECT data_time as obsTime,station_code as staID, value ",
                " FROM ",tableName,
                " WHERE (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )",
                " AND  (station_code in ('",join(staIDs,"','"),"') ) and element_type = '",varNameS[n_var],"' and isValid = true"
                ]);
        else
            # read all stations, becase ready by station will be too slow if there are too many stations
            sqlquery=join(["SELECT data_time as obsTime,station_code as staID, value ",
                " FROM ",tableName,
                " WHERE (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )",
                " and element_type = '",varNameS[n_var],"' and isValid = true"
                ]);
        end
        if disp;println("SQL Query: ", sqlquery);end
        response=DBInterface.execute(conn, sqlquery)
        # close(conn)
        data_aVar=DataFrame(response)
        data_aVar.obsTime=ZonedDateTime.(data_aVar.obsTime,tz"Asia/Shanghai")

        # extract the wanted data and put into obsData
        # Threads.@threads for n_t in eachindex(Times)
        for n_t in eachindex(Times)
            # if disp;println("processing station: ",n_sta," ",staIDs[n_sta]);end
            ind_all_currT=findall(data_aVar.obsTime.==Times[n_t])
            if length(ind_all_currT) >0
                q_validInd,indOut=ismember(data_aVar.staID[ind_all_currT],staIDs)
                obsData[n_t,n_var,indOut[q_validInd]]=replace!(data_aVar.value[ind_all_currT[q_validInd]],missing=>NaN)
            end
        end
    end
    close(conn)
    return obsData
end

"""

Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,5,2)
varNameS=["PM25","PM10","SO2","NO2","CO","O3"]
staInfo=getValidStaInfo(Times)
@time getQCData(Times,varNameS,staInfo.staID)
"""
function getQCAnalysis(Times,varNameS,staIDs;tableName="tc_prod_zhfx.t_station_data_qc",disp=false,connInfo = getConnInfo())
    conn=getConn(connInfo)
    Times,varNameS,staIDs=in2vec.([Times,varNameS,staIDs])
    if eltype(Times) == DateTime
        Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
    end

    
    # read the specified stations
    sqlquery=join(["SELECT * ",
        " FROM ",tableName,
        " WHERE (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )",
        " AND  (station_code in ('",join(staIDs,"','"),"') ) and (element_type in ('",join(varNameS,"','"),"') )",
        " order by data_time"
        ]);

    if disp;println("SQL Query: ", sqlquery);end
    response=DBInterface.execute(conn, sqlquery)
    close(conn)
    df_analysis=DataFrame(response)
    rename!(df_analysis,Dict(:station_code => :staID,:data_time => :obsTime,:element_type => :pollName))
    df_analysis.obsTime=ZonedDateTime.(df_analysis.obsTime,tz"Asia/Shanghai")
    for key in ["isValid","isOutRange","isPM25OutweightPM10"]
        df_analysis[!,key]=convert(Vector{Union{Missing,Bool}},df_analysis[!,key])
    end
    return df_analysis
end

# function getObsPoll(Times,varNameS,staIDs;tableName="tc_prod_zhfx.station_qx_data_original",disp=false)
#     Times,varNameS,staIDs=in2vec.([Times,varNameS,staIDs])
#     if eltype(Times) == DateTime
#         Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
#     end

#     obsData=fill(NaN,length(Times),length(varNameS),length(staIDs))
#     conn = getConn()
#     for n_sta in eachindex(staIDs)
#         if disp;println("processing station: ",n_sta," ",staIDs[n_sta]);end
#         sqlquery=join(["SELECT data_time as obsTime,",join(varNameS,","),
#             " FROM ",tableName,
#             " WHERE (station_code = '",staIDs[n_sta],"' ) ",
#             " AND (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' )"
#         ]);
#         response=DBInterface.execute(conn, sqlquery)
#         data_aSta=DataFrame(response)
#         data_aSta.obsTime=ZonedDateTime.(data_aSta.obsTime,tz"Asia/Shanghai")
#         ind=indexin(data_aSta.obsTime,Times)
#         q_validInd=map(!isnothing,ind)
#         obsData[ind[q_validInd],:,n_sta]=replace!(Matrix(data_aSta[q_validInd,varNameS]),missing=>NaN)
#     end
#     close(conn)
#     return obsData
# end

function query2df(sqlquery;connInfo = getConnInfo())
    # 
    conn=getConn(connInfo)
    response=DBInterface.execute(conn, sqlquery)
    close(conn)
    return DataFrame(response)
end


"""
using Dates
Times=DateTime(2023, 10, 1):Hour(1):DateTime(2023,11,1)
staID=getValidStaInfo(Times)
"""
function getValidStaInfo(Times;tableName="tc_prod_zhfx.station_qx_data_original",varName=nothing,connInfo = getConnInfo())

    if isnothing(varName)
        sqlquery=join(["SELECT station_code as staID,count(*) as count",
        " FROM ",tableName,
        " WHERE duration_code = 6 and data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'",
        " group by staID"  ]);
    elseif lowercase.(varName) == "all"
        sqlquery=join(["SELECT station_code as staID,count(*) as count",
            " FROM ",tableName,
            " WHERE  data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'",
            " And ( pm25 >= 0  or pm10 >=0 or so2 >=0 or no2 >=0 or co >=0 or o3>=0)",
            " group by staID"]);
    else
        sqlquery=join(["SELECT station_code as staID,count(*) as count",
        " FROM ",tableName,
        " WHERE duration_code = 6 and (data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"') ",
        " And ",varName," >= 0 ",
        " group by staID"    ]);
    end
    response=query2df(sqlquery;connInfo = connInfo)
    staInfo_registed=getStaInfo(;connInfo=connInfo)
    lia,lob=ismember(staInfo_registed.staID,response.staID)
    staInfo_valid=staInfo_registed[lia,:]
    staInfo_valid.count=response.count[lob[lia]]
    return staInfo_valid
end
getValidStaInfo(Times::AbstractVector{DateTime},varName::AbstractString;tableName="tc_prod_zhfx.station_qx_data_original") = getValidStaInfo(Times;tableName=tableName,varName=varName)

"""
using Dates
Times=DateTime(2023, 10, 1):Hour(1):DateTime(2023,11,1)
staID=getValidStaID(Times)
"""
function getValidStaID(Times;tableName="tc_prod_zhfx.station_qx_data_original",varName=nothing)
    if isnothing(varName)
        sqlquery=join(["SELECT distinct station_code as staID",
            " FROM ",tableName,
            " WHERE duration_code = 6 and  data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'",
             ]);
    elseif lowercase.(varName) == "all"
        sqlquery=join(["SELECT distinct station_code as staID",
            " FROM ",tableName,
            " WHERE duration_code = 6 and  data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'",
            " And ( pm25 >= 0  or pm10 >=0 or so2 >=0 or no2 >=0 or co >=0 or o3>=0)"
             ]);
        # println(sqlquery)
    else
        sqlquery=join(["SELECT distinct station_code as staID",
            " FROM ",tableName,
            " WHERE  duration_code = 6 and data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'",
            " And ",varName," >= 0 "
            ]);
    end
    staID=query2df(sqlquery)
    return staID
end
getValidStaID(Times::AbstractVector{DateTime},varName::AbstractString;tableName="tc_prod_zhfx.station_qx_data_original") =  getValidStaID(Times;tableName=tableName,varName=varName)
# function getValidStaID(Times,varName;tableName="tc_prod_zhfx.station_qx_data_original")
#     sqlquery=join(["SELECT distinct station_code as staID",
#         " FROM ",tableName,
#         " WHERE  data_time between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'",
#         " And ",varName," >= 0 "
#          ]);
#     staID=query2df(sqlquery)
#     return staID
# end


getStaInfo(staID::String) = getStaInfo([staID])
"""
staInfo=getStaInfo()
@time getStaInfo(staInfo.staID[1:2000])
staInfo=getStaInfo(staInfo.staID)
"""
function getStaInfo(staIDs::Union{Nothing,Array{String}}=nothing;tableName="tc_prod_zhfx.site_info",connInfo = getConnInfo(),types=nothing)
    conn=getConn(connInfo)
    if isnothing(staIDs)
        # sqlquery="SELECT base_site.site_code, base_site.site_name, base_site.lng, base_site.lat, base_site.site_type, dic_province.province_name, dic_city.city_name, dic_county.county_name
        # FROM tc_prod_zhfx.base_site 
        # left JOIN tc_prod_zhfx.dic_city ON base_site.city_code = dic_city.city_code
        # left JOIN tc_prod_zhfx.dic_province ON dic_city.province_code = dic_province.province_code
        # left JOIN tc_prod_zhfx.dic_county ON base_site.county_code = dic_county.county_code;";
        sqlquery="SELECT site_code,site_name,site_type,lng,lat,county_name,city_name,province_name FROM $tableName;";
        response=DBInterface.execute(conn, sqlquery)
        staInfo=DataFrame(response)
        # change field names
        rename!(staInfo,Dict(:site_code => :staID,:site_name => :staName,:site_type => :type,:lng => :lon,:province_name => :province,:city_name => :city,:county_name=>:county))
        # convert staID
        staInfo[!,:staID]=convert(Vector{String},staInfo[!,:staID])
        # convert lon lat to Float64
        for field in [:lon,:lat]
            replace!(staInfo[!,field],missing => NaN)
            staInfo[!,field]=convert(Vector{Float64},staInfo[!,field])
        end
    else
        if length(staIDs) < 1000
            sqlquery="SELECT site_code,site_name,site_type,lng,lat,county_name,city_name,province_name FROM $tableName
            where site_code in ('$(join(staIDs,"' , '"))');";
            response=DBInterface.execute(conn, sqlquery)
            staInfo=DataFrame(response)
            # change field names
            rename!(staInfo,Dict(:site_code => :staID,:site_name => :staName,:site_type => :type,:lng => :lon,:province_name => :province,:city_name => :city,:county_name=>:county))
            # convert staID
            staInfo[!,:staID]=convert(Vector{String},staInfo[!,:staID])
            # convert lon lat to Float64
            for field in [:lon,:lat]
                replace!(staInfo[!,field],missing => NaN)
                staInfo[!,field]=convert(Vector{Float64},staInfo[!,field])
            end
        else
            staInfo_all=getStaInfo(;connInfo=connInfo)
            staInfo=staInfo_all[in.(staInfo_all.staID,[staIDs]),:]
        end
    end
    for field in [:staName,:type,:province,:city,:county]
        replace!(staInfo[!,field],missing => "") # replace missing with ""
        staInfo[!,field]=string.(staInfo[!,field]) # set types to string
    end
    if !isnothing(types)# extract types            
        staInfo=staInfo[in.(staInfo.type, [types]),:];
    end
    close(conn)
    return staInfo
end
# # testing code 
# columnnames={"obsTime","staID","lon","lat","PM25","PM10","CO","NO2","SO2","O3"};
# varNameS=["PM25","PM10","CO","NO2","SO2","O3"]
# Times=DateTime(2015, 8, 24):Hour(1):DateTime(2015, 8, 29, 23)
# obsData=getObsPoll(Times,varNameS,["110000041"])


function insertDF(df;tableName="tc_prod_zhfx.t_station_data_qc",connInfo = getConnInfo())    
    conn=getConn(connInfo)
    data2write=(;pairs(eachcol(df))...)
    keynames=string.(keys(data2write))
    sqlquery= "INSERT INTO " * tableName * " ( " * join(keynames,", ") * " ) VALUES ( " * join(["?" for n in eachindex(keynames)]," , ") * ")"
    
    println("SQL Query: ", sqlquery)
    stmt = DBInterface.prepare(conn,sqlquery)
    DBInterface.executemany(stmt,data2write)
    
    # for n_row in eachrow(df)
    #     data2write=(;pairs(n_row)...)
    #     sqlquery= "INSERT INTO " * tableName * " ( " * join(keynames,", ") * " ) VALUES ( " * join(["?" for n in eachindex(keynames)]," , ") * ")"
    #     DBInterface.execute(conn, sqlquery)
    # end
    
    close(conn)
end

"""
Times=vcat(DateTime(2014) , DateTime(2014,6) , DateTime(2015) , collect(DateTime(2015,3):Month(3):DateTime(2020)) , collect(DateTime(2020,2):Month(2):DateTime(2022)) , collect(DateTime(2022,1):Month(1):DateTime(2030)) )
printPartion(Times)
"""
function printPartion(Times)
    for (n,t) in enumerate(Times)
        println("PARTITION p", n -1 ," VALUES LESS THAN (TO_DAYS('" ,Dates.format(t,"Y-mm-dd") ,"')),")
    end
    println("PARTITION p",length(Times)," VALUES LESS THAN MAXVALUE")
end
# function in2vec(in)
#     if  in isa StepRange;       in=collect(in);        end
#     if !(in isa AbstractVector);      in=[in];          end
#     return in
# end
end # module ZZDBM

