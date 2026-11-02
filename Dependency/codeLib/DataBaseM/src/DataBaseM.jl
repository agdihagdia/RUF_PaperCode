module DataBaseM
using MySQL,DataFrames,Dates,DBInterface,TimeZones
export getObsPoll,getStaInfo,getValidStaID,getValidStaInfo,query2df,insertDF
include("carefulImportPoll.jl")
# include("plotDataBase.jl")

"""

updates:
2023-04-24: use ZonedDateTime in getObsPoll

"""


"""



Sample Usage:
Times=DateTime(2022,5):Hour(1):DateTime(2022,6)
varNameS=["PM25","PM10","SO2","NO2","CO","O3"]
staInfo=getStaInfo()
getObsPoll(Times,varNameS,staInfo.staID)
"""
function getObsPoll(Times,varNameS,staIDs;tableName="verified.pollutant")
    Times,varNameS,staIDs=in2vec.([Times,varNameS,staIDs])
    if eltype(Times) == DateTime
        Times=ZonedDateTime.(Times,tz"Asia/Shanghai")
    end

    obsData=fill(NaN,length(Times),length(varNameS),length(staIDs))
    conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"])
    for n_sta in eachindex(staIDs)
        sqlquery=join(["SELECT obsTime,",join(varNameS,","),
            " FROM ",tableName,
            " WHERE ((staID = '",staIDs[n_sta],"' ) ",
            " AND (obsTime between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"' ))"
        ]);
        response=DBInterface.execute(conn, sqlquery)
        data_aSta=DataFrame(response)
        data_aSta.obsTime=ZonedDateTime.(data_aSta.obsTime,tz"Asia/Shanghai")
        ind=indexin(data_aSta.obsTime,Times)
        q_validInd=map(!isnothing,ind)
        obsData[ind[q_validInd],:,n_sta]=replace!(Matrix(data_aSta[q_validInd,varNameS]),missing=>NaN)
    end
    close(conn)
    return obsData
end

function query2df(sqlquery)
    conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"])
    response=DBInterface.execute(conn, sqlquery)
    close(conn)
    return DataFrame(response)
end

function getValidStaInfo(Times;validRatio=0,tableName="original.pollutant")
    staInfo_registed=getStaInfo()
    staInfo=getValidStaID(Times;validRatio=validRatio,tableName=tableName)
    
    # add empty fields to staInfo_curr
    missingVec=Vector{Union{String,Missing}}(missing,length(staInfo.staID))
    staInfo.staName = copy(missingVec)
    staInfo.province = copy(missingVec)
    staInfo.city = copy(missingVec)
    staInfo.type = copy(missingVec)

    # update staInfo for national stations
    ind=indexin(staInfo.staID,staInfo_registed.staID)
    q2update=.!isnothing.(ind)
    fields2update=["lon","lat","staName","city","province","type"]
    staInfo[q2update,fields2update] .= staInfo_registed[ind[q2update],fields2update]

    # exclude station without lon or lat
    q_valid=  .!ismissing.(staInfo.lon) .& .!ismissing.(staInfo.lat)
    if sum(q_valid) != size(staInfo,1)
        @warn "removing $(size(staInfo,1) - sum(q_valid)) staInfo without lon or lat"
        staInfo=staInfo[q_valid,:]
    end
    return staInfo

end
# function getValidStaInfo(Times;validRatio=0,tableName="verified.pollutant_julia")
#     sqlquery=join(["SELECT $(tableName).staID,stainfo.staName,stainfo.province,stainfo.city,stainfo.lon,stainfo.lat,count(*)  ",
#         " FROM ",tableName, " join verified.stainfo on stainfo.staID = $(tableName).staID ",
#         " WHERE obsTime between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'  group by staID"
#         ]);
#         # println(sqlquery)
#     staInfo=query2df(sqlquery)
#     rename!(staInfo,:"count(*)" => :count)
#     return staInfo
# end
# getValidStaID = getValidStaInfo
function getValidStaID(Times;validRatio=0,tableName="original.pollutant",getLonLat=true)
    if getLonLat # slow version, with lonlat
        sqlquery=join(["SELECT distinct staID,lon,lat,count(*) ",
            " FROM ",tableName,
            " WHERE obsTime between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'  group by staID"
            ]);
        staInfo=query2df(sqlquery)
        rename!(staInfo,:"count(*)" => :count)
    else # quick version, withou lonlat (6 times faster)
        sqlquery=join(["SELECT distinct staID,count(*)",
            " FROM ",tableName,
            " WHERE obsTime between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'  group by staID"
        ]);
        staInfo=query2df(sqlquery)
        missingVec=Vector{Union{Float64,Missing}}(missing,length(staInfo.staID))
        staInfo.lon = copy(missingVec)
        staInfo.lat = copy(missingVec)
    end
    return staInfo
end

function getValidStaInfo(Times,varName;tableName="original.pollutant")
    staInfo_registed=getStaInfo()
    staInfo_valid=getValidStaID(Times,varName;tableName=tableName)
    return staInfo_registed[in.(staInfo_registed.staID,[staInfo_valid.staID]),:]
end

function getValidStaID(Times,varName;tableName="original.pollutant")
        sqlquery=join(["SELECT staID,lon,lat,count(*) ",
            " FROM ",tableName,
            " WHERE $varName is not null and ",
            " obsTime between '",Dates.format(Times[1],"Y-m-d H:M:S"),"' and '",Dates.format(Times[end],"Y-m-d H:M:S"),"'  group by staID"
            ]);
        staInfo=query2df(sqlquery)
        rename!(staInfo,:"count(*)" => :count)
    return staInfo
end
getStaInfo(staID::String) = getStaInfo([staID])

function getStaInfo(staIDs::Union{Nothing,Array{String}}=nothing;tableName="original.stainfo_all")
    conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"])
    if isnothing(staIDs)
        sqlquery="SELECT * FROM $(tableName)";
        response=DBInterface.execute(conn, sqlquery)
        staInfo=DataFrame(response)
    else
        # for future Implementation
        sqlquery=join(["SELECT * FROM $(tableName) WHERE staID in ('"   join(staIDs,"','")    "')"]);
        response=DBInterface.execute(conn, sqlquery)
        staInfo=DataFrame(response)
    end
    for field in ["lon","lat"]
        # replace!(staInfo[:,field],missing => NaN)
        staInfo[!,field]=convert(Vector{Float64},staInfo[!,field])
    end
    close(conn)
    return staInfo
end
# # testing code 
# columnnames={"obsTime","staID","lon","lat","PM25","PM10","CO","NO2","SO2","O3"};
# varNameS=["PM25","PM10","CO","NO2","SO2","O3"]
# Times=DateTime(2015, 8, 24):Hour(1):DateTime(2015, 8, 29, 23)
# obsData=getObsPoll(Times,varNameS,["110000041"])


function insertDF(df;tableName="verified.pollutant")    
    conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],multi_statements=true)

    data2write=(;pairs(eachcol(df))...)
    keynames=string.(keys(data2write))
    sqlquery= "INSERT INTO " * tableName * " ( " * join(keynames,", ") * " ) VALUES ( " * join(["?" for n in eachindex(keynames)]," , ") * ")"
    stmt = DBInterface.prepare(conn,sqlquery)
    DBInterface.executemany(stmt,data2write)
    
    close(conn)
end

function printPartion(Times)
    for (n,t) in enumerate(Times)
        println("PARTITION p", n -1 ," VALUES LESS THAN (TO_DAYS('" ,Dates.format(t,"Y-mm-dd") ,"')),")
    end
    println("PARTITION p",length(Times)," VALUES LESS THAN MAXVALUE")
end
"""
Times=vcat(DateTime(2014) , DateTime(2014,6) , DateTime(2015) , collect(DateTime(2015,3):Month(3):DateTime(2020)) , collect(DateTime(2020,2):Month(2):DateTime(2022)) , collect(DateTime(2022,1):Month(1):DateTime(2030)) )
printPartion(Times)
"""

function in2vec(in)
    if  in isa StepRange;       in=collect(in);        end
    if !(in isa AbstractVector);      in=[in];          end
    return in
end
end # module DataBaseM
