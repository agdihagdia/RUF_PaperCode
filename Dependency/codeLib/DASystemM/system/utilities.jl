using Revise,DASystemM,CSV,utilities,NaNStatistics
using TimeZones,NaqpmsUtilities
# 2022-12-09 To ensure the integraty of the afterDA file (datFile_afterDA), the OI data are write the the tmp (datFile_afterDA_tmp)
# 2022-12-13 Add support for designating workDir to call the OI from other folder
# 2022-12-13 Add support for supressing prints by disp=false
function extract2BeforeDA(Times; replace=true,disp=true,  workDir=pwd(),addSand=false,  grd_ori=gradData("$workDir/test.oriNaqpms.ctl";disp=false),
        datFun_ori=nothing,polls2extract=["PM25","PM10","SO2","NO2","CO","O3"])
    if Times isa ZonedDateTime;Times=[Times];end
    # ---------------- set datFun_ori
    function datFun_ori_default(t)
        t=astimezone(t,tz"UTC")
        str_startDay=Dates.format(t-Hour(12),"yyyymmdd")
        return @sprintf("%s/naqpmsFolder/%s/out/data/%s12/testd1.%s.grd",workDir,str_startDay,str_startDay,Dates.format(t,"yyyymmddHH"))
    end
    if isnothing(datFun_ori);        datFun_ori=datFun_ori_default; end
    # ---------------- get grd_out
    grd_out=gradData("$workDir/beforeDA/pollutants.ctl",grd_ori.lonlatsize,polls2extract,fill(1,size(polls2extract));disp=false, datFileName="testd1.%y4%m2%d2%h2.grd");
    # ----------------- extract
    Threads.@threads for t_curr in Times
        grd_ori_curr=copy(grd_ori)
        grd_out_curr=copy(grd_out)
        grd_ori_curr.datFileName=datFun_ori(t_curr)
        if !isfile(grd_ori_curr.datFileName); println("Input file not ready, skipped: " * grd_ori_curr.datFileName);continue;end
        setTime!(grd_out_curr,astimezone(t_curr,tz"Asia/Shanghai"))
        if !replace && isfile(grd_out_curr.datFileName); 
            if disp; println("Output file exist, skipped: " * grd_out_curr.datFileName);  end
            continue;
        end
        if disp;    println("Extracting ...: " * grd_ori_curr.datFileName); end
        for poll in polls2extract
            data2write=toMass!(combinePollFromGrd(grd_ori_curr,poll;addSand=addSand),poll)
            write(grd_out_curr,data2write,poll)
        end
    end
end

function getStaObsFromTxt(t;disp=true, workDir=pwd(),
    file_staInfo="$workDir/staInfo_Nation.csv",
    file_obs=@sprintf("%s/obs/%s0000.txt",workDir,astimezone(t,tz"Asia/Shanghai"))
    )
    # read staInfo
    staInfo=CSV.read(file_staInfo,DataFrame,types=[String,Float64,Float64])
    if !isfile(file_obs);
        if disp;println("Observation file not ready: ",file_obs);end
        return ;
    end
    # get observations
    staObs=txt2staData(file_obs,staInfo)
    return staObs
end

# function generateEnState(beforeDA::AbstractArray{T};meanRatio=0.2,eleRatio=0.2,member=9)
#     return generateSpread(beforeDA;meanRatio=meanRatio,eleRatio=eleRatio,member=member) .+ beforeDA
# end
function generateSpread(beforeDA::AbstractArray{T};meanRatio=0.3,eleRatio=0.1,member=9) where T
    meanRatio=convert(T,meanRatio)
    eleRatio=convert(T,eleRatio)
    spread= ( meanRatio*mean(beforeDA) .+ eleRatio .* beforeDA ) .* 
                reshape(convert.(T,range(-1,1;length=member)),1,1,member)
    return spread
end