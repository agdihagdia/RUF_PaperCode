pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "codeLib")))
using Revise,GradDataM,Dates,TimeZones,Statistics,ZZDBM,DASystemM,PlotFunM
includet("$(pathof(DASystemM))/../../system/"*"rapidUpdatedForecast_2d.jl")

function runARUF(Times_RUF;pollNameS=["PM10"],overwrite=false,disp=true,animate=true,staObs=nothing,par=true)

    tStart_RUF=Times_RUF[1]
    # ======================== check wether can skip ===========================
    rufFolder="DATA/RUF/$(Dates.format(astimezone(tStart_RUF,tz"Asia/Shanghai"),"yyyymmddHH"))"
    if !isdir(rufFolder); mkdir(rufFolder); end
    # check the info file --------------------
    tStart_Base=getTStartBase(Times_RUF)
    if isnothing(staObs);        staObs=getStaObs_aT(tStart_RUF);   end
    file_info="$rufFolder/info.toml"
    info_new=Dict("tStart"=>Dates.format(astimezone(tStart_RUF,tz"Asia/Shanghai"),"yyyymmddHH"),
        "tStart_Base"=>Dates.format(astimezone(tStart_Base,tz"UTC"),"yyyymmddHH"),
        "pollNameS"=>pollNameS,
        "valid_obs"=>  vec(sum(.!isnan.(read(staObs;Times=tStart_RUF,varNames=pollNameS));dims=(1,3))),
        "isDone"=>false,"isDone_inputs"=>false,"isDone_aICs"=>false,"isDone_RUF"=>false,"isDone_animate"=>false)
    function writeInfo(fileName,info)
        open(fileName, "w") do io
            TOML.print(io, info)
        end
    end    
    if sum(info_new["valid_obs"]) < 1
        println("observation not ready, skip")
        return
    end
    if isfile(file_info)
        info_old=TOML.parsefile(file_info)
        if all(isequal(info_old[key],info_new[key]) for key in ["tStart","tStart_Base","pollNameS","valid_obs"])
            for key in ["isDone","isDone_inputs","isDone_aICs","isDone_RUF","isDone_animate"]
                info_new[key]=info_old[key] # inherite the  
            end
        else
            info_old=nothing # set to nothing if the forecast info is no the same
        end
    else 
        info_old=nothing 
    end
    if !overwrite && !isnothing(info_old)  && info_old["isDone"] 
        if disp;println("Forecast already exist with info: ",info_old);end
        return ;# nothing updated, skip
    else
        if disp;println("Starting forecast with info: ",info_new);end
    end
    # ================== prepare RUF Folder =================
    if overwrite || !info_new["isDone_inputs"]
        if disp;println("======================== preparing input folder...");end
        prepareRUFinputs(Times_RUF;overwrite=true,disp=disp)
        info_new["isDone_inputs"]=true;        writeInfo(file_info,info_new)
    end
    # ================== do the IC assimilation ================
    if overwrite ||  !info_new["isDone_aICs"]
        if disp;println("======================== Doing OI...");end
        AsimilatingIC(tStart_RUF;disp=disp,replace=true,pollNameS=pollNameS,staObs=staObs)
        info_new["isDone_aICs"]=true;        writeInfo(file_info,info_new)
    end
    # ================== do the RUF =======================
    if overwrite ||  !info_new["isDone_RUF"]
        if disp;println("======================== Doing RUF...");end
        forecast_RUF2d(Times_RUF,pollNameS;replace=true,disp=disp,nonnagative=true,isDecay=true,par=par,
            dx=9e3,dy=9e3 )
        info_new["isDone_RUF"]=true;        writeInfo(file_info,info_new)
    end
    # plot the forecast ---------------------
    if animate && ( overwrite ||  !info_new["isDone_RUF"])
        if disp;println("======================== generating animation...");end
        animateAFore(Times_RUF;rufFolder=rufFolder)
        info_new["isDone_animate"]=true;        writeInfo(file_info,info_new)
    end
    info_new["isDone"]=true;        writeInfo(file_info,info_new)
end
function animateAFore(Times_RUF;rufFolder="DATA/RUF/$(Dates.format(astimezone(Times_RUF[1],tz"Asia/Shanghai"),"yyyymmddHH"))")
    lonM,latM=load("datagrid.jld2","lonM","latM")
    grd=gradData("$rufFolder/ruf.ctl")
    data2plot=cat([read(setTime!(grd,astimezone(t,tz"Asia/Shanghai")),"PM10",1) for t in Times_RUF]...;dims=3)
    animateSandStorm(data2plot,Times_RUF,lonM,latM;dpi=100,
        outFileName="$rufFolder/dust storm started at $(Dates.format(tStart_RUF,"yyyy-mm-dd HH")).mp4")
end

cd("/home/zsm/IETM/RUF_2d")
# staObs=staDATA("../qc202504.mat")
staObs=staDATA("/home/zsm/数据库运维-Julia/data/original/ori202504.mat")
# ================= prepare base folder ======================
# tStart_baseAll=ZonedDateTime.(DateTime(2023,4,1,12):Day(1):DateTime(2023,4,15,12),tz"UTC")
# tStart_baseAll=ZonedDateTime.(DateTime(2025,4,11,12):Day(1):DateTime(2025,4,13,12),tz"UTC")#zsm2025/5/27
# for tStart in tStart_baseAll
#     # tStart = tStart_baseAll[2]
#     Times = tStart .+ Hour.(0:100)
#     varOut=["U","V","PM10"]
#     varIn=["U","V","PM10"]  
#     prepareBaseFore(tStart,Times;overwrite=false,disp=true,varOut=varOut,varIn=varIn)
# end

# ========================== do the RUF =========================
# tStart_RUF_all=ZonedDateTime.(DateTime(2023,4,2):Hour(1):DateTime(2023,4,16),tz"Asia/Shanghai") # ***********************
# tStart_RUF_all=ZonedDateTime.(DateTime(2025,4,19,1):Hour(1):DateTime(2025,4,21,0),tz"Asia/Shanghai") #zsm2025/5/27
tStart_RUF_all=ZonedDateTime(DateTime(2025,4,11,18),tz"Asia/Shanghai")
foreHour_RUF=Hour(3*24);
if !isempty(ARGS)
    iChunk=parse(Int,ARGS[1])
    numChunk=parse(Int,ARGS[2])
    tStart_RUF_all=tStart_RUF_all[iChunk:numChunk:end]
end
# Threads.@threads 
# @time for n_start in eachindex(tStart_RUF_all)
    # tStart_RUF=tStart_RUF_all[n_start]
    tStart_RUF=tStart_RUF_all
    Times_RUF=tStart_RUF .+ collect(Hour(0):Hour(1):foreHour_RUF)
    staObs_curr=extract(staObs;Times=[tStart_RUF])
    runARUF(Times_RUF;pollNameS=["PM10"],overwrite=true,disp=true,animate=false,staObs=staObs_curr,par=true)
# end

# using NCDatasets
# ds=NCDataset("/sh2/home/3clear/wuhj/paper_IETMsand/IETM_3d/DATA/baseFolder/2023040912/original/NAQPMS.2023041006.out.grd02.nc" )
