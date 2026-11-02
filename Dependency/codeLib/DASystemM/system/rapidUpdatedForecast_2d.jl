using Revise,Printf,Dates,TimeZones,MAT,DelimitedFiles,DBInterface,MySQL,NCDatasets
using DASystemM,StaDataTool,LocalizerM,GradDataM,Statistics,StaDataM,NaqpmsUtilities,DataFrames,TOML#,Infiltrator

function log(x...)
    println(now(), " ", join(x, " ")...)
    flush(stdout)
end


# ===================================== Edit this for each implement ==============================
# ------------------- Prepare the CTM inputs -------------------------------------
function naqFileFun(tStart,t_curr;domain=2)
    tStartStr=Dates.format(astimezone(tStart,tz"UTC"),"yyyymmddHH")
    t_currStr=Dates.format(astimezone(t_curr,tz"UTC"),"yyyymmddHH")
    fileName="/sh2/home/3clear/bedrock1/runs/$tStartStr/china_3_27km_v05/naqp/apops/out_nc/NAQPMS.$t_currStr.out.grd0$domain.nc"
    return fileName
end
function checkGrdFinished(grd,Times;disp=false)
    for t in Times
        grd=setTime!(copy(grd),t)
        if !isfile(grd.datFileName)
            if disp;println("File not exist: ",grd.datFileName);end
            return false
        end
    end
    return true
end
function prepareBaseFore(tStart,Times;baseFileFun=naqFileFun,projectFolder=pwd(),overwrite=false,disp=false,
    varOut=["U","V","W","DZ","TEMP","PLEV","PM10"],
    varIn=["U","V","W","DZ","TEMP","PLEV",["PM10","aer0201","aer0202","aer0203","aer0204"]]
    )
    tStartStr=Dates.format(astimezone(tStart,tz"UTC"),"yyyymmddHH")
    oriFolder=dirname(baseFileFun(tStart,tStart))
    baseFolder="$projectFolder/DATA/baseFolder/$tStartStr"
    outFileFun(t)="$baseFolder/ctm.$(Dates.format(astimezone(t,tz"UTC"),"yyyymmddHH")).grd"
    baseCTL="$baseFolder/ctm.ctl"
    staFileName_Base="$baseFolder/staData.$tStartStr.mat"

    # skip if original folder not exist
    if !isdir(oriFolder) 
        @info "oriFolder not exist: $oriFolder"
        return 
    end
    # skip if already exist
    if !overwrite && isfile(baseCTL) && isfile(staFileName_Base) && checkGrdFinished(gradData(baseCTL),astimezone.(Times,tz"UTC");disp=disp)
        @info "baseFolder already finished: $baseFolder"
        return 
    end

    if !isdir(baseFolder);mkdir(baseFolder);end
    # link original folder
    file_link="$baseFolder/original"
    if islink(file_link); rm(file_link); end
    symlink(oriFolder,file_link; dir_target = true)
    grd=nothing
    surfaceFore_all=Vector{Any}(undef,length(Times))
    for (n_t,t_curr) in enumerate(Times)
        # ------------------- read data
        # if disp;println(" read data ...........");end
        baseFile=baseFileFun(tStart,t_curr)
        outFile=outFileFun(t_curr)
        # if isfile(outFile) && !overwrite
        #     if disp;println("File already exist: $outFile");end
        #     continue
        # end
        ds=NCDataset(baseFile )
        if disp;println(" extracting: $baseFile ...");end
        data_all=Vector{Any}(undef,length(varOut))
        for (n_var,var) in enumerate(varOut)
            varSub=varIn[n_var]
            if varSub isa AbstractString
                data_all[n_var]=ds[varSub][:,:,:]
            else
                data_all[n_var]=sum(collect(ds[varSubSub][:,:,:] for varSubSub in varSub))
            end
        end
        surfaceFore_all[n_t]=cat(collect(M[:,:,1] for M in data_all)...;dims=4)
        # ------------------ write ctl and get grd
        if isnothing(grd)            
            nxny=size(data_all[1])[1:2]
            zS=collect(size(x,3) for x in data_all)
            writeCtl(baseCTL,nxny,varOut,zS;datFileName=raw"ctm.%y4%m2%d2%h2.grd")
            grd=gradData(baseCTL)
        end
        # ------------------ write data 
        # if disp;println(" write data ...........");end
        grd.datFileName=outFile
        for (n_var,var) in enumerate(varOut)
            write(grd,data_all[n_var],var)
        end
    end
    # extract staData ------------------
    surfaceFore_all=cat(surfaceFore_all...;dims=3)    
    leads=Dates.value.(Dates.Hour.(Times .- Times[1]))
    extractStaData(staFileName_Base,surfaceFore_all,Times[1],leads,varOut)
end
# ------------------------------------- prepare the obs --------------------------------------
function getConn(tag="directYBS")
    if isequal(tag,"directYBS")
        conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],port=9030)
    elseif isequal(tag,"downYBS")
        conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],port=39030)
    end
    return conn
end
# ------------ prepare staObs ----------------
function getStaObs_aT(t;pollName="PM10")
    conn = DBInterface.connect(MySQL.Connection,ENV["RUF_DB_HOST"],ENV["RUF_DB_USER"],ENV["RUF_DB_PASSWORD"],port=9030)
    tStr=Dates.format(astimezone(t,tz"Asia/Shanghai"),"yyyy-mm-dd HH:MM:SS")
    sqlquery="SELECT lon,lat,$(pollName)_1h as $(pollName),dtable.code as staID,datadate FROM  v00183clear.tw_relation_station_info AS rtable  LEFT JOIN v00183clear.t_sc_pollutant_s_hour as dtable on  rtable.CODE = dtable.CODE where rtable.air_met=1 and datadate='$tStr'"
    response=DBInterface.execute(conn, sqlquery)
    data_aT=DataFrame(response)
    q_invalid=ismissing.(data_aT[!,pollName]) .| (data_aT[!,pollName] .< 0)
    data_aT=data_aT[.!q_invalid,:]
    staInfo=DataFrame(lon=convert.(Float64,data_aT.lon),lat=convert.(Float64,data_aT.lat),staID=convert.(String,data_aT.staID))
    staDataAll=Float64.(reshape(data_aT[:,pollName],1,1,:,1))
    staObs=staDATA([t],[pollName],staInfo,staDataAll)
    return staObs
end

"""
prepare RUF inputs -----------
make directory:
    DATA/RUF/  DATA/RUF/yyyymmddHH 

"""
function prepareRUFinputs(Times;overwrite=false,disp=false)
    baseFileFun(tStart,t)="DATA/baseFolder/$(Dates.format(astimezone(tStart,tz"UTC"),"yyyymmddHH"))/ctm.$(Dates.format(astimezone(t,tz"UTC"),"yyyymmddHH")).grd"
    tStart_Base=getTStartBase(Times;baseFileFun=baseFileFun)
    baseFolderCTM=dirname(baseFileFun(tStart_Base,tStart_Base))
    rufFolder="DATA/RUF/$(Dates.format(astimezone(Times[1],tz"Asia/Shanghai"),"yyyymmddHH"))"
    rufInputFolder="$rufFolder/inputs"
    if !isdir(rufFolder);mkdir(rufFolder);end
    if !isdir(rufInputFolder);mkdir(rufInputFolder);end
    inputFileFun(t)="$rufInputFolder/input.$(Dates.format(astimezone(t,tz"Asia/Shanghai"),"yyyymmddHH")).grd"
    if !isdir(rufInputFolder);mkdir(rufInputFolder);end
    # write inputs.ctl
    grdBase=gradData(baseFolderCTM.*"/ctm.ctl")    
    if disp;    log("Using baseline forecast: $baseFolderCTM");  end
    writeCtl(rufInputFolder*"/input.ctl",grdBase.lonlatsize,grdBase.varName_all,grdBase.nZ_all;datFileName="input.%y4%m2%d2%h2.grd")
    # link mesh input data
    for t in Times
        oriFile=baseFileFun(tStart_Base,t)# TimeZone of UTC
        outfile=inputFileFun(t) # TimeZone of Asia/Shanghai
        if isfile(outfile) && !overwrite
            # if disp;println("File already exist: $outfile");end
        else
            if islink(outfile) || isfile(outfile);rm(outfile);end
            symlink("../../../../"*oriFile,outfile) # link file
        end
    end
    # link staData
    tStartStr_Base=Dates.format(astimezone(tStart_Base,tz"UTC"),"yyyymmddHH")
    oriMatFile="DATA/baseFolder/$tStartStr_Base/staData.$tStartStr_Base.mat"
    outMatFile="$rufFolder/inputs/staData.mat"
    if islink(outMatFile) || isfile(outMatFile) ;rm(outMatFile);end
    symlink("../../../../"*oriMatFile,outMatFile) # link file

    # write output ctl    
    q_outVar= .![x in ["U","V","W","DZ","TEMP","PLEV"] for x in uppercase.(grdBase.varName_all)]
    writeCtl(rufFolder*"/ruf.ctl",grdBase.lonlatsize,grdBase.varName_all[q_outVar],grdBase.nZ_all[q_outVar];datFileName="ruf.%y4%m2%d2%h2.grd")
end
# ======================================= Common functions (Do not edit) =======================================
function getTStartBase(Times_RUF;
    baseFileFun=(tStart,t)->"DATA/baseFolder/$(Dates.format(astimezone(tStart,tz"UTC"),"yyyymmddHH"))/ctm.$(Dates.format(astimezone(t,tz"UTC"),"yyyymmddHH")).grd")
    tStart_ruf=astimezone(Times_RUF[1],tz"UTC")
    tStart_base_candidates = floor(tStart_ruf,Day) + Hour(12) .+Day.([0,-1,-2])
    if Hour(tStart_ruf) < Hour(12);   tStart_base_candidates = tStart_base_candidates - Day(1); end 
    # find the tStart_base with available base forecast
    tStartBase=ZonedDateTime(1990,tz"UTC")
    for t in tStart_base_candidates
        fileBase_final=baseFileFun(t,Times_RUF[end])
        # println(fileBase_final)
        if isfile(fileBase_final)
            println("Base files ready: ",fileBase_final)
            tStartBase=t
            break
        end
    end
    if tStartBase == ZonedDateTime(1990,tz"UTC")
        return nothing
    else
        return tStartBase
    end
end
function adjust3DbySurface(IC3d,aIC_surface)
    conMin=mean(aIC_surface)/1000
    aIC_surface[aIC_surface .< conMin] .= conMin
    IC_surface=IC3d[:,:,1];IC_surface[IC_surface .< conMin] .= conMin
    assimilatedRatio= aIC_surface ./ IC_surface
    replace!(assimilatedRatio, NaN => 1.0)
    aIC3d=IC3d[:,:,:] .* assimilatedRatio
    # limit the maximum of aIC
    surfaceMax=maximum(aIC3d[:,:,1,1])
    aIC3d[aIC3d .> surfaceMax] .= surfaceMax
    return aIC3d
end
includet("$(pathof(DASystemM))/../../system/"*"EnOI_anHour.jl")
function AsimilatingIC(t;disp=false,replace=false,pollNameS=["PM25","PM10","SO2","NO2","CO","O3"],
    staObs=getStaObs_aT(t))
    t_Beijing=astimezone(t,tz"Asia/Shanghai")
    if disp;    println("--------------- Generating Assimilated ICs ----------------"); end
    reanalysisFolder="DATA/EnOI"
    rufFolder="DATA/RUF/$(Dates.format(t_Beijing,"yyyymmddHH"))"
    grdIn=gradData(rufFolder*"/inputs/input.ctl")
    grdIn=setTime!(grdIn,t_Beijing)
    # extract beforeDA ---------------------
    grdBefore=gradData("$reanalysisFolder/beforeDA/beforeDA.ctl",grdIn.lonlatsize,pollNameS,fill(1,size(pollNameS));datFileName="beforeDA.%y4%m2%d2%h2.grd");# create the ctl file
    grdBefore=setTime!(grdBefore,t_Beijing)
    for pollName in pollNameS
        dataBefore=read(grdIn,pollName,1)
        write(grdBefore,dataBefore,pollName)
    end
    # write obs ------------
    write2mat(staObs,"$reanalysisFolder/obs/staData.$(Dates.format(t_Beijing,"yyyymmddHH")).mat")
    # do the OI --------------------------
    EnOI_anHour(t;polls2OI=pollNameS,replace=replace ,disp=disp  , workDir=reanalysisFolder,  staObs=staObs)    
    # construct and write 3d assimilated ICs ----------------
    grdAfter=setTime!(gradData("$reanalysisFolder/afterDA/afterDA.ctl"),t_Beijing)
    grdOut=setTime!(gradData(rufFolder*"/ruf.ctl"),t_Beijing)
    for pollName in pollNameS
        aIC_surface=read(grdAfter,pollName,1)
        IC3d=read(grdIn,pollName)
        aIC3d=adjust3DbySurface(IC3d,aIC_surface)
        if get_nZ(grdOut,pollName) == size(aIC3d,3)
            write(grdOut,aIC3d,pollName)
        else
            error("The number of layers of IC3d and grdOut should be the same")
        end
    end
    # write the info file ----------------
    file_info="$reanalysisFolder/afterDA/info."*Dates.format(t_Beijing,"yyyymmddHH")
    if !isfile(file_info) || replace
        info=Dict("pollNameS"=>pollNameS,
            "valid_obs"=>  vec(sum(.!isnan.(read(staObs;Times=t,varNames=pollNameS));dims=(1,3))))
        open(file_info, "w") do io
            TOML.print(io, info)
        end
    end
end

function forecast_RUF2d(Times_all,pollNameS;replace=false,disp=true,nonnagative=true,isDecay=true,par=true,
    dx=15e3,dy=15e3,
    rufFolder="DATA/RUF/$(Dates.format(astimezone(Times_all[1],tz"Asia/Shanghai"),"yyyymmddHH"))",
    )
    if disp
        log(" ================= Rapid Updated $pollNameS Forecast starting from $(Dates.format(astimezone(Times_all[1],tz"Asia/Shanghai"),"yyyymmddTHH")) =============")
    end
    # set folders and file names ------------
    
    grdIn=gradData(rufFolder*"/inputs/input.ctl")
    grdOut=gradData(rufFolder*"/ruf.ctl") # ruf.ctl is prepared by AsimilatingIC, and the asslimilated IC is already in the first file

    staFileName_RUF="$rufFolder/staData.$(Dates.format(astimezone(Times_all[1],tz"Asia/Shanghai"),"yyyymmddHH")).mat"
    # ==================  exit if not ready or finished =====================
    # assimilated IC not exist -------------
    grdOut=setTime!(grdOut,astimezone(Times_all[1],tz"Asia/Shanghai"));
    if !isfile(grdOut.datFileName)
        if disp;        log("assimilated IC not exist, passing... ",grdOut.datFileName); end
        return
    end
    # finished -------------
    grdOut=setTime!(grdOut,astimezone(Times_all[end],tz"UTC"))
    file_end=grdOut.datFileName
    if  !replace && isfile(file_end) && isfile(staFileName_RUF)
        if disp;        log("Forecast already finished, passing... ",file_end); end
        return
    end
    # ================== configuration =====================
    if !(pollNameS isa Vector);pollNameS=[pollNameS];end
    pollNameS=uppercase.(pollNameS)
    # set parallel -------------
    if par
        # set advection paralell if there are many threads
        if Threads.nthreads() > 1.5* length(pollNameS)
            par_var,par_adv=false,true
        else
            par_var,par_adv=true,false
        end
    else
        par_var,par_adv=false,false
    end
    # ================= preparation =====================
    # read all input data --------------------
    u_all=Array{Float32,3}(undef,grdIn.lonlatsize...,length(Times_all))
    v_all=Array{Float32,3}(undef,grdIn.lonlatsize...,length(Times_all))
    conBase_all=Array{Float32,4}(undef,grdIn.lonlatsize...,length(Times_all),length(pollNameS))
    Threads.@threads for n_t in eachindex(Times_all)
        t=Times_all[n_t]
        grdCurr=setTime!(copy(grdIn),astimezone(t,tz"Asia/Shanghai"))
        u_all[:,:,n_t]=read(grdCurr,"U",1)
        v_all[:,:,n_t]=read(grdCurr,"V",1)
        for (n_var,pollName) in enumerate(pollNameS)
            conBase_all[:,:,n_t,n_var]=read(grdCurr,pollName,1)
        end
    end
    # prepare assimilated ICs and  output conRUF_all ---------
    conRUF_all=fill(NaN32,size(conBase_all));
    grdOut=setTime!(grdOut,astimezone(Times_all[1],tz"Asia/Shanghai"))    
    for (n_var,pollName) in enumerate(pollNameS)
        aIC=read(grdOut,pollName,1)
        if all(0 .<= aIC .< Inf)
            conRUF_all[:,:,1,n_var]=aIC
        else
            error("The assimilated ICs contain NaN or Inf values")
        end
    end
    
    # ================== Initial Error Forecast =====================
    if disp;    log("do the Initial Error Forecast ...");   end
    if par_var 
        Threads.@threads for n_var in eachindex(pollNameS)
        # for n_var in eachindex(pollNameS)
            pollName=pollNameS[n_var]
            if disp;    log("Starting the RUF for "*pollName*"......"); end
            if isDecay
                decayRatio=getDecayRatio(pollName)
            else
                decayRatio=1.0
            end
            aIC=conRUF_all[:,:,1,n_var]
            conRUF_all[:,:,:,n_var]=RUF_aVar(conBase_all[:,:,:,n_var],aIC,u_all,v_all,decayRatio;par=par_adv,disp=false,nonnagative=nonnagative,dx=dx,dy=dy)
        end
    else
        for n_var in eachindex(pollNameS)
            pollName=pollNameS[n_var]
            if disp;    log("Starting the RUF for "*pollName*"......"); end
            if isDecay
                decayRatio=getDecayRatio(pollName)
            else
                decayRatio=1.0
            end
            aIC=conRUF_all[:,:,1,n_var]
            conRUF_all[:,:,:,n_var]=RUF_aVar(conBase_all[:,:,:,n_var],aIC,u_all,v_all,decayRatio;par=par_adv,disp=false,nonnagative=nonnagative,dx=dx,dy=dy)
        end
    end
    # ================== Write the RUF Forecast =====================
    Threads.@threads for n_t in eachindex(Times_all)
    # for n_t in eachindex(Times_all)
        t_curr=Times_all[n_t]
        grd_curr=setTime!(copy(grdOut),astimezone(t_curr,tz"Asia/Shanghai"))
        for (n_var,pollName) in enumerate(pollNameS)
            write(grd_curr,conRUF_all[:,:,n_t,n_var],pollName)
        end
    end
    # extract StaData
    leads_RUF=Dates.value.(Dates.Hour.(Times_all .- Times_all[1]))
    extractStaData(staFileName_RUF,conRUF_all,Times_all[1],leads_RUF,pollNameS)
end
function getDecayRatio(pollName)
    # accroding to "Interpreting aerosol lifetimes using the GEOS-Chem model and constraints from radionuclide measurements"
    if uppercase(pollName) in ["PM25" ,"TPM25"]     ;    lifeTime=24*6 ;  
    elseif uppercase(pollName) in ["PM10" ,"TPM10"] ;    lifeTime=24*3 ;
    elseif uppercase(pollName) in ["SAND"] ;    lifeTime=24*7 ;
    elseif uppercase(pollName) in ["O3"]            ;    lifeTime=8 ;   
    end
    decayRatio=1-1/lifeTime
    return decayRatio
end

function RUF_aVar(conBase_all,aIC,u_all,v_all,decayRatio;par=true,disp=true,nonnagative=true,dx=15e3,dy=15e3)
    # set par = true to enable multi threading in advDiffAnHour
    sizelonlat=size(conBase_all)[1:2]
    numT=size(conBase_all,3)
    
    # prepare dx dy -------------------
    DX=dx*ones(sizelonlat)
    DY=dy*ones(sizelonlat)
    # do the IETM
    conRUF_all=fill(NaN,(size(conBase_all)));
    ErrFromInit=Array{eltype(conRUF_all),2}(undef,sizelonlat...)
    for n_t = 1:numT
        if disp; log(n_t);  end
        if n_t==1
            ErrFromInit=conBase_all[:,:,1]-aIC
            conRUF_all[:,:,1]=aIC
        else
            # do the IETM ----------------------
            u=u_all[:,:,n_t]
            v=v_all[:,:,n_t]
            ErrFromInit=decayRatio*advDiffAnHour(ErrFromInit, u, v, DX, DY;par=par);
            conBase_curr=conBase_all[:,:,n_t]
            conRUF_curr=conBase_curr .- ErrFromInit
            if nonnagative
                conRUF_curr[conRUF_curr.<0] .= 0
            end
            conRUF_all[:,:,n_t]=conRUF_curr
            ErrFromInit=conBase_curr .- conRUF_curr;
        end
    end
    return conRUF_all
end
function extractStaData(staFileName,conFore_all,tStart_RUF,leads,pollNameS;projFolder=pwd())
    # the dimensions of conFore_all should be [lon,lat,lead,poll]
    if eltype(leads) == Hour; leads=map(x->x.value,leads); end
    # get lonM latM 
    lonM,latM=load("$projFolder/datagrid.jld2","lonM","latM")
    # read staInfo
    staInfo=CSV.read("staInfo.csv", DataFrame)   
    q_inDomain=getStaInDomain(lonM,latM,staInfo[:,:lon],staInfo[:,:lat])
    # q_inDomain=getStaInDomain(lonM,latM,staInfo.lon,staInfo.lat)
    staInfo=staInfo[q_inDomain,:]
    # extract 
    oMapper=obsMapper(lonM,latM,staInfo.lon,staInfo.lat)
    staDataAll=reshape(
        mapslices(x->map2sta(oMapper,x),conFore_all,dims=(1,2)),
        1,length(staInfo.lon),length(leads),length(pollNameS))
    staDataAll=permutedims(staDataAll,[1,4,2,3])
    staRUF=staDATA(tStart_RUF,pollNameS,staInfo,leads,staDataAll)
    write2mat(staRUF,staFileName)
end
