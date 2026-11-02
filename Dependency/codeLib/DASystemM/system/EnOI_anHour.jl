using Revise,DASystemM,StaDataTool,MyPlotM,MAT,CSV,utilities,NaNStatistics,JLD2
using TimeZones,NaqpmsUtilities
# 2022-12-09 To ensure the integraty of the afterDA file (datFile_afterDA), the OI data are write the the tmp (datFile_afterDA_tmp)
# 2022-12-13 Add support for designating workDir to call the OI from other folder
# 2022-12-13 Add support for supressing prints by disp=false

includet((@__DIR__) * "/utilities.jl")

"""
================ procedure to implement ===============
1. make folder workDir/beforeDA workDir/beforeDA 
2. make datagrid.jld2 to workDir/ containing lonM latM 
3. make data in beforeDA workDir/beforeDA/pollutants.ctl and corresponding afterDA.yyyymmddHH.grd
4. make staObs
5. run the EnOI_anHour()
"""
function EnOI_anHour(t;polls2OI=["PM25","PM10","SO2","NO2","CO","O3"],replace=true ,disp=true,par=false, workDir=pwd(),  staObs=nothing)

    # ============== check input and output files, skip if need ===============
    # check staObs ----------
    if isnothing(staObs)
        staObs=getStaObsFromTxt(t)
    end
    # check Inputs grd ----------
    grd_beforeDA=setTime!(gradData("$workDir/beforeDA/beforeDA.ctl";disp=false),astimezone(t,tz"Asia/Shanghai"))
    if !isfile(grd_beforeDA.datFileName) ; println("Inputs not ready: " * grd_beforeDA.datFileName );
        return;
    end
    # check Outputs grd ----------
    grd_afterDA=gradData("$workDir/afterDA/afterDA.ctl",grd_beforeDA.lonlatsize,polls2OI,fill(1,size(polls2OI));disp=false, datFileName="afterDA.%y4%m2%d2%h2.grd")
    grd_afterDA=setTime!(grd_afterDA,astimezone(t,tz"Asia/Shanghai"))
    datFile_afterDA=grd_afterDA.datFileName
    if isfile(datFile_afterDA)
        if replace == true;
            @warn "Replacing : " * datFile_afterDA;
        else
            if disp; println("Output exist, skipping: " * datFile_afterDA ); end
            return 
        end
    end

    # ============== prepare data for DA ===============
    # get lonM latM
    grd=gradData("$workDir/datagrid.ctl";disp=false)
    lonM=read(grd,"lonM",1)
    latM=read(grd,"latM",1)
    # (lonM,latM)=load("$workDir/datagrid.jld2","lonM","latM")
    # quanlity assurance for observations ----------
    # removing pm10 larger than 2500
    n_pm10=indexin(["pm10"],lowercase.(staObs.varNames))[1]    
    if !isnothing(n_pm10)
        pm10=staObs.staDataAll[:,n_pm10,:]
        println("Limit $(sum(pm10 .> 2500)) pm10 not larger than 2500")
        pm10[pm10 .> 2500] .= 2500
        staObs.staDataAll[:,n_pm10,:] .= pm10
    end

    # prepare grds of beforeDA and afterDA --------
    datFile_afterDA_tmp=@sprintf("%s/afterDA/afterDA.%s.tmp",workDir,Dates.format(astimezone(t,tz"Asia/Shanghai"),"yyyymmddHH"))
    grd_afterDA.datFileName=datFile_afterDA_tmp

    # ============== do the DA and write to tmp file ==============
    if disp;    print("Doing the OI for : ");   end
    local data_afterDA
    # Threads.@threads for poll in polls2OI
    for poll in polls2OI
        if disp;    print(poll,"... "); end
        # prepare enState ----------
        println("prepare enState at $(now())")
        data_beforeDA=read(grd_beforeDA,poll,1)
        T=eltype(data_beforeDA)
        enMesh=data_beforeDA .+ generateSpread(data_beforeDA)
 
        # super observing ------------
        println("super observing at $(now())")
        staIDs_valid=staObs.staInfo.staID[.!isnan.(vec(read(staObs;Times=t,varNames=poll)))]
        staObs_M=extract(staObs;Times=t,varNames=poll,staIDs=staIDs_valid) # extract the stations with valid data
        staObs_M=mergeStaByGrid(staObs_M,lonM,latM)
        staInfo_M=staObs_M.staInfo
        obs_M=convert.(T,vec(read(staObs_M;Times=t,varNames=poll)))
        println("megering $(length(staIDs_valid)) stations to $(length(obs_M)) grid points")
        # extract staSim
        println("extract staSim")
        mp=obsMapper(lonM,latM,staInfo_M.lon,staInfo_M.lat;method="toGrid")
        staSim=hcat([map2sta(mp,enMesh[:,:,n]) for n in axes(enMesh,3)]...);
        staSim=convert.(T,staSim)

        # prepare localizer ----------
        println("prepare localizer at $(now())")
        spaceCoef_cutoff=convert(T,150e3)
        stalonlat=convert.(T,hcat(staInfo_M.lon,staInfo_M.lat))
        localizer=localizerM(hcat(lonM[:],latM[:]),stalonlat;spaceCoef_cutoff=spaceCoef_cutoff);       
        space_cutoff_varing =  1 .* getVaringSpaceCutOff(localizer;validStaNum=1,par=par)
        space_cutoff_varing[space_cutoff_varing .> 150e3] .= 150e3;
        space_cutoff_varing[space_cutoff_varing .< 15e3] .=15e3;
        
        # do the DA ----------
        println("do the DA at $(now())")
        enState=reshape(enMesh,length(lonM),:)
        enState_updated=EnKF(enState,obs_M,staSim,localizer,space_cutoff_varing;EnKFtype="DEnKF",obsErrFun=obs-> fill(convert(T,1.0),size(obs)),par=par);
        data_afterDA=reshape(mean(enState_updated,dims=2),size(lonM))
        q_nan=isnan.(data_afterDA); data_afterDA[q_nan]=data_beforeDA[q_nan];# fill the nan with beforeDA
        data_afterDA[data_afterDA.<0].=0
        write(grd_afterDA,data_afterDA,poll)
    end

    # mv the tmp file to final file
    mv(datFile_afterDA_tmp,datFile_afterDA;force=true)

    if disp;    print("\n");    end
    return data_afterDA
end