
using GradDataM,StaDataM,NaqpmsUtilities,utilities
export extractStaFromGrd

    # % modification:
    # % 2022-06-07 original version
    # 2023-03-02 force Times to be Vector
    # 2023-04-27 add option addSand
    
# extractStaFromGrd, but get the datFileName from ctlFile
function extractStaFromGrd(ctlFile,Times,lonM,latM,staInfo;nZ=1,varNames=[],disp=true,varType="direct",transUnit=false,addSand=true,version="default")
    # % extract forecast sta data from grd file

    # % input check ***********************
    Times, varNames = in2vec.([Times, varNames])
    grd=gradData(ctlFile;disp=disp)

    # % initialize ***************************
    oMapper=obsMapper(lonM,latM,staInfo.lon,staInfo.lat);
    staDataAll=fill(NaN,(length(Times),length(varNames),length(staInfo.lon)))

    # % extract *********************    
    Threads.@threads for n_t in eachindex(Times)
    # for n_t=1:length(Times)
        grd_curr=setTime!(copy(grd),Times[n_t])
        datFile=grd_curr.datFileName
        if disp;println("Extracting " * datFile * "...");end
        if !isfile(datFile)
            if disp; @warn "File not exist: $datFile ...";end
            continue
        end
        for (n_var,var_curr) in enumerate(varNames)
            # read from grd --------------------
            if varType == "direct"
                data_curr=read(grd_curr,var_curr,nZ)
            elseif varType == "combine"
                data_curr=combinePollFromGrd(grd_curr,var_curr,nZ;addSand=addSand,version=version)
            else
                @error "varType not found: $varType"
                data_curr=Array{Float32,2}(NaN,grd.lonlatsize)
            end
            # unit ---------------
            if transUnit == true && (uppercase.(var_curr) in ["SO2","NO2","CO","O3"])
                toMass!(data_curr,var_curr)
            end
            staDataAll[n_t,n_var,:]=map2sta(oMapper,data_curr);
        end
    end
    staForeOut=staDATA(Times,varNames,staInfo,[0],staDataAll)
    return staForeOut
end

function extractStaFromGrd(ctlFile,datFileNameFun,Times,lonM,latM,staInfo;nZ=1,varNames=[],disp=true,varType="direct",transUnit=false,addSand=true)
    # % extract forecast sta data from grd file

    # % input check ***********************
    if  !( varNames isa Vector);       varNames=[varNames];        end
    if  Times isa StepRange;       Times=collect(Times);        end
    grd=gradData(ctlFile;disp=disp)

    # % initialize ***************************
    oMapper=obsMapper(lonM,latM,staInfo.lon,staInfo.lat);
    staDataAll=fill(NaN,(length(Times),length(varNames),length(staInfo.lon)))

    # % extract *********************    
    Threads.@threads for n_t in eachindex(Times)
    # for n_t=1:length(Times)
        grd_curr=copy(grd)
        datFile=datFileNameFun(Times[n_t])
        grd_curr.datFileName=datFile
        if disp;println("Extracting " * datFile * "...");end
        if !isfile(datFile)
            if disp; @warn "File not exist: $datFile ...";end
            continue
        end
        for n_var in eachindex(varNames)
            var_curr=varNames[n_var]
            # read from grd --------------------
            if varType == "direct"
                data_curr=read(grd_curr,var_curr,nZ)
            elseif varType == "combine"
                data_curr=combinePollFromGrd(grd_curr,var_curr,nZ;addSand=addSand)
            else
                @error "varType not found: $varType"
                data_curr=Array{Float32,2}(NaN,grd.lonlatsize)
            end
            # unit ---------------
            if transUnit == true && (uppercase.(var_curr) in ["SO2","NO2","CO","O3"])
                toMass!(data_curr,var_curr)
            end
            staDataAll[n_t,n_var,:]=map2sta(oMapper,data_curr);
        end
    end
    staForeOut=staDATA(Times,varNames,staInfo,[0],staDataAll)
    return staForeOut
end

