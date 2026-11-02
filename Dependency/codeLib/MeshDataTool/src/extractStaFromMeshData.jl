
using  MeshDataM,StaDataTool

# extract forecast sta data from nc file with the format of meshData *****************************
# sample usage:
#     tMonth=datetime([2019 2 1],'TimeZone','Asia/Shanghai');
#     varNames={'pm25'};
#     outStaMatName='E:\dyml_biasCorrect\DSpostprocess\staPostprocess_pm25c_unbias';
#     fileNameFun=@(varName,nDayForecast) fullfile('E:','dyml_biasCorrect','DSpostprocess',[datestr(tMonth,'yyyymm') 'DSpost' num2str(nDayForecast) 'f_' varName '_unbias.nc']);
#     load staObs201902.mat staInfo
# parallel for multile file：
#     option.parType='byStart';%'byFile'
#     extractStaFromMeshData(fileName_all,...,option)
# modification:
    # 2021-10-12 original version
    # 2022-5-11 change the format of staInfo from dict to dataframe
    
function extractStaFromMeshData(fileName_all,lonM,latM,staInfo;nZs=[],leadTimes=[],varNames=[],tStarts=[],parType=false,disp=true)

    # % input check ***********************
    function in2vec(x)
        if x isa StepRange;            x=collect(x)
        elseif !(x isa Vector);       x=[x];        end
        return x
    end
    nZs,leadTimes,varNames,tStarts=in2vec.((nZs,leadTimes,varNames,tStarts))

    # % other inputs
    if fileName_all isa String;fileName_all=[fileName_all]; end
    meshObjIn=MeshDataM.meshDATA(fileName_all[1]);
    if isempty(nZs); nZs=meshObjIn.nZs; end
    if isempty(leadTimes); leadTimes=meshObjIn.leadTimes; end
    if isempty(varNames); varNames=meshObjIn.varNames; end
    if isempty(tStarts); tStarts=meshObjIn.tStarts; end
    if length(nZs)>1;@error "Only one nZs is supported, but the nZs are: " * string(nZs);end
    
    # % initialize ***************************
    oMapper=obsMapper(lonM,latM,staInfo.lon,staInfo.lat);
    staDataAll=fill(NaN,(length(tStarts),length(varNames),length(staInfo.lon),length(leadTimes)))

    # % extract *********************    
    if parType == "byFile"
        Threads.@threads for n_file in eachindex(fileName_all)
            meshObjIn=MeshDataM.meshDATA(fileName_all[n_file]);
            if disp;println("Extracting " * fileName_all[n_file] * "...");end
            lia_z,lia_l,lia_v,lia_t,nZs_aFile,lead_aFile,varNames_aFile,t_aFile=getLia_aFile(meshObjIn,nZs,leadTimes,varNames,tStarts);
            staData_aFile=extract_afile(meshObjIn,nZs_aFile,lead_aFile,varNames_aFile,t_aFile,oMapper;parType=parType,disp=disp);
            staDataAll[lia_t,lia_v,:,lia_l]=staData_aFile;
        end
    else
        for n_file in eachindex(fileName_all) # file itr---------
            meshObjIn=MeshDataM.meshDATA(fileName_all[n_file]);
            if disp;println("Extracting " * fileName_all[n_file] * "...");end
            lia_z,lia_l,lia_v,lia_t,nZs_aFile,lead_aFile,varNames_aFile,t_aFile=getLia_aFile(meshObjIn,nZs,leadTimes,varNames,tStarts);
            staData_aFile=extract_afile(meshObjIn,nZs_aFile,lead_aFile,varNames_aFile,t_aFile,oMapper;parType=parType,disp=disp)
            staDataAll[lia_t,lia_v,:,lia_l]=staData_aFile;
        end
    end
    staForeOut=staDATA(tStarts,varNames,staInfo,leadTimes,staDataAll)
    return staForeOut
end
function extract_afile(meshObj,nZs,leadTimes,varNames,tStarts,oMapper;parType=false,disp=true)
    sizeOut=(length(tStarts),length(varNames),length(oMapper.staLon),length(leadTimes))
    if 0==prod(sizeOut);   return [];    end
    staData_aFile=fill(NaN,sizeOut)
    if parType == "byLead" # most time is spend in reading nc, parallel do not help much 
        for n_start in eachindex(tStarts)
            tStart_curr=tStarts[n_start]
            dataReaded_aStart=read(meshObj;nZs=nZs,leadTimes=leadTimes,varNames=varNames,tStarts=tStart_curr)
            if disp;println();print("Extracting lead of: "); end
            Threads.@threads for n_lead in eachindex(leadTimes)
                if disp;print("$(leadTimes[n_lead]) "); end
                sta_aLeadStart=extractStaFromMesh(dataReaded_aStart[:,:,:,n_lead,:,:],oMapper)
                staData_aFile[n_start,:,:,n_lead]=permutedims(
                    reshape(sta_aLeadStart,length(oMapper.staLon),length(varNames))
                    ,(2,1))
            end
        end
    else
        for n_var in eachindex(varNames) # var itr---------
            varName_curr=varNames[n_var]
            dataReaded_aVar=read(meshObj;nZs=nZs,leadTimes=leadTimes,varNames=varName_curr,tStarts=tStarts)
            sta_aVar=fill(NaN,(length(tStarts),1,length(oMapper.staLon),length(leadTimes)))
            if disp;println("Extracting " * varName_curr * "..."); end
            if parType == "byStart"
                Threads.@threads for n_start=1:length(tStarts)#% start itr-----------
                    sta_aVar[n_start,1,:,:]=extractStaFromMesh(dataReaded_aVar[:,:,:,:,:,n_start],oMapper)
                end
            else
                for n_start=1:length(tStarts)#% start itr-----------
                    sta_aVar[n_start,1,:,:]=extractStaFromMesh(dataReaded_aVar[:,:,:,:,:,n_start],oMapper)
                end
            end
            staData_aFile[:,n_var,:,:]=sta_aVar;
        end
    end
    return staData_aFile
end

function extractStaFromMesh(dataReaded_avarNamestart,oMapper)
    sta_avarNamestart=mapslices(x->map2sta(oMapper,x),dataReaded_avarNamestart,dims=(1,2));
    return sta_avarNamestart
end
function getLia_aFile(meshObjIn,nZs,leadTimes,varNames,tStarts)
    lia_z=in.(nZs,[meshObjIn.nZs]);
    lia_l=in.(leadTimes,[meshObjIn.leadTimes]);
    lia_v=in.(varNames,[meshObjIn.varNames]);
    lia_t=in.(tStarts,[meshObjIn.tStarts]);
    nZs_aFile=nZs[lia_z];
    lead_aFile=leadTimes[lia_l];
    varNames_aFile=varNames[lia_v];
    t_aFile=tStarts[lia_t];
    return lia_z,lia_l,lia_v,lia_t,nZs_aFile,lead_aFile,varNames_aFile,t_aFile
end
