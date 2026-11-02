
module MeshDataM
export meshDATA,setCustomDimOrder
using NCDatasets,TimeZones,Dates,Base

# a class for read/write mesh data from/to nc file
# get and put parameters [... nZ, leadTime, variables, tStarts]
# there are three set of dimensions:
# default dimensions:[lon(1), lat(2), nZ(3), leadTime(4), variables(5), tStarts(6)]
# dimensions in nc file:[lon*lat, nZ, leadTime, variables, tStarts],(lon and lat are arranged as one dimension and set to a disk block to speed up read and write);
# custom dimensions: set the dimension of the input and output of this class by setOutsideDimOrder;
# Copyright: Huangjian Wu   All rights reserved.
# If you have any problem, please contact me: wuhuangjian@126.com

# ********************* update ************
# 2016-09-18     original version, 
# example usage: outData=meshObj.get(nZs,leadTimes,varNames,tStarts)
#                   meshObj.put(dataAll,nZs,leadTimes,varNames,tStarts)

struct meshStruct
    fileName::String
    sizeLonLat::Tuple
    nZs::Vector
    leadTimes::Vector
    varNames::Vector
    tStarts::Vector
    customDimOrder::Vector
end
# constrcut functions ***************************************
# construct function from file---------
function meshDATA(fileName::AbstractString)
    # read parameters from nc file
    ds=NCDataset(fileName,"r")
    nZs=ds["nZs"][:]
    sizeLonLat=(Int.(ds["sizeLonLat"][:])...,)
    leadTimes=ds["leadTimes"][:]
    varNames=split(String(vec(ds["varNames"][:])))
    # get tStart ************************************************************************************************
    tStarts_stringVec =ds["tStarts"][:,:] |> eachrow |> collect .|> String 
    if length(tStarts_stringVec[1])>19
        tStarts=ZonedDateTime.(tStarts_stringVec, "yyyy-mm-ddTHH:MM:SSzzzz")
    else
        tZone=ds["TimeZone"][:] |> vec |> String |> TimeZone
        tStarts = ZonedDateTime.(DateTime.(tStarts_stringVec,"yyyy-mm-dd HH:MM:SS"),tZone)
    end
    close(ds)
    return meshStruct(fileName,sizeLonLat,nZs,leadTimes,varNames,tStarts,[])
end

# create meshData file and write the data ---------
function meshDATA(fileName::AbstractString,sizeLonLat,nZs,leadTimes,varNames,tStarts,dataAll)
    meshObj= meshDATA(fileName,sizeLonLat,nZs,leadTimes,varNames,tStarts)
    Base.write(meshObj,dataAll)
    return meshObj
end
# create meshData file --------------------
function meshDATA(fileName::AbstractString,sizeLonLat,nZs,leadTimes,varNames,tStarts;chunksize=[prod(sizeLonLat),1,1,1,1])
    nZs,leadTimes,varNames,tStarts=in2vec.([nZs,leadTimes,varNames,tStarts])
    if isfile(fileName)
        @warn fileName * " already exist, would not change the file and will get parameters from it"
        return meshDATA(fileName)
    end
    if sizeLonLat isa Vector;sizeLonLat=(sizeLonLat...,);end
    if !(nZs isa Vector);nZs=[nZs];end
    if !(leadTimes isa Vector);leadTimes=[leadTimes];end
    if !(varNames isa Vector);varNames=[varNames];end
    if !(tStarts isa Vector);tStarts=[tStarts];end
    ds=NCDataset(fileName,"c",format= :netcdf4)
    # for nZs 
    defVar(ds,"nZs",convert.(Float64,nZs),("DnZs",))
    # for sizeLonLat 
    defVar(ds,"sizeLonLat",convert.(Float64,[sizeLonLat...,]),("info_sizeLonLat",))
    # for leadTimes 
    defVar(ds,"leadTimes",convert.(Float64,leadTimes),("DleadTimes",))
    # for varNames 
    varNames_str=join(varNames," ")
    defVar(ds,"varNames",reshape(collect(varNames_str),1,:),("one","info_varNames"))
    # for time tStarts
    tStart_vec=collect.(Dates.format.(tStarts, "yyyy-mm-ddTHH:MM:SSzzzz"))
    defVar(ds,"tStarts",permutedims(reduce(hcat,tStart_vec),[2,1]),("DtStarts","info_time"))
    # for dataAll  [lon lat nZ leadTime variables tStarts]
    ds.dim["Dlonlat"]=Int(prod(sizeLonLat))
    ds.dim["DvarNames"]=Int(length(varNames))
    if isnothing(chunksize)
        defVar(ds,"dataAll",Float32,("Dlonlat","DnZs","DleadTimes","DvarNames","DtStarts"),fillvalue=NaN)
    else
        defVar(ds,"dataAll",Float32,("Dlonlat","DnZs","DleadTimes","DvarNames","DtStarts"),fillvalue=NaN,chunksizes=chunksize)
    end
    close(ds)
    return meshStruct(fileName,sizeLonLat,nZs,leadTimes,varNames,tStarts,[])
end

function Base.read(obj::meshStruct;nZs=obj.nZs,leadTimes=obj.leadTimes,varNames=obj.varNames,tStarts=obj.tStarts)
    nZs,leadTimes,varNames,tStarts=in2vec.([nZs,leadTimes,varNames,tStarts])
    if eltype(tStarts) != eltype(obj.tStarts); @error "tStarts must in the format of $(eltype(obj.tStarts))";end
    # get the start and end of in the dataAll
    inds_z=getInds(nZs,obj.nZs)
    inds_l=getInds(Float32.(leadTimes),Float32.(obj.leadTimes))
    inds_v=getInds(uppercase.(varNames),uppercase.(obj.varNames))
    inds_t=getInds(tStarts,obj.tStarts)
    # print("z",zs,ze,"l",ls,le,"v",vs,ve,"t",ts,te)

    # read the data from nc file
    outData=Array{Union{Missing, Float32}}(undef,(prod(obj.sizeLonLat),length(inds_z),length(inds_l),length(inds_v),length(inds_t)))

    NCDataset(obj.fileName) do ds
        # outData=ds["dataAll"][:,inds_z,inds_l,inds_v,inds_t]
        for nz in eachindex(inds_z)
            iz_curr,nz_curr=getin(inds_z,nz)
            for nl in eachindex(inds_l)
                il_curr,nl_curr=getin(inds_l,nl)
                for nv in eachindex(inds_v)
                    iv_curr,nv_curr=getin(inds_v,nv)
                    for nt in eachindex(inds_t)
                        it_curr,nt_curr=getin(inds_t,nt)

                        outData[:,nz_curr,nl_curr,nv_curr,nt_curr]=ds["dataAll"][:,iz_curr,il_curr,iv_curr,it_curr]

                        if inds_t isa UnitRange;   break;  end
                    end
                    if inds_v isa UnitRange;   break;  end
                end
                if inds_l isa UnitRange;   break;  end
            end
            if inds_z isa UnitRange;   break;  end
        end
    end
    outData=nomissing(outData,NaN)
    # reshape the data to ouput
    sizeRead=size(outData);
    sizeOut=(obj.sizeLonLat...,sizeRead[2:end]...)
    outData=reshape(outData,sizeOut);
    if !isempty(obj.customDimOrder)
        outData=permutedims(outData,obj.customDimOrder);
    end
    return outData
end
function Base.write(obj::meshStruct,data2write;nZs=obj.nZs,leadTimes=obj.leadTimes,varNames=obj.varNames,tStarts=obj.tStarts)

    # input check
    nZs,leadTimes,varNames,tStarts=in2vec.([nZs,leadTimes,varNames,tStarts])
    if  eltype(tStarts) != ZonedDateTime; @error "Times must in the format of ZonedDateTime";end

    # prepare the data as the default dimemsion
    if !isempty(obj.customDimOrder)
        ind_custom=sortperm(obj.customDimOrder)
        data2write=permutedims(data2write,ind_custom);
    end
    # prepare the data as the size in nc, combine lon and lat
    sizeInput=size(data2write)
    sizeInNC=(prod(sizeInput[1:2]),sizeInput[3:end]...)
    data2write=reshape(data2write,sizeInNC);
    
    # get the start and end of the data2write
    inds_z=getInds(nZs,obj.nZs)
    inds_l=getInds(Float32.(leadTimes),Float32.(obj.leadTimes))
    inds_v=getInds(uppercase.(varNames),uppercase.(obj.varNames))
    inds_t=getInds(tStarts,obj.tStarts)

    # read the data from nc file
    data2write=Float32.(data2write)
    NCDataset(obj.fileName,"a") do ds
        # ds["dataAll"][:,inds_z,inds_l,inds_v,inds_t]=Float32.(data2write)
        for nz in eachindex(inds_z)
            iz_curr,nz_curr=getin(inds_z,nz)
            for nl in eachindex(inds_l)
                il_curr,nl_curr=getin(inds_l,nl)
                for nv in eachindex(inds_v)
                    iv_curr,nv_curr=getin(inds_v,nv)
                    for nt in eachindex(inds_t)
                        it_curr,nt_curr=getin(inds_t,nt)
            
                        ds["dataAll"][:,iz_curr,il_curr,iv_curr,it_curr]=data2write[:,nz_curr,nl_curr,nv_curr,nt_curr]
        
                        if inds_t isa UnitRange;   break;  end
                    end
                    if inds_v isa UnitRange;   break;  end
                end
                if inds_l isa UnitRange;   break;  end
            end
            if inds_z isa UnitRange;   break;  end
        end
    end
end

function setCustomDimOrder(obj::meshStruct,order)
    # the first two dimension must be lon and lat, (1,2)
    obj.customDimOrder=order;
end

function getin(inds,n)
    if inds isa UnitRange;
        iDisk_curr ,nMemory_curr= inds, :; 
    else
        iDisk_curr ,nMemory_curr= inds[n], n; 
    end
    return iDisk_curr ,nMemory_curr
end
function getInds(sample,all)
    inds=indexin(sample,all)
    if length(inds)==1
        ind_out=inds[1]
    else
        if length(inds) == inds[end] - inds[1] +1
            # range
            ind_out=inds[1]:inds[end]
        else
            # list
            ind_out=inds
        end
    end
    return ind_out
end
function in2vec(in)
    if  in isa StepRange;       in=collect(in);        end
    if !(in isa AbstractVector);      in=[in];          end
    return in
end

end#module