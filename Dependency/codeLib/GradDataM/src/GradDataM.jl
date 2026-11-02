module GradDataM
export gradData,GradDataStr,read_aVar,write_aVar,changeDatFileName,writeCtl,setTime!,get_nZ,setDatFileName!
using Base,Dates
# sample usage:
# read multiple variables in a line:
#  lonM,latM=read([gradData("lonlatM_d2.ctl")],["lonM","latM"],1)
"""
updates:
2023-04-30 add check for the range of nZ for read(grd,varName,nZ)
2023-10-15 add function to write a layer Base.write(obj::GradDataStr,data2write::Array,varName::String,nZ::Int)
"""
mutable struct GradDataStr
    datFileName::String
    datFormatString::String
    lonlatsize::Tuple
    # # one variables with one z
    varName_all::Vector # cell for variable name, upper case
    nZ_all::Vector# vector for z of variables
    varStart_all::Vector
end
function changeDatFileName(obj::GradDataStr,fileName::AbstractString)
	obj.datFileName=fileName
end
function setDatFileName!(obj::GradDataStr,fileName::AbstractString)
	obj.datFileName=fileName
    return obj
end
function setTime!(obj::GradDataStr,t)
    datFileName=obj.datFormatString
    d2replace=Dict("%y4"=>"yyyy","%y2"=>"yy","%m2"=>"mm","%d2"=>"dd","%h2"=>"HH","%n2"=>"MM","%s2"=>"SS","%d3"=>"ddd","%d4"=>"dddd")
    for (k,v) in d2replace
        datFileName=replace(datFileName,k => Dates.format(t,v))
    end
    obj.datFileName=datFileName
    return obj
end
# function Base.length(grd::GradDataStr)
#     return 1
# end
# function Base.iterate(grd::GradDataStr)
#     return grd
# end
function gradData(fileName;disp=true)
    local grd
    open(fileName) do f
        # read the ctl and split
        lines_all=readlines(f)
        splitedLines_all=split.(lines_all)
        firstWord_all=uppercase.([x[1] for x in splitedLines_all])
        # get lonlatsize -------
        getInd(pattern)=findfirst(x->x==pattern,firstWord_all)
        if in("PDEF",firstWord_all)
            lonlatsize=parse.(Int,splitedLines_all[getInd("PDEF")][2:3])
        else
            if disp; @warn "PEDF not found in $fileName, get the lonlatsize from XDEF and YDEF"; end
            lon=parse(Int,splitedLines_all[getInd("XDEF")][2])
            lat=parse(Int,splitedLines_all[getInd("YDEF")][2])
            lonlatsize=[lon,lat]
        end
        lonlatsize=(lonlatsize...,)
        # get datFileName ---------------
        datFormatString=joinpath(dirname(fileName) , splitedLines_all[getInd("DSET")][2][2:end])
        # datFileName=joinpath(dirname(fileName),splitedLines_all[getInd("DSET")][2][2:end])
        # get varName_all, varStart_all and nZ_all------
        ind_varStart=getInd("VARS")+1
        ind_varEnd=getInd("ENDVARS")-1
        varName_all=firstWord_all[ind_varStart:ind_varEnd]
        nZ_all=[parse(Int,x[2]) for x in splitedLines_all[ind_varStart:ind_varEnd]]
        @. nZ_all[nZ_all==0]=1
        if 1==length(nZ_all)
            varStart_all=[0]
        else
            varStart_all=[0 ;prod(lonlatsize)*cumsum(nZ_all[1:end-1])]
        end

        grd=GradDataStr(datFormatString,datFormatString,lonlatsize,varName_all,nZ_all,varStart_all)
    end
    return grd
end
function Base.copy(grd::GradDataStr)
    return GradDataStr(grd.datFileName,grd.datFormatString,grd.lonlatsize,grd.varName_all,grd.nZ_all,grd.varStart_all)
end
function gradData(ctlName,nxny,varNameS,zS;disp=true,datFileName=splitext(basename(ctlName))[1]*".grd",varDiscription=fill("",length(varNameS)))
    # write the ctl file if not exist
    if !isfile(ctlName)
        writeCtl(ctlName,nxny,varNameS,zS;datFileName=datFileName,varDiscription=varDiscription)
    end
    return gradData(ctlName;disp=disp)
end
function writeCtl(ctlName,nxny,varNameS,zS;datFileName=splitext(basename(ctlName))[1]*".grd",varDiscription=fill("",length(varNameS)))
    # make the path if not exist
    pathStr=dirname(ctlName);  if !isdir(pathStr); mkpath(pathStr); end
    # write the file
    open(ctlName,"w") do f
        println(f,"DSET ^"*datFileName);
        println(f,"TITLE not specified ");
        println(f,"OPTIONS byteswapped  template");
        println(f,"UNDEF     -1.E20");
        println(f,"PDEF $(nxny[1])  $(nxny[2])  LCC 1. 1.0");
        println(f,"XDEF $(nxny[1]) LINEAR 1. 1.0");
        println(f,"YDEF $(nxny[2]) LINEAR 1. 1.0");
        println(f,"ZDEF   $(max(zS...)) LEVELS");
        println(f,"VARS  $(length(varNameS))");
        for n_var in eachindex(varNameS)
            println(f,"$(varNameS[n_var])    $(zS[n_var])  0  $(varDiscription[n_var])");
        end
        println(f,"ENDVARS");
     end
end

function get_nZ(obj::GradDataStr,varName::String)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    return obj.nZ_all[ind]
end
function read_aVar(obj::GradDataStr,varName::String)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    dataOut=Array{Float32}(undef,(obj.lonlatsize...,obj.nZ_all[ind]))
    open(obj.datFileName) do f
        seek(f,obj.varStart_all[ind]*4)
        read!(f,dataOut)
    end
    dataOut=ntoh.(dataOut)
    return dataOut
end
function Base.read(obj::GradDataStr,varName::String)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    if isnothing(ind)
        @error varName * " not found !!"
        return
    end
    dataOut=Array{Float32}(undef,(obj.lonlatsize...,obj.nZ_all[ind]))
    open(obj.datFileName) do f
        seek(f,obj.varStart_all[ind]*4)
        read!(f,dataOut)
    end
    dataOut=ntoh.(dataOut)
    return dataOut
end
function Base.read(obj::GradDataStr,varName::String,nZ::Int)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    if isnothing(ind)
        @error varName * " not found !!"
        return
    end
    if nZ > obj.nZ_all[ind] || nZ < 1
        @error "nZ=$nZ is not in the nZ range of $varName, which is 1 to $(obj.nZ_all[ind])"
        return
    end
    dataOut=Array{Float32}(undef,obj.lonlatsize)
    open(obj.datFileName) do f
        seek(f,(obj.varStart_all[ind]+prod(obj.lonlatsize)*(nZ-1))*4)
        read!(f,dataOut)
    end
    dataOut=ntoh.(dataOut)
    return dataOut
end
function Base.read(obj::GradDataStr,varName::String,nZ_all::Vector)
    if length(nZ_all)==1
        dataOut=Base.read(obj,varName,nZ_all[1])
    else
        dataOut=Array{Float32,3}(undef,(obj.lonlatsize...,length(nZ_all)))
        for n_z=1:length(nZ_all)
            dataOut[:,:,n_z]=Base.read(obj,varName,nZ_all[n_z])
        end
    end
    return dataOut
end
function Base.read(obj::GradDataStr,varName::String,nZ_range::UnitRange)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    if isnothing(ind)
        @error varName * " not found !!"
        return
    end
    if step(nZ_range) == 1
        dataOut=Array{Float32,3}(undef,(obj.lonlatsize...,length(nZ_range)))
        open(obj.datFileName) do f
            seek(f,(obj.varStart_all[ind]+prod(obj.lonlatsize)*(first(nZ_range)-1))*4)
            read!(f,dataOut)
        end
        dataOut .= ntoh.(dataOut)
    else
        @error "step(nZ_range) != 1 is not supported"
        return
    end
    return dataOut
end
function Base.read!(data2read,obj::GradDataStr,varName::String,nZ_range::UnitRange)
    # this version may be slower than the Base.read version
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    if isnothing(ind)
        @error varName * " not found !!"
        return
    end
    if step(nZ_range) == 1
        # dataOut=Array{Float32,3}(undef,(obj.lonlatsize...,length(nZ_range)))
        open(obj.datFileName) do f
            seek(f,(obj.varStart_all[ind]+prod(obj.lonlatsize)*(first(nZ_range)-1))*4)
            read!(f,data2read)
        end
        data2read .= ntoh.(data2read)
    else
        @error "step(nZ_range) != 1 is not supported"
        return
    end
    return data2read
end
function createFile(obj::GradDataStr)    
    open(obj.datFileName,"w") do f
        seek(f,(obj.varStart_all[end]+prod(obj.lonlatsize)*obj.nZ_all[end] -1)*4)
        Base.write(f,hton.(Float32.(NaN)))
    end
end
function Base.write(obj::GradDataStr,data2write::Array,varName::String)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    data2write=hton.(Float32.(data2write))
    if ! isfile(obj.datFileName);createFile(obj);end
    open(obj.datFileName,"r+") do f
        seek(f,obj.varStart_all[ind]*4)
        Base.write(f,data2write)
    end
end
function Base.write(obj::GradDataStr,data2write::Array,varName::String,nZ::Int)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    data2write=hton.(Float32.(data2write))
    if ! isfile(obj.datFileName);createFile(obj);end
    if isnothing(ind)
        @error varName * " not found !!"
        return
    end
    if nZ > obj.nZ_all[ind] || nZ < 1
        @error "nZ=$nZ is not in the nZ range of $varName, which is 1 to $(obj.nZ_all[ind])"
        return
    end
    open(obj.datFileName,"r+") do f
        seek(f,(obj.varStart_all[ind]+prod(obj.lonlatsize)*(nZ-1))*4)
        Base.write(f,data2write)
    end
end
function write_aVar(obj::GradDataStr,data2write::Array,varName::String)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    data2write=hton.(Float32.(data2write))
    if ! isfile(obj.datFileName);createFile(obj);end
    open(obj.datFileName,"r+") do f
        seek(f,obj.varStart_all[ind]*4)
        Base.write(f,data2write)
    end
end


# ===========================zhaozheng=============================
function Base.write(obj::GradDataStr,data2write::Array,varName::String,nZ_range::UnitRange)
    ind=findfirst(x->x==uppercase(varName),obj.varName_all)
    data2write=hton.(Float32.(data2write))
    if ! isfile(obj.datFileName);createFile(obj);end
    if isnothing(ind)
        @error varName * " not found !!"
        return
    end
    if step(nZ_range) == 1
        # println("writing ",nZ_range," of ",varName," in ",obj.datFileName,"  ",obj.varStart_all[ind])
        open(obj.datFileName,"r+") do f
            seek(f,(obj.varStart_all[ind]+prod(obj.lonlatsize)*(first(nZ_range)-1))*4)
            Base.write(f,data2write)
        end
    else
        @error "step(nZ_range) != 1 is not supported"
        return
    end
end

end