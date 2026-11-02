module utilities
include("boundMap.jl")
include("obsMapper.jl")
include("geoDist.jl")
include("pollUtilities.jl")
include("imageSlicing.jl")
include("peaks.jl")
include("upscaleMapper.jl")
include("unitTransformation.jl")
using CSV,DataFrames,NaNStatistics,Statistics
export ismember,getExtractExpr,in2vec,readStaInfo,getStatistics,cartesianWin2polar,meshgrid,daymean
function ismember(a,b)
    lob=indexin(a,b)
    lia=map(!isnothing,lob)
    return lia,lob
end
function cartesianWin2polar(u,v)
    # % u is west wind;v is south wind
    # % 0 is north win; pi/2 is east wind
    # % 2017-10-19 fix a bug, change +pi to 2*pi
    
    radian=angle.(-(v .+ u .* 1im));
    radian[radian .< 0] = radian[radian .<0 ] .+ 2*pi;
    speed=sqrt.(u.^2+v.^2);
    return radian,speed
end
function daymean(hourlyData::AbstractArray)
    # 获取每天平均值
    # first dim should be time
    # error if fisrt dim can not be divided by 24
    @assert size(hourlyData)[1]%24==0 "first dim should be hourly data, and can be divided by 24"
    sizeIn=size(hourlyData)
    sizeOut=(Int(sizeIn[1]/24),sizeIn[2:end]...)
    return reshape(nanmean(reshape(hourlyData,24,:),dims=1),sizeOut)
end
# extract variables in a dict to Main space
function getExtractExpr(d::Dict)
    expr = quote end
    for (k, v) in d
        push!(expr.args, :($(Symbol(k)) = $v))
    end
    return expr
end
#  针对矩阵的处理
function getStatistics(reference::AbstractMatrix{T}, predict::AbstractMatrix{T}) where T<:Number
    # 获取每一列的统计结果
    dataAll = [getStatistics(reference[:, i], predict[:, i]) for i in axes(reference, 2)]
    
    # 使用 keys() 获取字段名，并重组成一个含有 Vector 的 NamedTuple
    ks = keys(dataAll[1])
    return NamedTuple{ks}( ( [getfield(d, k) for d in dataAll] for k in ks ) )
end

#  针对向量的处理
function getStatistics(reference::AbstractVector,predict::AbstractVector)
    q_valid=.!isnan.(reference) .& .!isnan.(predict)
    pmean=mean(predict[q_valid])
    rmean=mean(reference[q_valid])
    pstd=std(predict[q_valid])
    rstd=std(reference[q_valid])
    cc=nancor(predict[q_valid],reference[q_valid])
    rmse=sqrt(mean((predict[q_valid]-reference[q_valid]).^2))
    crmsd=sqrt(mean((predict[q_valid] .- pmean .- (reference[q_valid].-rmean)).^2))
    nmb=mean(predict[q_valid]-reference[q_valid])./mean(reference[q_valid])#Mean normalized bias
    mae=mean(abs.(predict[q_valid]-reference[q_valid]))
    nmae=mae./mean(reference[q_valid])#Normalized mean absolute error
    return (rmse=rmse,mae=mae,corr=cc,bias=pmean-rmean,pstd=pstd,rstd=rstd,crmsd=crmsd,nmb=nmb,nmae=nmae)
end
function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

"""
sample usage:

Times,varNames,leadTimes,staIDs=in2vec.([Times,varNames,leadTimes,staIDs])
"""
function in2vec(in)
    if  in isa StepRange;       in=collect(in);        end
    if !(in isa AbstractVector);      in=[in];          end
    return in
end

function readStaInfo(fileName)
    staInfo=CSV.read(fileName,DataFrame)
    staInfo.staID=string.(staInfo.staID)
    for col in ["city","note","province","stationName","stationType"]
        if col in names(staInfo)
            staInfo[!,col]=string.(staInfo[!,col])
        end
    end
    # staInfo.note=string.(staInfo.note)
    return staInfo
end


export readAObsTxt,writeAObsTxt

function readAObsTxt(fileName::String;par=true)
    if par;ntask=Threads.nthreads();else;ntask=1;end
    dataInFile=CSV.read(fileName,DataFrame,ntasks=ntask,types=[String,Float32,Float32,Float32,Float32,Float32,Float32,Float32,Float32,Float32,Float32,Float32],header=["staID","lon","lat","pm25","pm10","co","no2","so2","o3","o38","ext1","ext2"])
    # replace nagavtive values with NaN for fields of  pm25     pm10     co       no2      so2      o3       o38  
    for col in ["pm25","pm10","co","no2","so2","o3","o38","ext1","ext2"]
        dataInFile[dataInFile[!,col].<0,col] .= NaN
    end
    return dataInFile
end

function writeAObsTxt(fileName,df)
    for col in ["o38","ext1","ext2"]
        df[!,col] .= -999
    end
    # conver NaN to -999
    for col in ["pm25","pm10","co","no2","so2","o3"]
        # df[ismissing.(df[!,col]),col] .= -999
        df[isnan.(df[!,col]),col] .= -999
    end
    CSV.write(fileName,df;delim='\t',writeheader=false)
end

end # module

