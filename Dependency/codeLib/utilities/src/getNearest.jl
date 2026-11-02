
using Distances
export getNearestData,getNearestIJ,getNearestLonLat,getNearest,getUnstructNearest

function getUnstructNearest(unLon,unLat,unPoll,lonV,latV)
    unLonLat_points=eachrow(hcat(unLon,unLat))
    lonLatV_points=eachrow(hcat(lonV,latV))
    getDist(unLonLat_points,lonlat)= Haversine().(unLonLat_points,[lonlat])
    idxV=fill(0,size(lonV))
    for i in eachindex(idxV)
        idxV[i]=argmin(getDist(unLonLat_points,lonLatV_points[i]))
    end
    pollV=unPoll[idxV,:]
    return pollV
end

function getNearest(lonM::AbstractArray,latM::AbstractArray,lon::AbstractVector,lat::AbstractVector;outType="lonlat",disp=true)
    output_all= [getNearest(lonM,latM,lon[n],lat[n];outType=outType,disp=disp) for n in eachindex(lon)]
    if lowercase.(outType)=="lonlat"
        out_lon=[output_all[n][1] for n in eachindex(output_all)]
        out_lat=[output_all[n][2] for n in eachindex(output_all)]
        return out_lon,out_lat
    elseif uppercase.(outType)=="IJ"
        out_I=[output_all[n][1] for n in eachindex(output_all)]
        out_J=[output_all[n][2] for n in eachindex(output_all)]
        return out_I,out_J
    elseif lowercase.(outType) in ["ind","linearind"]
        return output_all
    else
        error("outType should be lonlat or IJ or ind")
    end
end
function getNearest(lonM::AbstractArray,latM::AbstractArray,lon::Real,lat::Real;outType="lonlat",disp=true)
    I,J=getNearestIJ(lonM,latM,lon,lat;disp=disp)
    if lowercase.(outType)=="lonlat"
        return lonM[I,J],latM[I,J]
    elseif uppercase.(outType)=="IJ"
        return I,J
    elseif lowercase.(outType)=="ind"
        return CartesianIndex(I,J)
    elseif lowercase.(outType)=="linearind"
        return I + (J-1)*size(lonM,1)
    else
        error("outType should be lonlat or IJ or ind")
    end
end
function getNearestLonLat(lonM::AbstractArray,latM::AbstractArray,lon,lat)
    I,J=getNearestIJ(lonM,latM,lon,lat)
    ind=CartesianIndex.(I, J)
    return lonM[ind],latM[ind]
end

# 1. 提供一个内联的纯标量距离计算函数，完全剥离数组分配
# 返回相对距离（单调递增即可用于比较），省略了乘以地球半径，进一步提升速度
@inline function fast_haversine_metric(lon1::Real, lat1::Real, lon2::Real, lat2::Real)
    rad = π / 180.0
    φ1 = lat1 * rad
    φ2 = lat2 * rad
    Δφ = (lat2 - lat1) * rad
    Δλ = (lon2 - lon1) * rad
    
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    # 只要比较大小，直接返回 a 即可（因为 a 与距离正相关，省去 atan 和 sqrt 操作！）
    return a 
end

function getNearestIJ(lonM::AbstractArray, latM::AbstractArray, lon::Real, lat::Real; disp=true)
    size1, size2 = size(lonM)
    i = size1 ÷ 2 + 1
    j = size2 ÷ 2 + 1
    i_pre, j_pre = 0, 0
    
    while !(i == i_pre && j == j_pre)
        i_pre, j_pre = i, j
        
        # 1. 显式 for 循环搜索最佳的 i (固定 j)
        min_dist_i = Inf
        best_i = i
        for k in axes(lonM, 1)
            dist = fast_haversine_metric(lonM[k, j], latM[k, j], lon, lat)
            if dist < min_dist_i
                min_dist_i = dist
                best_i = k
            end
        end
        i = best_i
        
        # 2. 显式 for 循环搜索最佳的 j (固定 i)
        min_dist_j = Inf
        best_j = j
        for k in axes(lonM, 2)
            dist = fast_haversine_metric(lonM[i, k], latM[i, k], lon, lat)
            if dist < min_dist_j
                min_dist_j = dist
                best_j = k
            end
        end
        j = best_j
    end
    
    if disp && (i == 1 || i == size1 || j == 1 || j == size2)
        println("boundary sta occur for lon $lon lat $lat")
    end
    
    return i, j
end
# 内存消耗较高，且速度慢，已弃
# function getNearestIJ(lonM::AbstractArray,latM::AbstractArray,lon::Real,lat::Real;disp=true)
#     # only support numeric lon, lat
#     sizeGrid=size(lonM)
#     i,j=Int.(round.(sizeGrid ./2))
#     i_pre,j_pre=0,0
#     getDist(lonV,latV,lon,lat)= Haversine().(eachrow(hcat(lonV,latV)),[[lon,lat]])
#     while !(i==i_pre && j==j_pre)
#         i_pre,j_pre=i,j
#         i=argmin(getDist(lonM[:,j],latM[:,j],lon,lat))
#         j=argmin(getDist(lonM[i,:],latM[i,:],lon,lat))
#     end
#     if disp && ((i==0 | i==sizeGrid[1]) .| (j==0 | j==sizeGrid[2]))
#         println("boundary sta occur for lon $lon lat $lat")
#     end
#     return i,j
# end

function getNearestIJ(lonM::AbstractArray, latM::AbstractArray, lonV::AbstractVector, latV::AbstractVector; disp=true, par=false)
    # 修正点：显式创建标准的 Vector{Int}，而不是保留 ChainedVector 结构
    # 这样 n 就会退化为普通的 Int 索引，不会触发 ChainedVectorIndex 的校验
    I = Vector{Int}(undef, length(lonV))
    J = Vector{Int}(undef, length(latV))
    
    if par
        Threads.@threads for n in 1:length(I)
            # 使用整数 n 进行索引是安全的
            I[n], J[n] = getNearestIJ(lonM, latM, lonV[n], latV[n]; disp=disp)
        end
    else
        for n in 1:length(I)
            I[n], J[n] = getNearestIJ(lonM, latM, lonV[n], latV[n]; disp=disp)
        end
    end
    
    return I, J
end