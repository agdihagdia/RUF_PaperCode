using Distances
export obsMapper,map2sta,getStaInDomain, obsNearestMapper
import Interpolations as Itp
import ScatteredInterpolation as SItp

include("getNearest.jl")
# Sample Usage:
#     using MAT
#     vars = matread("lonlatM.mat")
#     lonM,latM=(vars[k] for k in ["lonM","latM"])
#     yM=@. (lonM-latM).^2
#     v=matread("staInfoNation.mat")
#     staInfo=v["staInfo"]
#     staLon=staInfo["lon"]
#     staLat=staInfo["lat"]
#     mp=obsMapper(lonM,latM,vec(staLon),vec(staLat))
#     y=map(mp,yM)
#     # using BenchmarkTools
#     # @benchmark y=map(mp,yM)
#     using PyPlot
#     #using Plots; pyplot()
#     pcolor(lonM,latM,yM, cmap="jet",vmin=2000,vmax=8000)
#     scatter(vec(staLon),vec(staLat),s=6,c=vec(y),cmap="jet",vmin=2000,vmax=8000)
#     colorbar()
#     gcf()

# updates:
#     2022-09-29 add getIJ_ingrid to cancualte the nearest IJ without using PyCall
#     2022-09-29 update obsMapper eliminate the dependence of PyCall which might call error
#     2022-11-29 add functions getNearestLonLat
#     2022-11-29 add support for vetor lonV latV in getNearestLonLat
#     2022-11-29 add support for method "toGrid" in obsMapper
#     2022-11-29 add a fied method in struct obsMapper
#     2025-08-21 change interpolation method to inverse distance weight in obsMapper "toSta" to avoid strange error. SItp.Shepard()


struct obsNearestMapper
    lonM::Matrix
    latM::Matrix
    staLon::AbstractVector
    staLat::AbstractVector
    staIdx::AbstractVector
end
function obsNearestMapper(lonM, latM,staLon, staLat;dispProgress=false,disp=true,par=true)
    # getNearest(lonM::AbstractArray,latM::AbstractArray,lon::AbstractVector,lat::AbstractVector;outType="lonlat",disp=true)
    if par
        ind_allChunk=Vector{Any}(undef,length(staLon))
        if dispProgress; p = Progress(length(ind_allChunk)); end
        Threads.@threads for j in eachindex(staLon)
            if dispProgress; ProgressMeter.next!(p); end
            ind_allChunk[j]=getNearest(lonM,latM,staLon[j],staLat[j];outType="ind",disp=disp)
        end
        staIdx=vcat(ind_allChunk...)
    else
        staIdx = getNearest(lonM,latM,staLon,staLat;outType="ind",disp=disp)
    end
    return obsNearestMapper(lonM, latM,staLon, staLat,  staIdx)
end

function map2sta(mapper::obsNearestMapper,gridData)
    return gridData[mapper.staIdx]
end

struct obsMapper
    lonM::Matrix
    latM::Matrix
    staLon::AbstractVector
    staLat::AbstractVector
    staX::AbstractVector
    staY::AbstractVector
    method::String
end

function getNearestLonLat(mapper::obsMapper)
    if mapper.method == "toGrid"
        inds=CartesianIndex.(mapper.staX,mapper.staY)
        lonV=mapper.lonM[inds]
        latV=mapper.latM[inds]
    else # toSta
        lonV,latV=getNearestLonLat(mapper.lonM,mapper.latM,mapper.staLon,mapper.staLat)
    end
    return lonV,latV
end
function obsMapper(lonM,latM,staLon::AbstractVector,staLat::AbstractVector;method="toSta",par=false)
    sizeGrid=size(lonM)

    if method=="toSta"
        staX,staY=getNearestIJ(lonM,latM,staLon,staLat;disp=false,par=par) # checked pass
        staX_float,staY_float=Float64.(staX),Float64.(staY)
        # function calXY_aSta(lonM,latM,lon,lat,intX::Int,intY::Int,sizeGrid)
        #     i,j=intX,intY
        #     ij_around_all=[[ii,jj] for ii=max(1,i-1):min(i+1,sizeGrid[1]) for jj=max(1,j-1):min(j+1,sizeGrid[2]) ]
        #     ij_around_all=hcat(ij_around_all...)'
        #     lon_around=[lonM[ij[1],ij[2]] for ij in eachrow(ij_around_all)]
        #     lat_around=[latM[ij[1],ij[2]] for ij in eachrow(ij_around_all)]
        #     lonlat_around_all=hcat(lon_around,lat_around)
        #     xy_float=SItp.evaluate( 
        #         SItp.interpolate(SItp.InverseMultiquadratic(), Float64.(lonlat_around_all'), Float64.(ij_around_all) )
        #           ,[lon;lat])
        #     return xy_float
        # end
        # 局部平面最小二乘法 (Local LS Fitting) 进行插值
        # 使用最小二乘法拟合 3x3 邻域，求取平滑的 u, v
        function calXY_LS_smooth(lonM, latM, target_lon, target_lat, intI, intJ, sizeGrid)
            # 1. 准备 3x3 邻域数据，使用标量累加避免分配
            sum_lon = 0.0; sum_lat = 0.0
            sum_i = 0.0;   sum_j = 0.0
            sum_ii = 0.0;  sum_jj = 0.0; sum_ij = 0.0
            sum_i_lon = 0.0; sum_j_lon = 0.0; sum_i_lat = 0.0; sum_j_lat = 0.0
            
            count = 0
            for dj in -1:1, di in -1:1
                ii, jj = intI + di, intJ + dj
                # 边界检查
                if 1 <= ii <= sizeGrid[1] && 1 <= jj <= sizeGrid[2]
                    l_on, l_at = lonM[ii, jj], latM[ii, jj]
                    # 这里以 intI, intJ 为中心计算偏移，增加数值稳定性
                    x, y = Float64(di), Float64(dj)
                    
                    sum_lon += l_on; sum_lat += l_at
                    sum_i += x; sum_j += y
                    sum_ii += x*x; sum_jj += y*y; sum_ij += x*y
                    sum_i_lon += x*l_on; sum_j_lon += y*l_on
                    sum_i_lat += x*l_at; sum_j_lat += y*l_at
                    count += 1
                end
            end

            # 2. 构建正规方程组求解局部梯度 (雅可比矩阵)
            # 这是一个简单的 3x3 矩阵运算，求 [∂lon/∂i, ∂lon/∂j, lon_offset]
            # 以及 [∂lat/∂i, ∂lat/∂j, lat_offset]
            # 为了性能，直接硬编码 Cramer 法则求解 2x2 子系统（忽略偏移项的耦合）
            
            # 简化版：中心化后的 LS 求解
            det = sum_ii * sum_jj - sum_ij^2
            if abs(det) < 1e-12
                return [Float64(intI), Float64(intJ)] # 退化情况
            end

            # 计算局部梯度 (Jacobian)
            dLon_di = (sum_i_lon * sum_jj - sum_ij * sum_j_lon) / det
            dLon_dj = (sum_ii * sum_j_lon - sum_i_lon * sum_ij) / det
            dLat_di = (sum_i_lat * sum_jj - sum_ij * sum_j_lat) / det
            dLat_dj = (sum_ii * sum_j_lat - sum_i_lat * sum_ij) / det
            
            # 计算中心点残差
            lon_0 = sum_lon / count
            lat_0 = sum_lat / count
            
            # 3. 求解最终的浮点偏移 (u, v)
            # 方程组: 
            # dLon_di * u + dLon_dj * v = target_lon - lon_0
            # dLat_di * u + dLat_dj * v = target_lat - lat_0
            b1 = target_lon - lon_0
            b2 = target_lat - lat_0
            
            det_final = dLon_di * dLat_dj - dLon_dj * dLat_di
            u = (b1 * dLat_dj - dLon_dj * b2) / det_final
            v = (dLon_di * b2 - b1 * dLat_di) / det_final
            
            # 4. 返回结果
            # 因为是以 intI 为中心拟合，结果需要加回 intI
            return [intI + u, intJ + v]
        end
        intStaX,intStaY=Int.(staX),Int.(staY)
        if par 
            Threads.@threads for n_sta in eachindex(staLon)
                staX_float[n_sta],staY_float[n_sta]=calXY_LS_smooth(lonM,latM,staLon[n_sta],staLat[n_sta],intStaX[n_sta],intStaY[n_sta],sizeGrid)
            end
        else
            for n_sta in eachindex(staLon)
                staX_float[n_sta],staY_float[n_sta]=calXY_LS_smooth(lonM,latM,staLon[n_sta],staLat[n_sta],intStaX[n_sta],intStaY[n_sta],sizeGrid)
            end
        end
        staX,staY=staX_float,staY_float
    elseif method == "toGrid"
        staX,staY=getNearestIJ(lonM,latM,staLon,staLat)
    end
    return obsMapper(lonM,latM,staLon,staLat,staX,staY,method)
end
function getNearestData(gridLon,gridLat,meshData,lonV,latV)
    getInd(xBase,x) =argmin(abs.(xBase.-x))
    I=(x->getInd(gridLon,x)).(lonV)
    J=(x->getInd(gridLat,x)).(latV)
    outV=[meshData[i,j] for (i,j) in zip(I,J)]
    return outV
end
function getStaInDomain(lonM,latM,staLon::AbstractVector,staLat::AbstractVector)
    
    staIdx = getNearest(lonM,latM,staLon,staLat;outType="ind",disp=false)
    isInnerM=trues(size(lonM))
    isInnerM[:,[1,end]] .= false
    isInnerM[[1,end],:] .= false
    q_inDomain=isInnerM[staIdx]
    return q_inDomain
end

function map2sta(mapper::obsMapper,gridData)
    if mapper.method == "toGrid"
        outData=gridData[CartesianIndex.(mapper.staX,mapper.staY)]
    else # toSta
        itp = Itp.extrapolate(Itp.interpolate(gridData, Itp.BSpline(Itp.Linear())),NaN)
        outData=itp.(mapper.staX,mapper.staY)
    end
    return outData
end
