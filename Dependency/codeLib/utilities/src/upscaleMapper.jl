using ProgressMeter
export upscaleMapper,map2coarse
struct upscaleMapper
    lonM_fine::AbstractArray
    latM_fine::AbstractArray
    lonM_coarse::AbstractArray
    latM_coarse::AbstractArray
    idx_fineInCoarse::AbstractArray
    method::String
end
function upscaleMapper(lonM_fine, latM_fine, lonM_coarse, latM_coarse; method="avearea",dispProgress=true,disp=true,par=true)
    # getNearest(lonM::AbstractArray,latM::AbstractArray,lon::AbstractVector,lat::AbstractVector;outType="lonlat",disp=true)
    if par
        ind_allChunk=Vector{Any}(undef,size(lonM_fine,2))
        if dispProgress; p = Progress(length(ind_allChunk)); end
        Threads.@threads for j in axes(lonM_fine,2)
            if dispProgress; ProgressMeter.next!(p); end
            ind_allChunk[j]=getNearest(lonM_coarse,latM_coarse,lonM_fine[:,j],latM_fine[:,j];outType="ind",disp=disp)
        end
        idx_fineInCoarse=hcat(ind_allChunk...)
    else
        idx_fineInCoarse = getNearest(lonM_coarse,latM_coarse,lonM_fine,latM_fine;outType="ind",disp=disp)
    end
    return upscaleMapper(lonM_fine, latM_fine, lonM_coarse, latM_coarse, idx_fineInCoarse, method)
end

# 1. 内核函数：专门处理 CartesianIndex，类型完全确定
function _map_kernel!(data_coarse::Array{Float64, N}, counts::Array{Int, N}, 
    idx_fineInCoarse::AbstractVector{CartesianIndex{N}}, data_fine::AbstractVector{Float64}) where N

    # 遍历每一个精细网格点
    # @inbounds 移除边界检查（确保 idx_fineInCoarse 不会越界）
    # @simd 允许编译器进行向量化指令优化
    @inbounds @simd for i in eachindex(idx_fineInCoarse)
        idx_coarse = idx_fineInCoarse[i] # 这里编译器明确知道它是 CartesianIndex{2}
        val = data_fine[i]
        if !isnan(val)
            # 直接使用 CartesianIndex 索引，Julia对此优化极佳
            data_coarse[idx_coarse] += val
            counts[idx_coarse] += 1
        end
    end

    # 归一化步骤：融合循环，避免额外分配
    @inbounds for i in eachindex(data_coarse)
        c = counts[i]
        if c > 0
            data_coarse[i] /= c
        else
            data_coarse[i] = NaN
        end
    end

    return data_coarse
end

# 2. 外部接口：建立“函数屏障”，让编译器在传入内核前识别类型
function map2coarse(mapper::upscaleMapper, data_fine)
    # 提取数组，即便 mapper 结构体类型不稳定，
    # 传入 _map_kernel! 后，参数就会变为具体的强类型
    idx_fineInCoarse = mapper.idx_fineInCoarse
    sizeOut = size(mapper.lonM_coarse)

    # 预分配内存（复用 data_coarse 存储累加和，最后转为平均值）
    data_coarse = zeros(Float64, sizeOut)
    counts = zeros(Int, sizeOut)

    # 调用内核
    return _map_kernel!(data_coarse, counts, vec(idx_fineInCoarse), Float64.(vec(data_fine)))
end
