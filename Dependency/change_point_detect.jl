using StaDataM,Dates,TimeZones,NaNStatistics,Statistics
using Distributions
function smooth_data(data::Vector{Float64}, window_size::Int = 5)
    n = length(data)
    smoothed_data = zeros(Float64, n)
    
    half_window = div(window_size, 2)  # 窗口的半宽度

    for i in 1:n
        # 定义窗口范围
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        
        # 计算窗口内的均值
        window_data = data[start_idx:end_idx]
        smoothed_data[i] = nanmean(window_data)
    end
    
    return smoothed_data
end
function compute_derivative(data::Vector{Float64})
    derivative = fill(NaN, length(data))
    #右导数，不计算最后一点
    for i in 1:length(data) - 1
        derivative[i] = data[i+1] - data[i]
    end
    return derivative
end
function detect_starting_points(data::Vector{Float64}, derivative_threshold::Float64, concentration_threshold::Float64)
    data=smooth_data(data)
    derivative = compute_derivative(data)   
    valid_derivative = filter(!isnan, derivative)
    if isempty(valid_derivative)
        return Int64[], Int64[]  # 所有导数都是 NaN，直接返回
    end
    derivative_max = maximum(valid_derivative)
    derivative_min = minimum(valid_derivative)

    # 计算归一化导数
    # normal_der = (derivative .- derivative_min) ./ (derivative_max .- derivative_min + 1e-6)
    normal_der = (derivative ) ./ (derivative_max .- derivative_min + 1e-6)
    #计算整个序列浓度最大值
    data_max=maximum(data)
    if data_max > concentration_threshold
    # q_start= (normal_der .> derivative_threshold) .& (derivative .> 0) .& (data .> adaptive_threshold)
        q_start= (normal_der .> derivative_threshold) .& (derivative .> 0) 
    else
        return Int64[], normal_der
    end
    return findall(q_start),normal_der
end
function filter_starting_points(
    obs::AbstractVector,
    forecast::AbstractVector,                  # ← 修正拼写
    forecast_starts::AbstractVector{<:Integer},
    obs_starts::AbstractVector{<:Integer};
    search_radius::Int = 3
)
    # 空候选：直接返回 (-1, -1)
    if isempty(obs_starts) || isempty(forecast_starts)
        return (-1, -1)
    end

    # 半径下限为 1，避免 0 造成空切片
    r = max(1, min(search_radius, length(obs_starts), length(forecast_starts)))

    # 取前 r 个候选
    cand_obs = obs_starts[1:r]
    cand_fc  = forecast_starts[1:r]

    # 边界与 NaN 保护（NaN 当成 +Inf，不会被选中）
    @assert all(1 .<= cand_obs .<= length(obs)) "obs_starts contain out-of-bounds indices."
    @assert all(1 .<= cand_fc  .<= length(forecast)) "forecast_starts contain out-of-bounds indices."

    vals_obs = @view obs[cand_obs]
    vals_fc  = @view forecast[cand_fc]

    # 将 NaN 替换为 +Inf（只在比较时用，不改原数据）
    safe_obs = map(x -> isnan(x) ? Inf : x, vals_obs)
    safe_fc  = map(x -> isnan(x) ? Inf : x, vals_fc)

    k_obs = findmin(safe_obs)[2]   # 返回 (最小值, 位置) 的第 2 项是位置
    k_fc  = findmin(safe_fc)[2]

    return (obs_starts[k_obs], forecast_starts[k_fc])
end
function judge_forcast(obs,ietm,obs_point,ietm_point,tolerence)
    obs_firstpoint,forecast_firstpoint=filter_starting_points(obs,ietm,ietm_point,obs_point)
    if obs_firstpoint == -1 || forecast_firstpoint == -1
        return 0
    end
    if abs(obs_firstpoint - forecast_firstpoint) <= tolerence
        return 1
    else
        return 0
    end
end
function bool2leadtime(judge_array)
    judge_array_curr=copy(judge_array)
    rows,cols = size(judge_array_curr)
    for i in 1:rows
        for j in 1:cols
            if judge_array_curr[i,j] == 1
                judge_array_curr[i,j] = j
            else
                judge_array_curr[i,j] = 0
            end
        end
    end
    leadtimes_arr=maximum(judge_array_curr,dims=2)
    return leadtimes_arr
end
