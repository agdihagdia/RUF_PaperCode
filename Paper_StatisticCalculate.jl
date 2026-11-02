using StaDataM,Dates,TimeZones
using DataFrames, Statistics
using utilities: getStatistics

EVA_TIME = ZonedDateTime.(DateTime(2025,4,9):Hour(1):DateTime(2025,4,21),tz"Asia/Shanghai")
RUF_LEADS = collect(1:12)
BASE_LEADS = collect(1:24) # 预报时效，不是数组下标；不包含0 h。
Strat_RUF = collect((first(EVA_TIME)-Hour(maximum(RUF_LEADS))):Hour(1):(last(EVA_TIME)-Hour(minimum(RUF_LEADS))))
Strat_Base = ZonedDateTime.(DateTime(2025,4,8,12):Day(1):DateTime(2025,4,20,12),tz"UTC")

OBS_FILEPATH = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
RUF_ALLStation_FileFunction(t) = "/home/zsm/IETM/RUF_DryDeposition/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"
Base_FileFunction(t) = "/media/dataLL/zsm/paper_data/baseFolder/$(Dates.format(t,"yyyymmddHH"))/staData.$(Dates.format(t,"yyyymmddHH")).mat"
BaseFilePath = Base_FileFunction.(Strat_Base)
RUF_FilePath = RUF_ALLStation_FileFunction.(Strat_RUF)

# 第1维均为时间，第2维均为同一顺序的NationIDs。
# 以有效预报时间匹配，禁止通过截取相同长度来“对齐”。
function align_time_rows(values, source_times, target_times; label="forecast")
    size(values, 1) == length(source_times) || error("$label: 时间与数据行数不一致")
    length(unique(source_times)) == length(source_times) || error("$label: 有效时间重复")
    time_index = Dict(t => i for (i, t) in enumerate(source_times))
    absent = filter(t -> !haskey(time_index, t), target_times)
    isempty(absent) || error("$label: 缺少$(length(absent))个评估时刻，首个为$(first(absent))")
    rows = [time_index[t] for t in target_times]
    @assert source_times[rows] == target_times
    return values[rows, :]
end

#OBS
println("读取观测；评估范围：", first(EVA_TIME), " — ", last(EVA_TIME))
obs_staDATA = staDATA(OBS_FILEPATH)
staInfo = obs_staDATA.staInfo
NationInfo = staInfo[map(x -> isequal(x, 'B') || isequal(x, "B"), staInfo.type), :]
NationIDs = NationInfo.staID
isempty(NationIDs) && error("未找到国控站点")
length(unique(NationIDs)) == length(NationIDs) || error("国控站点ID重复")
obs = read(obs_staDATA; Times=EVA_TIME, varNames="PM10", staIDs=NationIDs, leadTimes=[0])[:,1,:,1]

#Base data
all(isfile, BaseFilePath) || error("部分Base起报文件不存在")
println("读取Base：", length(BaseFilePath), "个逐日起报文件，时效1—24 h")
BaseForecast_staDATA = staDATA(BaseFilePath; Times=Strat_Base,
    varNames=["PM10"], staInfo=NationInfo, leadTimes=BASE_LEADS, disp=true)
BaseForecast = read(BaseForecast_staDATA; Times=Strat_Base, varNames="PM10",
    staIDs=NationIDs, leadTimes=BASE_LEADS)[:,1,:,:]
# 获取数组维度
Ntime, Nsta, Nlead = size(BaseForecast)
# 转换为嵌套数组
BaseForecast_nested = [
    [collect(@view BaseForecast[i, j, :]) for j in 1:Nsta]
    for i in 1:Ntime
]
BaseForecast_long = [
    vcat((BaseForecast_nested[i][j] for i in 1:Ntime)...)
    for j in 1:Nsta
]
# Base目录名是UTC起报时间，转为北京时间后加预报时效。
# 每天20时使用前一批次24 h，21时才开始使用新批次1 h。
Base_valid_times = vcat([
    astimezone(t, tz"Asia/Shanghai") .+ Hour.(BASE_LEADS)
    for t in Strat_Base
]...)
Base_aligned = align_time_rows(
    reduce(hcat, BaseForecast_long), Base_valid_times, EVA_TIME; label="Base")

#RUF
all(isfile, RUF_FilePath) || error("部分RUF起报文件不存在")
println("读取RUF：", length(RUF_FilePath), "个逐小时起报文件")
RUF_staDATA = staDATA(RUF_FilePath; Times=Strat_RUF,
    varNames=["PM10"], staInfo=NationInfo, leadTimes=RUF_LEADS, disp=true)
RUF = read(RUF_staDATA; Times=Strat_RUF, varNames="PM10",
    staIDs=NationIDs, leadTimes=RUF_LEADS)[:,1,:,:]

# 固定时效L：有效时间t对应起报时间t-L，而不是直接拿相同起报时间比较。
RUF_aligned = Array{Float64}(undef, length(EVA_TIME), length(NationIDs), length(RUF_LEADS))
for (i, lead) in enumerate(RUF_LEADS)
    RUF_valid_times = Strat_RUF .+ Hour(lead)
    RUF_aligned[:, :, i] = align_time_rows(
        RUF[:, :, i], RUF_valid_times, EVA_TIME; label="RUF $(lead) h")
end
@assert size(obs) == size(Base_aligned) == size(RUF_aligned)[1:2]

# 所有时效和两种预报使用完全相同的有效“时刻—站点”样本。
# 只剔除缺测/非有限值，不更改浓度值，也不增加浓度阈值筛选。
valid_common = isfinite.(obs) .& isfinite.(Base_aligned)
for i in eachindex(RUF_LEADS)
    valid_common .&= isfinite.(@view RUF_aligned[:, :, i])
end
N_valid = count(valid_common)
N_valid >= 2 || error("共同有效样本不足，无法计算统计参数")
println("对齐完成：", length(EVA_TIME), "小时；",
    count(vec(any(valid_common; dims=1))), "个有共同有效样本的国控站；", N_valid, "个样本")

# 展平成共同有效的站点—小时向量，计算总体统计；矩阵直接传入会变成逐站统计。
# getStatistics(reference, predict)：观测在前，bias为预报减观测。
obs_curr = obs[valid_common]
Base_curr = Base_aligned[valid_common]
RUF_statistics = Dict{Int,NamedTuple}()
Base_statistics = Dict{Int,NamedTuple}()
statistics_rows = NamedTuple[]
for (i, lead) in enumerate(RUF_LEADS)
    RUF_curr = RUF_aligned[:, :, i][valid_common]
    @assert length(RUF_curr) == length(Base_curr) == length(obs_curr)
    RUF_statistics[lead] = getStatistics(obs_curr, RUF_curr)
    Base_statistics[lead] = getStatistics(obs_curr, Base_curr)
    ruf_stat, base_stat = RUF_statistics[lead], Base_statistics[lead]
    reduction = base_stat.rmse == 0 ? NaN : 100 * (1 - ruf_stat.rmse / base_stat.rmse)
    push!(statistics_rows, (;
        lead_h=lead, n_pairs=N_valid, rmse_reduction_pct=reduction,
        (Symbol("RUF_", k) => v for (k, v) in pairs(ruf_stat))...,
        (Symbol("Base_", k) => v for (k, v) in pairs(base_stat))...,
    ))
end

# 保留全部9项统计参数；可在REPL中查看Statistics_Result和上面的两个字典。
Statistics_Result = DataFrame(statistics_rows)
Mean_RMSE_reduction_pct = mean(Statistics_Result.rmse_reduction_pct)
show(stdout, MIME"text/plain"(), Statistics_Result[:, [
    :lead_h, :n_pairs, :RUF_rmse, :Base_rmse, :RUF_corr, :Base_corr, :rmse_reduction_pct,
]])
println("\n前", length(RUF_LEADS), "个预报时效的平均RMSE降幅：", Mean_RMSE_reduction_pct, "%")
