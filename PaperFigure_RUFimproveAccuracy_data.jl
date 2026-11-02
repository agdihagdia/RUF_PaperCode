# Compatibility layer for reusing errortable_zsm3.jl with the current
# staDATAStru files.  The current StaDataM multi-file constructor is not
# compatible with times::Vector{ZonedDateTime}, so cycles are read one by one.

const NEW_OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const NEW_FC_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const NEW_STATS_FILE = joinpath(get(ENV, "RUF_PAPER_STATS_DIR", "/home/zsm/IETM/code/RUF_PaperStatisticCalculation/output"), "figure6_city_metrics.csv")

finite_mean_new(values) = begin
    selected = Float64.(values[isfinite.(values)])
    isempty(selected) ? NaN : mean(selected)
end

staObs = staDATA(NEW_OBS_FILE)
staInfo = staObs.staInfo
# MATLAB cells with missing city metadata are decoded as empty String vectors.
# Keep only scalar, non-empty city names so the original broadcast comparisons work.
cities = [city for city in unique(staInfo.city) if city isa AbstractString && !isempty(city)]
var2plot = "PM10"
PlaCar = ccrs.PlateCarree()

leadtimes = collect(1:72)
leadtimes24 = collect(1:24)
leadtimes2plot = [1, 6, 12, 24]
cities2plot = ["太原", "北京", "南京", "海口"]
tStart_all_obs = ZonedDateTime.(DateTime(2025,4,11):Hour(1):DateTime(2025,4,21), tz"Asia/Shanghai")
tStart_One = ZonedDateTime(DateTime(2025,4,11), tz"Asia/Shanghai")

println("Loading observations for errortable_zsm3 reuse")
obs = Vector{Vector{Float64}}(undef, length(cities))
for (city_index, city) in enumerate(cities)
    values = vec(nanmean(read(staObs, Times=tStart_all_obs, varNames="PM10", staIDs=get_staIDs(city)), dims=3))
    values = Float64.(values)
    values[.!isfinite.(values) .| (values .< 0) .| (values .>= 2000)] .= NaN
    obs[city_index] = values
end

println("Loading audited Figure 6 metrics")
metric_table = CSV.read(NEW_STATS_FILE, DataFrame)
RMSE_IETM = fill(NaN, length(cities), 24)
R_IETM = fill(NaN, length(cities), 24)
city_lookup = Dict(string(city) => index for (index, city) in enumerate(cities))
for row in eachrow(metric_table)
    city = string(row.city)
    frequency = Int(row.update_frequency_hour)
    haskey(city_lookup, city) || continue
    1 <= frequency <= 24 || continue
    RMSE_IETM[city_lookup[city], frequency] = Float64(row.rmse)
    R_IETM[city_lookup[city], frequency] = Float64(row.r)
end

IETM_Nlead = [Float64[] for _ in 1:length(cities), _ in 1:24]
all_start_times_format = [DateTime[] for _ in 1:length(cities), _ in 1:24]
total_len = length(tStart_all_obs)
fixed_start = first(tStart_all_obs) - Hour(1)

template_path = joinpath(NEW_FC_ROOT, "2025041118", "staData.2025041118.mat")
template = staDATA(template_path)
template_cities = [value isa AbstractVector ? join(string.(value)) : strip(string(value)) for value in template.staInfo.city]
positions_by_city = Dict(city => findall(==(city), template_cities) for city in cities2plot)

println("Building representative update-frequency series by single-file reads")
for frequency in leadtimes2plot
    series = Dict(city => fill(NaN, total_len) for city in cities2plot)
    starts = DateTime[]
    filled = 0
    block = 0
    while filled < total_len
        cycle = fixed_start + Hour(block * frequency)
        push!(starts, DateTime(cycle))
        take = min(frequency, total_len - filled)
        cycle_stamp = Dates.format(cycle, "yyyymmddHH")
        path = joinpath(NEW_FC_ROOT, cycle_stamp, "staData.$cycle_stamp.mat")
        if isfile(path)
            object = staDATA(path)
            for city in cities2plot
                positions = positions_by_city[city]
                for offset in 1:take
                    values = view(object.staDataAll, 1, 1, positions, offset)
                    series[city][filled + offset] = finite_mean_new(values)
                end
            end
        else
            println("  missing cycle retained as NaN: ", cycle_stamp)
        end
        filled += take
        block += 1
    end
    for city in cities2plot
        city_index = city_lookup[city]
        IETM_Nlead[city_index, frequency] = series[city]
        all_start_times_format[city_index, frequency] = copy(starts)
    end
    println("  completed U", frequency, "h")
end

println("New-data compatibility preparation complete")
