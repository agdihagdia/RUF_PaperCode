using CSV
using Dates
using DataFrames
using GradDataM
using JLD2
using NaNStatistics
using Serialization
using StaDataM
using Statistics
using TimeZones

const COMPARE_OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const COMPARE_NATION_ROOT = "/home/zsm/IETM/RUF_DryDeposition_Nation/Ruf_forcast"
const COMPARE_ALL_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const COMPARE_DATAGRID = "/home/zsm/IETM/RUF_2d/datagrid.jld2"
const COMPARE_METRICS = joinpath(get(ENV, "RUF_PAPER_STATS_DIR", "/home/zsm/IETM/code/RUF_PaperStatisticCalculation/output"), "figure3_network_city_metrics.csv")
const COMPARE_TS = joinpath(get(ENV, "RUF_PAPER_STATS_DIR", "/home/zsm/IETM/code/RUF_PaperStatisticCalculation/output"), "figure3_network_snapshot_ts_lead2.csv")
const COMPARE_CACHE = "/home/zsm/IETM/figure/PapreFigure/20260817/figure_compare_NationAndAllsta_prepared_cache_20260817.bin"
const COMPARE_BJT = tz"Asia/Shanghai"
const COMPARE_STAMP = dateformat"yyyymmddHH"

compare_clean(value) = value === missing ? String() : value isa AbstractVector ? join(string.(value)) : strip(string(value))
compare_city(value) = compare_clean(value)
compare_stamp(value::DateTime) = Dates.format(value, COMPARE_STAMP)

function compare_finite_mean(values)
    finite_values = Float64.(values[isfinite.(values)])
    return isempty(finite_values) ? NaN : mean(finite_values)
end

println("Loading observation and station metadata for Nation/All-stations comparison")
staObs = staDATA(COMPARE_OBS_FILE)
template_path = joinpath(COMPARE_ALL_ROOT, "2025041108", "staData.2025041108.mat")
template = staDATA(template_path)
full_station_ids = compare_clean.(template.staInfo.staID)
full_city_labels = compare_city.(template.staInfo.city)
valid_station_indices = findall(!isempty, full_city_labels)
staInfo = template.staInfo[valid_station_indices, :]
city_labels = full_city_labels[valid_station_indices]
sta_lon = Float64.(staInfo.lon)
sta_lat = Float64.(staInfo.lat)

cities = ["延安", "太原", "大同", "忻州"]
city_groups = [findall(==(city), city_labels) for city in cities]
any(isempty, city_groups) && error("one or more Figure 3 cities have no forecast stations")

tStart_obs = ZonedDateTime.(DateTime(2025, 4, 11):Hour(1):DateTime(2025, 4, 13), COMPARE_BJT)
valid_times = DateTime.(tStart_obs)
snapshot_times = [DateTime(2025, 4, 11, hour) for hour in (10, 12, 14, 16)]

observation_ids = compare_clean.(staObs.staInfo.staID)
observation_position = Dict(id => index for (index, id) in enumerate(observation_ids))
selected_station_ids = full_station_ids[valid_station_indices]
all(haskey(observation_position, id) for id in selected_station_ids) ||
    error("forecast station IDs are absent from the observation file")
selected_observation_positions = [observation_position[id] for id in selected_station_ids]
pm10_position = findfirst(==("PM10"), uppercase.(compare_clean.(staObs.varNames)))
pm10_position === nothing && error("PM10 is absent from the observation file")
observation_time_lookup = Dict(
    DateTime(astimezone(time, COMPARE_BJT)) => index for (index, time) in enumerate(staObs.times)
)
observation_indices = [observation_time_lookup[time] for time in valid_times]
observation_field = Float64.(
    staObs.staDataAll[observation_indices, pm10_position, selected_observation_positions, 1]
)
observation_field[
    .!isfinite.(observation_field) .| (observation_field .< 0) .| (observation_field .>= 2000)
] .= NaN

obs = Vector{Vector{Float64}}(undef, length(cities))
for (city_index, positions) in enumerate(city_groups)
    obs[city_index] = vec(nanmean(@view(observation_field[:, positions]), dims=2))
end
obs_data2plot = [
    copy(vec(@view observation_field[findfirst(==(time), valid_times), :])) for time in snapshot_times
]
lonM, latM = load(COMPARE_DATAGRID, "lonM", "latM")

if isfile(COMPARE_CACHE)
    println("Loading cached Nation/All-stations lead-2 series and grids")
    cached = deserialize(COMPARE_CACHE)
    ietm3d = cached.ietm3d
    ietm3d_allsta = cached.ietm3d_allsta
    Nation_data2plot = cached.Nation_data2plot
    Allsta_data2plot = cached.Allsta_data2plot
    missing_cycles = cached.missing_cycles
else
    println("Building Nation/All-stations lead-2 series from single forecast files")
    ietm3d = [fill(NaN, length(valid_times)) for _ in cities, _ in 1:72]
    ietm3d_allsta = [fill(NaN, length(valid_times)) for _ in cities, _ in 1:72]
    missing_cycles = NamedTuple[]
    for (time_index, valid_time) in enumerate(valid_times)
        cycle = valid_time - Hour(2)
        cycle_text = compare_stamp(cycle)
        for (experiment, root, target) in (
            ("Nation-only", COMPARE_NATION_ROOT, ietm3d),
            ("All-stations", COMPARE_ALL_ROOT, ietm3d_allsta),
        )
            path = joinpath(root, cycle_text, "staData.$cycle_text.mat")
            if !isfile(path)
                push!(missing_cycles, (experiment=experiment, cycle_bjt=cycle_text))
                continue
            end
            object = staDATA(path)
            compare_clean.(object.staInfo.staID) == full_station_ids ||
                error("$experiment station ordering changed at $cycle_text")
            lead_position = findfirst(==(2), Int.(object.leadTimes))
            lead_position === nothing && error("lead 2 absent at $experiment/$cycle_text")
            field = Float64.(object.staDataAll[1, 1, valid_station_indices, lead_position])
            field[.!isfinite.(field) .| (field .< 0)] .= NaN
            for (city_index, positions) in enumerate(city_groups)
                target[city_index, 2][time_index] = compare_finite_mean(@view field[positions])
            end
        end
        (time_index == 1 || time_index % 12 == 0 || time_index == length(valid_times)) &&
            println("  processed valid time $time_index/$(length(valid_times)): $(compare_stamp(valid_time))")
    end

    Nation_data2plot = Vector{Matrix{Float32}}(undef, length(snapshot_times))
    Allsta_data2plot = Vector{Matrix{Float32}}(undef, length(snapshot_times))
    println("Reading lead-2 model grids for four snapshots")
    for (time_index, valid_time) in enumerate(snapshot_times)
        cycle = valid_time - Hour(2)
        cycle_text = compare_stamp(cycle)
        for (root, target) in (
            (COMPARE_NATION_ROOT, Nation_data2plot),
            (COMPARE_ALL_ROOT, Allsta_data2plot),
        )
            ctl = joinpath(root, cycle_text, "ruf.ctl")
            isfile(ctl) || error("missing grid control file: $ctl")
            grid = gradData(ctl)
            setTime!(grid, ZonedDateTime(valid_time, COMPARE_BJT))
            target[time_index] = Float32.(read(grid, "PM10", 1))
        end
    end
    serialize(COMPARE_CACHE, (
        ietm3d=ietm3d,
        ietm3d_allsta=ietm3d_allsta,
        Nation_data2plot=Nation_data2plot,
        Allsta_data2plot=Allsta_data2plot,
        missing_cycles=missing_cycles,
    ))
    println("Saved comparison cache: $COMPARE_CACHE")
end

metric_table = CSV.read(COMPARE_METRICS, DataFrame)
function compare_metric(experiment, city, column)
    selected = metric_table[(metric_table.experiment .== experiment) .& (metric_table.city .== city), :]
    nrow(selected) == 1 || error("missing audited metric for $experiment/$city")
    return Float64(selected[1, column])
end
Nation_city_rmse = [compare_metric("Nation-only", city, :rmse) for city in cities]
Nation_city_r = [compare_metric("Nation-only", city, :r) for city in cities]
All_city_rmse = [compare_metric("All-stations", city, :rmse) for city in cities]
All_city_r = [compare_metric("All-stations", city, :r) for city in cities]

ts_table = CSV.read(COMPARE_TS, DataFrame)
function compare_ts(experiment, valid_time)
    stamp = parse(Int, compare_stamp(valid_time))
    selected = ts_table[(ts_table.experiment .== experiment) .& (ts_table.valid_time_bjt .== stamp), :]
    nrow(selected) == 1 || error("missing audited TS for $experiment/$stamp")
    return Float64(selected.ts[1])
end
Nation_TS = [compare_ts("Nation-only", time) for time in snapshot_times]
All_TS = [compare_ts("All-stations", time) for time in snapshot_times]

println(
    "Nation/All-stations preparation complete; missing station cycles: ",
    length(missing_cycles),
)
