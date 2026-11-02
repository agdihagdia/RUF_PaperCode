using Dates
using NaNStatistics
using Serialization
using StaDataM
using Statistics
using TimeZones

const FIG3_OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const FIG3_FC_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const FIG3_BJT = tz"Asia/Shanghai"
const FIG3_STAMP = dateformat"yyyymmddHH"
const FIG3_CACHE_FILE = "/home/zsm/IETM/figure/PapreFigure/20260817/figure3_v2_prepared_cache_20260817.bin"

figure3_clean(value) = value === missing ? String() : value isa AbstractVector ? join(string.(value)) : strip(string(value))
figure3_city(value) = figure3_clean(value)
figure3_stamp(value::DateTime) = Dates.format(value, FIG3_STAMP)

function figure3_finite_mean(values)
    finite_values = Float64.(values[isfinite.(values)])
    return isempty(finite_values) ? NaN : mean(finite_values)
end

println("Loading observations and station metadata for figure3_v2 reuse")
staObs = staDATA(FIG3_OBS_FILE)
template_path = joinpath(FIG3_FC_ROOT, "2025041118", "staData.2025041118.mat")
template = staDATA(template_path)

full_forecast_ids = figure3_clean.(template.staInfo.staID)
template_city_labels = figure3_city.(template.staInfo.city)
valid_station_indices = findall(!isempty, template_city_labels)
staInfo = template.staInfo[valid_station_indices, :]
city_labels = template_city_labels[valid_station_indices]
unique_cities = unique(city_labels)
cities2detect = unique_cities
city_lookup = Dict(city => index for (index, city) in enumerate(cities2detect))
city_groups = [findall(==(city), city_labels) for city in cities2detect]

var2plot = "PM10"
PlaCar = ccrs.PlateCarree()
tStart_all = ZonedDateTime.(DateTime(2025, 4, 8):Hour(1):DateTime(2025, 4, 21), FIG3_BJT)
tStart_all_obs = ZonedDateTime.(DateTime(2025, 4, 11):Hour(1):DateTime(2025, 4, 14), FIG3_BJT)
valid_times = DateTime.(tStart_all_obs)

observation_ids = figure3_clean.(staObs.staInfo.staID)
observation_position = Dict(id => index for (index, id) in enumerate(observation_ids))
selected_ids = full_forecast_ids[valid_station_indices]
all(haskey(observation_position, id) for id in selected_ids) ||
    error("forecast station IDs are absent from the observation file")
selected_observation_positions = [observation_position[id] for id in selected_ids]
pm10_position = findfirst(==("PM10"), uppercase.(figure3_clean.(staObs.varNames)))
pm10_position === nothing && error("PM10 is absent from the observation file")
observation_time_lookup = Dict(
    DateTime(astimezone(time, FIG3_BJT)) => index for (index, time) in enumerate(staObs.times)
)
observation_indices = [observation_time_lookup[time] for time in valid_times]
observation_field = Float64.(
    staObs.staDataAll[observation_indices, pm10_position, selected_observation_positions, 1]
)
observation_field[
    .!isfinite.(observation_field) .| (observation_field .< 0) .| (observation_field .>= 2000)
] .= NaN

obs_array = Vector{Vector{Float64}}(undef, length(cities2detect))
for (city_index, positions) in enumerate(city_groups)
    obs_array[city_index] = vec(nanmean(@view(observation_field[:, positions]), dims=2))
end

leadtimes = collect(1:72)
ietm3d_array = [fill(NaN, length(valid_times)) for _ in cities2detect, _ in leadtimes]
first_valid = first(valid_times)
last_valid = last(valid_times)
cycle_times = collect((first_valid - Hour(72)):Hour(1):(last_valid - Hour(1)))

loaded_cycles = 0
missing_cycles = String[]
if isfile(FIG3_CACHE_FILE)
    println("Loading cached 1-72 h city series")
    cached = deserialize(FIG3_CACHE_FILE)
    ietm3d_array = cached.ietm3d_array
    loaded_cycles = cached.loaded_cycles
    missing_cycles = cached.missing_cycles
else
println("Building 1-72 h city series from single forecast files")
for (cycle_index, cycle) in enumerate(cycle_times)
    cycle_text = figure3_stamp(cycle)
    path = joinpath(FIG3_FC_ROOT, cycle_text, "staData.$cycle_text.mat")
    if !isfile(path)
        push!(missing_cycles, cycle_text)
        println("  missing cycle retained as NaN: $cycle_text")
        continue
    end
    object = staDATA(path)
    figure3_clean.(object.staInfo.staID) == full_forecast_ids ||
        error("station ordering changed at $cycle_text")
    Int.(object.leadTimes) == collect(0:72) || error("unexpected lead times at $cycle_text")
    field = Float64.(object.staDataAll[1, 1, valid_station_indices, :])
    field[.!isfinite.(field) .| (field .< 0)] .= NaN

    lead_min = max(1, Int(Dates.value(first_valid - cycle) ÷ 3_600_000))
    lead_max = min(72, Int(Dates.value(last_valid - cycle) ÷ 3_600_000))
    columns = (lead_min + 1):(lead_max + 1)
    for (city_index, positions) in enumerate(city_groups)
        city_means = vec(nanmean(@view(field[positions, columns]), dims=1))
        for (offset, lead) in enumerate(lead_min:lead_max)
            time_index = Int(Dates.value((cycle + Hour(lead)) - first_valid) ÷ 3_600_000) + 1
            ietm3d_array[city_index, lead][time_index] = city_means[offset]
        end
    end
    global loaded_cycles += 1
    (cycle_index == 1 || cycle_index % 24 == 0 || cycle_index == length(cycle_times)) &&
        println("  processed cycle $cycle_index/$(length(cycle_times)): $cycle_text")
end
serialize(FIG3_CACHE_FILE, (
    ietm3d_array=ietm3d_array,
    loaded_cycles=loaded_cycles,
    missing_cycles=missing_cycles,
))
println("Saved prepared city-series cache: $FIG3_CACHE_FILE")
end

include(joinpath(@__DIR__, "Dependency", "change_point_detect.jl"))
derivative_threshold = 0.1
concentration_threshold = 200.0
obs_point = Vector{Vector{Int64}}(undef, length(cities2detect))
ietm_point = Matrix{Vector{Int64}}(undef, length(cities2detect), length(leadtimes))
bool_ietm = zeros(Int, length(cities2detect), length(leadtimes))
for city_index in eachindex(cities2detect)
    obs_point[city_index], _ = detect_starting_points(
        obs_array[city_index], derivative_threshold, concentration_threshold
    )
    for lead in leadtimes
        ietm_point[city_index, lead], _ = detect_starting_points(
            ietm3d_array[city_index, lead], derivative_threshold, concentration_threshold
        )
        bool_ietm[city_index, lead] = judge_forcast(
            obs_array[city_index], ietm3d_array[city_index, lead],
            obs_point[city_index], ietm_point[city_index, lead], 6,
        )
    end
end
advance_times = vec(bool2leadtime(bool_ietm))

struct Figure3ForecastCache end
const IETMdata = Figure3ForecastCache()

function readData(::Figure3ForecastCache, Time, leadtime, varname, city)
    city_index = city_lookup[string(city)]
    return copy(ietm3d_array[city_index, Int(leadtime)])
end

println(
    "Figure 3 new-data preparation complete: ", loaded_cycles, "/", length(cycle_times),
    " cycles loaded; ", length(missing_cycles), " missing",
)
