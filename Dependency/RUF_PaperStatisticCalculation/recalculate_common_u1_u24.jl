#!/usr/bin/env julia

pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "codeLib")))

using CSV
using DataFrames
using Dates
using NaNStatistics
using StaDataM
using Statistics
using TimeZones

const BJT = tz"Asia/Shanghai"
const STAMP = dateformat"yyyymmddHH"
const OBS_FILE = get(ENV, "OBS_FILE", "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat")
const FC_ROOT = get(ENV, "FC_ROOT", "/home/zsm/IETM/RUF_DryDeposition")
const OUTPUT_DIR = get(ENV, "OUTPUT_DIR", joinpath(@__DIR__, "output"))
const MAIN_START = DateTime(2025, 4, 11, 0)
const MAIN_END = DateTime(2025, 4, 21, 0)

clean(value) = value === missing ? "" : strip(string(value))
stamp(value) = Dates.format(value, STAMP)

function paired_metrics(model, observation)
    valid = isfinite.(model) .& isfinite.(observation)
    n = count(valid)
    n == 0 && return (n=0, rmse=NaN, r=NaN)
    m = Float64.(model[valid])
    o = Float64.(observation[valid])
    r = n >= 2 && std(m) > 0 && std(o) > 0 ? cor(m, o) : NaN
    return (n=n, rmse=sqrt(mean(abs2, m .- o)), r=r)
end

function main()
    main_times = collect(MAIN_START:Hour(1):MAIN_END)
    missing_cycle_labels = string.(
        CSV.read(joinpath(OUTPUT_DIR, "missing_forecast_cycles.csv"), DataFrame).cycle_bjt,
    )
    missing_cycles = DateTime.(missing_cycle_labels, STAMP)
    missing_u1_valid_times = Set(cycle + Hour(1) for cycle in missing_cycles if MAIN_START <= cycle + Hour(1) <= MAIN_END)
    common_time_mask = [!(time in missing_u1_valid_times) for time in main_times]

    observations = staDATA(OBS_FILE)
    pm10_position = findfirst(==("PM10"), uppercase.(observations.varNames))
    first_cycle = MAIN_START - Hour(1)
    first_path = joinpath(FC_ROOT, stamp(first_cycle), "staData.$(stamp(first_cycle)).mat")
    template = staDATA(first_path)
    forecast_ids = clean.(template.staInfo.staID)
    observation_position = Dict(clean(id) => index for (index, id) in enumerate(observations.staInfo.staID))
    obs_positions = [observation_position[id] for id in forecast_ids]
    obs_lookup = Dict(DateTime(astimezone(time, BJT)) => index for (index, time) in enumerate(observations.times))
    obs_pm10 = Float64.(observations.staDataAll[:, pm10_position, obs_positions, 1])
    obs_pm10[.!isfinite.(obs_pm10) .| (obs_pm10 .< 0) .| (obs_pm10 .>= 2000)] .= NaN

    city_labels = clean.(template.staInfo.city)
    city_names = sort(filter(!isempty, unique(city_labels)))
    city_groups = [findall(==(city), city_labels) for city in city_names]
    obs_city = fill(NaN, length(main_times), length(city_names))
    for (city_index, positions) in enumerate(city_groups)
        values = obs_pm10[[obs_lookup[time] for time in main_times], positions]
        obs_city[:, city_index] .= vec(nanmean(values, dims=2))
    end

    # Exact historical U24 segmentation: cycles at 2025041023 + k*24 h,
    # using leads 0:23 and assigning them to the following 24 verification hours.
    u24 = fill(NaN, length(main_times), length(city_names))
    filled = 0
    block = 0
    while filled < length(main_times)
        cycle = first_cycle + Hour(block * 24)
        path = joinpath(FC_ROOT, stamp(cycle), "staData.$(stamp(cycle)).mat")
        isfile(path) || error("U24 source cycle missing: $path")
        object = staDATA(path)
        field = Float64.(object.staDataAll[1, 1, :, :])
        take = min(24, length(main_times) - filled)
        for (city_index, positions) in enumerate(city_groups)
            u24[(filled + 1):(filled + take), city_index] .=
                vec(nanmean(@view(field[positions, 1:take]), dims=1))
        end
        filled += take
        block += 1
    end

    original = CSV.read(joinpath(OUTPUT_DIR, "figure6_city_metrics.csv"), DataFrame)
    u1_rows = original[original.update_frequency_hour .== 1, :]
    u1_lookup = Dict(row.city => row for row in eachrow(u1_rows))
    rows = NamedTuple[]
    for (city_index, city) in enumerate(city_names)
        common_u24 = paired_metrics(u24[common_time_mask, city_index], obs_city[common_time_mask, city_index])
        u1 = u1_lookup[city]
        same_count = Int(u1.n) == common_u24.n
        rmse_reduction = 100 * (common_u24.rmse - Float64(u1.rmse)) / common_u24.rmse
        r_improvement = abs(common_u24.r) > 1e-12 ? 100 * (Float64(u1.r) - common_u24.r) / abs(common_u24.r) : NaN
        push!(rows, (
            city=city, common_hours=common_u24.n, u1_rmse=Float64(u1.rmse),
            u24_common_rmse=common_u24.rmse, rmse_reduction_percent=rmse_reduction,
            u1_r=Float64(u1.r), u24_common_r=common_u24.r,
            r_improvement_percent=r_improvement, count_check_passed=same_count,
        ))
    end
    result = DataFrame(rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure4_hourly_vs_daily_common_hours.csv"), result)
    rmse_valid = isfinite.(result.rmse_reduction_percent)
    r_valid = isfinite.(result.r_improvement_percent)
    summary = DataFrame([(
        evaluated_cities=nrow(result), common_hours_per_complete_city=count(common_time_mask),
        all_count_checks_passed=all(result.count_check_passed),
        cities_rmse_reduction_over_50=count(result.rmse_reduction_percent[rmse_valid] .> 50),
        city_percent_rmse_reduction_over_50=100 * count(result.rmse_reduction_percent[rmse_valid] .> 50) / count(rmse_valid),
        cities_r_improvement_over_30=count(result.r_improvement_percent[r_valid] .> 30),
        city_percent_r_improvement_over_30=100 * count(result.r_improvement_percent[r_valid] .> 30) / count(r_valid),
    )])
    CSV.write(joinpath(OUTPUT_DIR, "figure4_hourly_vs_daily_common_hours_summary.csv"), summary)
    show(stdout, MIME("text/plain"), summary); println()
end

main()
