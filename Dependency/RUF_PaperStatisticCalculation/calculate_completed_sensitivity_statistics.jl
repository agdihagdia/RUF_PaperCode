#!/usr/bin/env julia

pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "codeLib")))

using CSV
using DataFrames
using Dates
using StaDataM
using Statistics
using TimeZones

include(joinpath(@__DIR__, "..", "change_point_detect.jl"))

const BJT = tz"Asia/Shanghai"
const STAMP = dateformat"yyyymmddHH"
const OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const ALL_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const NATION_ROOT = "/home/zsm/IETM/RUF_DryDeposition_Nation/Ruf_forcast"
const RUF2D_FILE = "/home/zsm/IETM/RUF_2d/DATA/RUF/2025041118/staData.2025041118.mat"
const RUF3D_FILE = "/home/zsm/IETM/RUF_DryDeposition/2025041118/staData.2025041118.mat"
const OUTPUT_DIR = get(ENV, "OUTPUT_DIR", joinpath(@__DIR__, "output"))
const DUST_THRESHOLD = 200.0
const GRADIENT_THRESHOLD = 0.1
const OBS_UPPER = 2000.0
const FIG3_CITIES = ["延安", "太原", "大同", "忻州"]
const FIG5_CITIES = ["太原", "洛阳", "常德", "海口"]

clean(value) = value === missing ? "" : strip(string(value))
stamp(value::DateTime) = Dates.format(value, STAMP)
hours(first_time::DateTime, last_time::DateTime) = collect(first_time:Hour(1):last_time)

function finite_values(values)
    return Float64.(values[isfinite.(values)])
end

function finite_mean(values)
    selected = finite_values(values)
    return isempty(selected) ? NaN : mean(selected)
end

function finite_max(values)
    selected = finite_values(values)
    return isempty(selected) ? NaN : maximum(selected)
end

function paired_metrics(model, observation)
    valid = isfinite.(model) .& isfinite.(observation)
    n = count(valid)
    n == 0 && return (n=0, rmse=NaN, r=NaN, bias=NaN)
    m = Float64.(model[valid])
    o = Float64.(observation[valid])
    r = n >= 2 && std(m) > 0 && std(o) > 0 ? cor(m, o) : NaN
    return (n=n, rmse=sqrt(mean(abs2, m .- o)), r=r, bias=mean(m .- o))
end

function categorical_counts(observation, forecast)
    valid = isfinite.(observation) .& isfinite.(forecast)
    observed_event = valid .& (observation .> DUST_THRESHOLD)
    forecast_event = valid .& (forecast .> DUST_THRESHOLD)
    hit = count(observed_event .& forecast_event)
    false_alarm = count(valid .& .!observed_event .& forecast_event)
    miss = count(valid .& observed_event .& .!forecast_event)
    denominator = hit + false_alarm + miss
    ts = denominator == 0 ? NaN : hit / denominator
    return (valid=count(valid), hit=hit, false_alarm=false_alarm, miss=miss, ts=ts)
end

function city_means_time_station(data::AbstractMatrix, groups)
    result = fill(NaN, size(data, 1), length(groups))
    for (city_index, positions) in enumerate(groups)
        for time_index in axes(data, 1)
            result[time_index, city_index] = finite_mean(view(data, time_index, positions))
        end
    end
    return result
end

function city_means_station_lead(data::AbstractMatrix, groups)
    result = fill(NaN, size(data, 2), length(groups))
    for (city_index, positions) in enumerate(groups)
        for lead_index in axes(data, 2)
            result[lead_index, city_index] = finite_mean(view(data, positions, lead_index))
        end
    end
    return result
end

function forecast_matrix(object)
    size(object.staDataAll, 1) == 1 || error("unexpected forecast time dimension")
    size(object.staDataAll, 2) == 1 || error("unexpected forecast variable dimension")
    return Float64.(object.staDataAll[1, 1, :, :])
end

function definitive_arrival(values)
    starts, _ = detect_starting_points(Float64.(values), GRADIENT_THRESHOLD, DUST_THRESHOLD)
    isempty(starts) && return 0
    candidates = starts[1:min(3, length(starts))]
    safe_values = [isfinite(values[index]) ? Float64(values[index]) : Inf for index in candidates]
    return candidates[argmin(safe_values)]
end

function main()
    mkpath(OUTPUT_DIR)
    println("loading observations: ", OBS_FILE)
    observations = staDATA(OBS_FILE)
    pm10_position = findfirst(==("PM10"), uppercase.(clean.(observations.varNames)))
    pm10_position === nothing && error("PM10 is absent from observations")

    all_template = staDATA(RUF3D_FILE)
    forecast_ids = clean.(all_template.staInfo.staID)
    length(unique(forecast_ids)) == length(forecast_ids) || error("forecast staID is not unique")
    observation_position = Dict(clean(id) => index for (index, id) in enumerate(observations.staInfo.staID))
    all(haskey(observation_position, id) for id in forecast_ids) || error("forecast stations are absent from observations")
    matched_observation_positions = [observation_position[id] for id in forecast_ids]
    obs_pm10 = Float64.(observations.staDataAll[:, pm10_position, matched_observation_positions, 1])
    obs_pm10[.!isfinite.(obs_pm10) .| (obs_pm10 .< 0) .| (obs_pm10 .>= OBS_UPPER)] .= NaN
    observation_time_lookup = Dict(
        DateTime(astimezone(time, BJT)) => index
        for (index, time) in enumerate(observations.times)
    )

    station_cities = clean.(all_template.staInfo.city)
    station_provinces = clean.(all_template.staInfo.province)
    all_cities = sort(filter(!isempty, unique(station_cities)))
    all_city_groups = [findall(==(city), station_cities) for city in all_cities]
    all_city_lookup = Dict(city => index for (index, city) in enumerate(all_cities))
    obs_city = city_means_time_station(obs_pm10, all_city_groups)

    bth_mask = station_provinces .∈ Ref(Set(["北京", "天津", "河北"]))
    bth_cities = sort(filter(!isempty, unique(station_cities[bth_mask])))
    selected_cities = unique(vcat(bth_cities, FIG3_CITIES))
    selected_groups = [findall(==(city), station_cities) for city in selected_cities]
    selected_lookup = Dict(city => index for (index, city) in enumerate(selected_cities))

    # The shared cycle cache enforces identical valid samples for both experiments.
    cycle_first = DateTime(2025, 4, 8, 0)
    cycle_last = DateTime(2025, 4, 13, 23)
    common_cycles = DateTime[]
    nation_cache = Dict{DateTime,Matrix{Float64}}()
    all_cache = Dict{DateTime,Matrix{Float64}}()
    missing_rows = NamedTuple[]
    requested_cycles = hours(cycle_first, cycle_last)
    for (cycle_index, cycle) in enumerate(requested_cycles)
        cycle_stamp = stamp(cycle)
        nation_path = joinpath(NATION_ROOT, cycle_stamp, "staData.$cycle_stamp.mat")
        all_path = joinpath(ALL_ROOT, cycle_stamp, "staData.$cycle_stamp.mat")
        nation_exists = isfile(nation_path)
        all_exists = isfile(all_path)
        if !(nation_exists && all_exists)
            push!(missing_rows, (
                cycle_bjt=cycle_stamp,
                nation_available=nation_exists,
                allstations_available=all_exists,
            ))
            continue
        end
        nation_object = staDATA(nation_path)
        all_object = staDATA(all_path)
        clean.(nation_object.staInfo.staID) == forecast_ids || error("Nation station order changed at $cycle_stamp")
        clean.(all_object.staInfo.staID) == forecast_ids || error("All-stations station order changed at $cycle_stamp")
        nation_cache[cycle] = city_means_station_lead(forecast_matrix(nation_object), selected_groups)
        all_cache[cycle] = city_means_station_lead(forecast_matrix(all_object), selected_groups)
        push!(common_cycles, cycle)
        cycle_index % 12 == 0 && println("loaded shared network cycles ", cycle_index, "/", length(requested_cycles))
    end
    CSV.write(joinpath(OUTPUT_DIR, "completed_network_missing_cycles.csv"), DataFrame(missing_rows))
    isempty(common_cycles) && error("no common Nation/all-stations cycles")
    println("common network cycles=", length(common_cycles), " requested=", length(requested_cycles))

    # Figure 3 city time-series statistics, fixed true lead 2 h, common samples only.
    fig3_times = hours(DateTime(2025, 4, 11, 0), DateTime(2025, 4, 13, 0))
    fig3_rows = NamedTuple[]
    comparison_rows = NamedTuple[]
    for city in FIG3_CITIES
        selected_index = selected_lookup[city]
        all_city_index = all_city_lookup[city]
        observation_values = [obs_city[observation_time_lookup[time], all_city_index] for time in fig3_times]
        nation_values = fill(NaN, length(fig3_times))
        all_values = fill(NaN, length(fig3_times))
        for (position, valid_time) in enumerate(fig3_times)
            cycle = valid_time - Hour(2)
            haskey(nation_cache, cycle) || continue
            nation_values[position] = nation_cache[cycle][3, selected_index]
            all_values[position] = all_cache[cycle][3, selected_index]
        end
        common_valid = isfinite.(observation_values) .& isfinite.(nation_values) .& isfinite.(all_values)
        nation_metric = paired_metrics(nation_values[common_valid], observation_values[common_valid])
        all_metric = paired_metrics(all_values[common_valid], observation_values[common_valid])
        push!(fig3_rows, (experiment="Nation-only", city=city, lead_hour=2, n=nation_metric.n, rmse=nation_metric.rmse, r=nation_metric.r))
        push!(fig3_rows, (experiment="All-stations", city=city, lead_hour=2, n=all_metric.n, rmse=all_metric.rmse, r=all_metric.r))
        rmse_reduction = 100 * (nation_metric.rmse - all_metric.rmse) / nation_metric.rmse
        r_improvement = 100 * (all_metric.r - nation_metric.r) / nation_metric.r
        push!(comparison_rows, (
            city=city,
            common_n=nation_metric.n,
            nation_rmse=nation_metric.rmse,
            allstations_rmse=all_metric.rmse,
            rmse_reduction_percent=rmse_reduction,
            nation_r=nation_metric.r,
            allstations_r=all_metric.r,
            r_relative_improvement_percent=r_improvement,
        ))
    end
    fig3_metrics = DataFrame(fig3_rows)
    fig3_comparison = DataFrame(comparison_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure3_network_city_metrics.csv"), fig3_metrics)
    CSV.write(joinpath(OUTPUT_DIR, "figure3_network_city_improvement.csv"), fig3_comparison)

    # Figure 3 snapshot TS values, fixed true lead 2 h and common city masks.
    snapshot_rows = NamedTuple[]
    for valid_time in [DateTime(2025, 4, 11, hour) for hour in (10, 12, 14, 16)]
        cycle = valid_time - Hour(2)
        cycle_stamp = stamp(cycle)
        nation_path = joinpath(NATION_ROOT, cycle_stamp, "staData.$cycle_stamp.mat")
        all_path = joinpath(ALL_ROOT, cycle_stamp, "staData.$cycle_stamp.mat")
        isfile(nation_path) && isfile(all_path) || error("snapshot common cycle missing: $cycle_stamp")
        nation_object = staDATA(nation_path)
        all_object = staDATA(all_path)
        nation_city = city_means_station_lead(forecast_matrix(nation_object), all_city_groups)[3, :]
        all_city_forecast = city_means_station_lead(forecast_matrix(all_object), all_city_groups)[3, :]
        observation_city = obs_city[observation_time_lookup[valid_time], :]
        common_valid = isfinite.(observation_city) .& isfinite.(nation_city) .& isfinite.(all_city_forecast)
        nation_count = categorical_counts(observation_city[common_valid], nation_city[common_valid])
        all_count = categorical_counts(observation_city[common_valid], all_city_forecast[common_valid])
        improvement = 100 * (all_count.ts - nation_count.ts) / nation_count.ts
        for (experiment, counts) in (("Nation-only", nation_count), ("All-stations", all_count))
            push!(snapshot_rows, (
                valid_time_bjt=stamp(valid_time),
                lead_hour=2,
                experiment=experiment,
                common_valid_cities=count(common_valid),
                hit=counts.hit,
                false_alarm=counts.false_alarm,
                miss=counts.miss,
                ts=counts.ts,
                allstations_relative_ts_improvement_percent=experiment == "All-stations" ? improvement : NaN,
            ))
        end
    end
    snapshot_metrics = DataFrame(snapshot_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure3_network_snapshot_ts_lead2.csv"), snapshot_metrics)

    # BTH warning comparison using the same arrival algorithm and common cycle cache.
    fig4_times = hours(DateTime(2025, 4, 11, 0), DateTime(2025, 4, 14, 0))
    bth_rows = NamedTuple[]
    for city in bth_cities
        selected_index = selected_lookup[city]
        all_city_index = all_city_lookup[city]
        observation_values = [obs_city[observation_time_lookup[time], all_city_index] for time in fig4_times]
        observed_arrival_index = definitive_arrival(observation_values)
        experiment_max = Dict("Nation-only" => 0, "All-stations" => 0)
        for (experiment, cache) in (("Nation-only", nation_cache), ("All-stations", all_cache))
            max_warning = 0
            forecast_arrival_at_max = 0
            for lead in 1:72
                model = fill(NaN, length(fig4_times))
                for (position, valid_time) in enumerate(fig4_times)
                    cycle = valid_time - Hour(lead)
                    haskey(cache, cycle) || continue
                    model[position] = cache[cycle][lead + 1, selected_index]
                end
                forecast_arrival_index = definitive_arrival(model)
                if observed_arrival_index > 0 && forecast_arrival_index > 0 &&
                        abs(observed_arrival_index - forecast_arrival_index) <= 6
                    max_warning = lead
                    forecast_arrival_at_max = forecast_arrival_index
                end
            end
            experiment_max[experiment] = max_warning
            push!(bth_rows, (
                city=city,
                experiment=experiment,
                observed_event=observed_arrival_index > 0,
                observed_arrival_bjt=observed_arrival_index == 0 ? "" : stamp(fig4_times[observed_arrival_index]),
                max_warning_lead_hour=max_warning,
                forecast_arrival_at_max_bjt=forecast_arrival_at_max == 0 ? "" : stamp(fig4_times[forecast_arrival_at_max]),
            ))
        end
        push!(bth_rows, (
            city=city,
            experiment="All-minus-Nation",
            observed_event=observed_arrival_index > 0,
            observed_arrival_bjt=observed_arrival_index == 0 ? "" : stamp(fig4_times[observed_arrival_index]),
            max_warning_lead_hour=experiment_max["All-stations"] - experiment_max["Nation-only"],
            forecast_arrival_at_max_bjt="",
        ))
    end
    bth_metrics = DataFrame(bth_rows)
    CSV.write(joinpath(OUTPUT_DIR, "bth_network_warning_comparison.csv"), bth_metrics)

    # Figure 5 unified single-cycle RUF2D/RUF3D comparison at 2025041118, leads 0:72.
    ruf2d = staDATA(RUF2D_FILE)
    ruf3d = staDATA(RUF3D_FILE)
    ids2d = clean.(ruf2d.staInfo.staID)
    ids3d = clean.(ruf3d.staInfo.staID)
    ids2d == ids3d == forecast_ids || error("RUF2D/RUF3D station order differs")
    collect(Int.(ruf2d.leadTimes)) == collect(0:72) || error("RUF2D leads are not 0:72")
    collect(Int.(ruf3d.leadTimes)) == collect(0:72) || error("RUF3D leads are not 0:72")
    data2d = forecast_matrix(ruf2d)
    data3d = forecast_matrix(ruf3d)
    fig5_start = DateTime(2025, 4, 11, 18)
    fig5_times = [fig5_start + Hour(lead) for lead in 0:72]
    fig5_rows = NamedTuple[]
    for city in FIG5_CITIES
        group = all_city_groups[all_city_lookup[city]]
        series2d = vec(city_means_station_lead(data2d, [group])[:, 1])
        series3d = vec(city_means_station_lead(data3d, [group])[:, 1])
        observation = [obs_city[observation_time_lookup[time], all_city_lookup[city]] for time in fig5_times]
        mean2d = finite_mean(series2d)
        max2d = finite_max(series2d)
        mean3d = finite_mean(series3d)
        max3d = finite_max(series3d)
        obs_mean = finite_mean(observation)
        obs_max = finite_max(observation)
        mean_error2d = abs(mean2d - obs_mean)
        mean_error3d = abs(mean3d - obs_mean)
        peak_error2d = abs(max2d - obs_max)
        peak_error3d = abs(max3d - obs_max)
        push!(fig5_rows, (
            city=city,
            cycle_bjt=stamp(fig5_start),
            leads="0-72",
            valid_period_bjt="$(stamp(first(fig5_times)))-$(stamp(last(fig5_times)))",
            ruf2d_mean=mean2d,
            ruf2d_max=max2d,
            ruf3d_mean=mean3d,
            ruf3d_max=max3d,
            observation_mean=obs_mean,
            observation_max=obs_max,
            ruf2d_mean_absolute_error=mean_error2d,
            ruf3d_mean_absolute_error=mean_error3d,
            mean_error_reduction_percent=100 * (mean_error2d - mean_error3d) / mean_error2d,
            ruf2d_peak_absolute_error=peak_error2d,
            ruf3d_peak_absolute_error=peak_error3d,
            peak_error_reduction_percent=100 * (peak_error2d - peak_error3d) / peak_error2d,
        ))
    end
    fig5_metrics = DataFrame(fig5_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure5_ruf2d_ruf3d_observation_stats_2025041118.csv"), fig5_metrics)

    # Compact summary used to update manuscript text.
    snapshot_improvements = snapshot_metrics[
        snapshot_metrics.experiment .== "All-stations",
        :allstations_relative_ts_improvement_percent,
    ]
    bth_differences = bth_metrics[
        (bth_metrics.experiment .== "All-minus-Nation") .& bth_metrics.observed_event,
        :max_warning_lead_hour,
    ]
    downstream = fig5_metrics[fig5_metrics.city .∈ Ref(Set(["常德", "海口"])), :]
    summary = DataFrame(
        metric = [
            "network_common_cycles",
            "figure3_mean_city_rmse_reduction_percent",
            "figure3_mean_city_r_relative_improvement_percent",
            "figure3_snapshot_ts_relative_improvement_min_percent",
            "figure3_snapshot_ts_relative_improvement_max_percent",
            "bth_observed_city_count",
            "bth_warning_difference_min_hour",
            "bth_warning_difference_max_hour",
            "beijing_warning_difference_hour",
            "tianjin_warning_difference_hour",
            "figure5_downstream_mean_error_reduction_min_percent",
            "figure5_downstream_mean_error_reduction_max_percent",
            "figure5_peak_error_reduction_max_percent",
            "figure5_ruf2d_downstream_peak_max",
        ],
        value = [
            length(common_cycles),
            finite_mean(fig3_comparison.rmse_reduction_percent),
            finite_mean(fig3_comparison.r_relative_improvement_percent),
            minimum(finite_values(snapshot_improvements)),
            maximum(finite_values(snapshot_improvements)),
            length(bth_differences),
            minimum(finite_values(bth_differences)),
            maximum(finite_values(bth_differences)),
            only(bth_metrics[(bth_metrics.city .== "北京") .& (bth_metrics.experiment .== "All-minus-Nation"), :max_warning_lead_hour]),
            only(bth_metrics[(bth_metrics.city .== "天津") .& (bth_metrics.experiment .== "All-minus-Nation"), :max_warning_lead_hour]),
            minimum(finite_values(downstream.mean_error_reduction_percent)),
            maximum(finite_values(downstream.mean_error_reduction_percent)),
            maximum(finite_values(fig5_metrics.peak_error_reduction_percent)),
            maximum(finite_values(downstream.ruf2d_max)),
        ],
    )
    CSV.write(joinpath(OUTPUT_DIR, "completed_sensitivity_summary.csv"), summary)

    println("wrote completed sensitivity outputs to ", OUTPUT_DIR)
    show(stdout, MIME("text/plain"), summary)
    println()
end

main()
