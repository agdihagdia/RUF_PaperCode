#!/usr/bin/env julia

pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "codeLib")))

using CSV
using DataFrames
using Dates
using NaNStatistics
using StaDataM
using Statistics
using TimeZones

include(joinpath(@__DIR__, "..", "change_point_detect.jl"))

const BJT = tz"Asia/Shanghai"
const STAMP = dateformat"yyyymmddHH"
const OBS_FILE = get(ENV, "OBS_FILE", "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat")
const FC_ROOT = get(ENV, "FC_ROOT", "/home/zsm/IETM/RUF_DryDeposition")
const OUTPUT_DIR = get(ENV, "OUTPUT_DIR", joinpath(@__DIR__, "output"))
const OBS_UPPER = parse(Float64, get(ENV, "OBS_UPPER", "2000"))
const DUST_THRESHOLD = parse(Float64, get(ENV, "DUST_THRESHOLD", "200"))
const GRADIENT_THRESHOLD = parse(Float64, get(ENV, "GRADIENT_THRESHOLD", "0.1"))

# These periods reproduce the periods used in the manuscript's original figure scripts.
const MAIN_START = DateTime(2025, 4, 11, 0)
const MAIN_END = DateTime(2025, 4, 21, 0)
const FIG1_START = MAIN_START
const FIG1_END = MAIN_END
const FIG3_START = DateTime(2025, 4, 11, 0)
const FIG3_END = DateTime(2025, 4, 13, 0)
const FIG4_START = DateTime(2025, 4, 11, 0)
const FIG4_END = DateTime(2025, 4, 14, 0)
const FIG5_OBS_START = DateTime(2025, 4, 11, 0)
const FIG5_OBS_END = DateTime(2025, 4, 15, 0)
const FIG5_CYCLE = DateTime(2025, 4, 11, 15)
const APP_START = DateTime(2025, 4, 9, 0)
const APP_END = DateTime(2025, 4, 15, 23)
const CYCLE_START = DateTime(2025, 4, 8, 0)
const CYCLE_END = DateTime(2025, 4, 21, 23)

const FIG3_CITIES = ["延安", "太原", "大同", "忻州"]
const FIG5_CITIES = ["太原", "洛阳", "常德", "海口"]
const APP_CITIES = ["太原", "北京", "南京", "海口"]
const FIG4_REPRESENTATIVE_CITIES = [
    "呼和浩特", "北京", "太原", "洛阳", "襄阳", "宜昌", "常德", "桂林", "茂名", "海口",
]

clean(value) = value === missing ? "" : strip(string(value))
stamp(value::DateTime) = Dates.format(value, STAMP)
hours(first_time, last_time) = collect(first_time:Hour(1):last_time)

function paired_metrics(model, observation)
    valid = isfinite.(model) .& isfinite.(observation)
    count_valid = count(valid)
    count_valid == 0 && return (n=0, rmse=NaN, r=NaN, bias=NaN)
    m = Float64.(model[valid])
    o = Float64.(observation[valid])
    r = count_valid >= 2 && std(m) > 0 && std(o) > 0 ? cor(m, o) : NaN
    return (n=count_valid, rmse=sqrt(mean(abs2, m .- o)), r=r, bias=mean(m .- o))
end

function finite_mean(values)
    selected = Float64.(values[isfinite.(values)])
    return isempty(selected) ? NaN : mean(selected)
end

function finite_max(values)
    selected = Float64.(values[isfinite.(values)])
    return isempty(selected) ? NaN : maximum(selected)
end

function safe_pct(numerator, denominator)
    denominator == 0 && return NaN
    return 100 * numerator / denominator
end

function obs_time_lookup(observations)
    return Dict(
        DateTime(astimezone(time, BJT)) => index
        for (index, time) in enumerate(observations.times)
    )
end

function city_means(data::AbstractMatrix, station_groups)
    output = fill(Float32(NaN), size(data, 1), length(station_groups))
    for (city_index, positions) in enumerate(station_groups)
        values = Float64.(data[:, positions])
        valid = isfinite.(values)
        sums = vec(sum(ifelse.(valid, values, 0.0), dims=2))
        counts = vec(sum(valid, dims=2))
        mask = counts .> 0
        output[mask, city_index] .= Float32.(sums[mask] ./ counts[mask])
    end
    return output
end

function station_type_inventory(info, matched_ids)
    types = clean.(info.type)
    ids = clean.(info.staID)
    matched = ids .∈ Ref(matched_ids)
    rows = NamedTuple[]
    for station_type in sort(unique(types))
        mask = types .== station_type
        push!(rows, (
            station_type=station_type,
            observation_count=count(mask),
            matched_forecast_count=count(mask .& matched),
        ))
    end
    return DataFrame(rows)
end

function fixed_lead_city_series(city_fc, cycle_lookup, city_index, valid_times, lead)
    result = fill(NaN, length(valid_times))
    for (position, valid_time) in enumerate(valid_times)
        cycle = valid_time - Hour(lead)
        cycle_position = get(cycle_lookup, cycle, 0)
        cycle_position == 0 && continue
        result[position] = city_fc[cycle_position, lead + 1, city_index]
    end
    return result
end

function fixed_lead_national_series(national_fc, cycle_lookup, valid_times, lead)
    result = fill(NaN, length(valid_times))
    for (position, valid_time) in enumerate(valid_times)
        cycle_position = get(cycle_lookup, valid_time - Hour(lead), 0)
        cycle_position == 0 && continue
        result[position] = national_fc[cycle_position, lead + 1]
    end
    return result
end

"""Reproduce errortable_zsm3.jl: start at evaluation_start-1 h and take leads 0:F-1."""
function legacy_update_series(city_fc, cycle_lookup, city_index, frequency, valid_times)
    output = fill(NaN, length(valid_times))
    fixed_start = first(valid_times) - Hour(1)
    filled = 0
    block = 0
    while filled < length(valid_times)
        cycle = fixed_start + Hour(block * frequency)
        cycle_position = get(cycle_lookup, cycle, 0)
        take = min(frequency, length(valid_times) - filled)
        if cycle_position != 0
            output[(filled + 1):(filled + take)] .= Float64.(city_fc[cycle_position, 1:take, city_index])
        end
        filled += take
        block += 1
    end
    return output
end

function definitive_arrival(values)
    starts, _ = detect_starting_points(Float64.(values), GRADIENT_THRESHOLD, DUST_THRESHOLD)
    isempty(starts) && return (candidate_count=0, first_candidate=0, definitive=0)
    radius = min(3, length(starts))
    candidates = starts[1:radius]
    safe = [isfinite(values[index]) ? Float64(values[index]) : Inf for index in candidates]
    chosen = candidates[argmin(safe)]
    return (candidate_count=length(starts), first_candidate=first(starts), definitive=chosen)
end

function categorical_counts(obs_values, fc_values, threshold)
    valid = isfinite.(obs_values) .& isfinite.(fc_values) .&
        (obs_values .>= 0) .& (obs_values .< OBS_UPPER) .& (fc_values .>= 0)
    obs_event = valid .& (obs_values .> threshold)
    fc_event = valid .& (fc_values .> threshold)
    hit = count(obs_event .& fc_event)
    false_alarm = count(valid .& .!obs_event .& fc_event)
    miss = count(valid .& obs_event .& .!fc_event)
    correct_negative = count(valid .& .!obs_event .& .!fc_event)
    return (valid=count(valid), hit=hit, false_alarm=false_alarm, miss=miss, correct_negative=correct_negative)
end

function main()
    mkpath(OUTPUT_DIR)
    println("Loading observation: $OBS_FILE")
    observations = staDATA(OBS_FILE)
    pm10_position = findfirst(==("PM10"), uppercase.(observations.varNames))
    pm10_position === nothing && error("PM10 not found in observation file")

    paths = Dict{DateTime,String}()
    for cycle in hours(CYCLE_START, CYCLE_END)
        path = joinpath(FC_ROOT, stamp(cycle), "staData.$(stamp(cycle)).mat")
        isfile(path) && (paths[cycle] = path)
    end
    isempty(paths) && error("no forecast station files found below $FC_ROOT")
    first_cycle = first(sort(collect(keys(paths))))
    template = staDATA(paths[first_cycle])
    forecast_ids = clean.(template.staInfo.staID)
    length(unique(forecast_ids)) == length(forecast_ids) || error("forecast station IDs are not unique")

    observation_ids = clean.(observations.staInfo.staID)
    observation_position = Dict(id => index for (index, id) in enumerate(observation_ids))
    all(haskey(observation_position, id) for id in forecast_ids) || error("forecast IDs are absent from observations")
    obs_positions = [observation_position[id] for id in forecast_ids]
    obs_pm10 = Float64.(observations.staDataAll[:, pm10_position, obs_positions, 1])
    obs_pm10[.!isfinite.(obs_pm10) .| (obs_pm10 .< 0) .| (obs_pm10 .>= OBS_UPPER)] .= NaN
    time_lookup = obs_time_lookup(observations)

    city_labels = clean.(template.staInfo.city)
    city_names = sort(filter(!isempty, unique(city_labels)))
    city_groups = [findall(==(city), city_labels) for city in city_names]
    city_lookup = Dict(city => index for (index, city) in enumerate(city_names))
    obs_city = city_means(obs_pm10, city_groups)
    obs_national = vec(Float64.(nanmean(obs_pm10, dims=2)))

    cycles = hours(CYCLE_START, CYCLE_END)
    cycle_position = Dict(cycle => index for (index, cycle) in enumerate(cycles))
    cycle_lookup = Dict{DateTime,Int}()
    city_fc = fill(Float32(NaN), length(cycles), 73, length(city_names))
    national_fc = fill(Float32(NaN), length(cycles), 73)

    # Appendix S3, national-control type B only.
    b_positions = findall(==("B"), clean.(template.staInfo.type))
    b_count = length(b_positions)
    b_hit = zeros(Int, b_count)
    b_false_alarm = zeros(Int, b_count)
    b_miss = zeros(Int, b_count)
    b_correct_negative = zeros(Int, b_count)
    b_valid = zeros(Int, b_count)

    station_pooled_sse = zeros(Float64, 24)
    station_pooled_n = zeros(Int, 24)
    loaded = 0
    for cycle in cycles
        path = get(paths, cycle, "")
        isempty(path) && continue
        object = staDATA(path)
        clean.(object.staInfo.staID) == forecast_ids || error("station ordering changed in $path")
        field = Float64.(object.staDataAll[1, 1, :, :])
        size(field) == (length(forecast_ids), 73) || error("unexpected forecast dimensions in $path: $(size(field))")
        field[.!isfinite.(field) .| (field .< 0)] .= NaN
        cpos = cycle_position[cycle]
        cycle_lookup[cycle] = cpos
        national_fc[cpos, :] .= Float32.(vec(nanmean(field, dims=1)))
        for (city_index, positions) in enumerate(city_groups)
            city_fc[cpos, :, city_index] .= Float32.(vec(nanmean(@view(field[positions, :]), dims=1)))
        end

        for lead in 1:24
            valid_time = cycle + Hour(lead)
            FIG1_START <= valid_time <= FIG1_END || continue
            obs_index = get(time_lookup, valid_time, 0)
            obs_index == 0 && continue
            fc = @view field[:, lead + 1]
            ob = @view obs_pm10[obs_index, :]
            valid = isfinite.(fc) .& isfinite.(ob)
            station_pooled_sse[lead] += sum(abs2, fc[valid] .- ob[valid])
            station_pooled_n[lead] += count(valid)
        end

        valid_time3 = cycle + Hour(3)
        if APP_START <= valid_time3 <= APP_END
            obs_index = get(time_lookup, valid_time3, 0)
            if obs_index != 0
                fc = @view field[b_positions, 4]
                ob = @view obs_pm10[obs_index, b_positions]
                valid = isfinite.(fc) .& isfinite.(ob)
                obs_event = valid .& (ob .> DUST_THRESHOLD)
                fc_event = valid .& (fc .> DUST_THRESHOLD)
                b_valid .+= valid
                b_hit .+= obs_event .& fc_event
                b_false_alarm .+= valid .& .!obs_event .& fc_event
                b_miss .+= valid .& obs_event .& .!fc_event
                b_correct_negative .+= valid .& .!obs_event .& .!fc_event
            end
        end
        loaded += 1
        (loaded == 1 || loaded % 24 == 0 || loaded == length(paths)) &&
            println("Loaded $loaded/$(length(paths)): $(stamp(cycle))")
    end

    missing_cycles = setdiff(cycles, sort(collect(keys(paths))))
    CSV.write(joinpath(OUTPUT_DIR, "missing_forecast_cycles.csv"), DataFrame(cycle_bjt=stamp.(missing_cycles)))
    CSV.write(
        joinpath(OUTPUT_DIR, "station_inventory.csv"),
        station_type_inventory(observations.staInfo, Set(forecast_ids)),
    )

    # Figure 1(c): original code first averages all stations at each valid hour,
    # then calculates RMSE over the national-mean time series.
    fig1_times = hours(FIG1_START, FIG1_END)
    fig1_obs = [obs_national[time_lookup[time]] for time in fig1_times]
    fig1_rows = NamedTuple[]
    for lead in 1:24
        model = fixed_lead_national_series(national_fc, cycle_lookup, fig1_times, lead)
        metric = paired_metrics(model, fig1_obs)
        station_rmse = station_pooled_n[lead] == 0 ? NaN : sqrt(station_pooled_sse[lead] / station_pooled_n[lead])
        push!(fig1_rows, (
            lead_hour=lead, valid_hours=metric.n, national_mean_timeseries_rmse=metric.rmse,
            national_mean_timeseries_r=metric.r, station_pooled_pairs=station_pooled_n[lead],
            station_pooled_rmse=station_rmse,
        ))
    end
    CSV.write(joinpath(OUTPUT_DIR, "figure1_lead_metrics.csv"), DataFrame(fig1_rows))

    # Figure 3 all-station experiment values. The national-only experiment is not available.
    fig3_times = hours(FIG3_START, FIG3_END)
    fig3_city_rows = NamedTuple[]
    for city in FIG3_CITIES
        ci = city_lookup[city]
        obs_values = [obs_city[time_lookup[time], ci] for time in fig3_times]
        model = fixed_lead_city_series(city_fc, cycle_lookup, ci, fig3_times, 2)
        metric = paired_metrics(model, obs_values)
        push!(fig3_city_rows, (city=city, lead_hour=2, n=metric.n, rmse=metric.rmse, r=metric.r))
    end
    CSV.write(joinpath(OUTPUT_DIR, "figure3_allstations_city_metrics.csv"), DataFrame(fig3_city_rows))

    fig3_snapshot_rows = NamedTuple[]
    for valid_time in [DateTime(2025, 4, 11, hour) for hour in (10, 12, 14, 16)]
        ci_cycle = get(cycle_lookup, valid_time, 0)
        obs_index = time_lookup[valid_time]
        model = ci_cycle == 0 ? fill(NaN, length(city_names)) : Float64.(city_fc[ci_cycle, 1, :])
        ob = Float64.(obs_city[obs_index, :])
        counts = categorical_counts(ob, model, DUST_THRESHOLD)
        ts = counts.hit / max(1, counts.hit + counts.false_alarm + counts.miss)
        push!(fig3_snapshot_rows, (
            valid_time_bjt=stamp(valid_time), actual_lead_hour=0, valid_cities=counts.valid,
            hit=counts.hit, false_alarm=counts.false_alarm, miss=counts.miss, ts=ts,
        ))
    end
    CSV.write(joinpath(OUTPUT_DIR, "figure3_allstations_snapshot_ts.csv"), DataFrame(fig3_snapshot_rows))

    # Figure 6: exact legacy segmentation used by errortable_zsm3.jl.
    main_times = hours(MAIN_START, MAIN_END)
    update_rows = NamedTuple[]
    update_series = Dict{Tuple{Int,Int},Vector{Float64}}()
    for (city_index, city) in enumerate(city_names)
        observation = [obs_city[time_lookup[time], city_index] for time in main_times]
        for frequency in 1:24
            model = legacy_update_series(city_fc, cycle_lookup, city_index, frequency, main_times)
            update_series[(city_index, frequency)] = model
            metric = paired_metrics(model, observation)
            push!(update_rows, (
                city=city, update_frequency_hour=frequency, n=metric.n,
                rmse=metric.rmse, r=metric.r, bias=metric.bias,
            ))
        end
    end
    update_city = DataFrame(update_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure6_city_metrics.csv"), update_city)
    summary_rows = NamedTuple[]
    for frequency in 1:24
        rows = update_city[update_city.update_frequency_hour .== frequency, :]
        push!(summary_rows, (
            update_frequency_hour=frequency,
            evaluated_cities=count(isfinite, rows.rmse),
            mean_city_rmse=finite_mean(rows.rmse),
            mean_city_r=finite_mean(rows.r),
        ))
    end
    update_summary = DataFrame(summary_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure6_frequency_summary.csv"), update_summary)
    CSV.write(
        joinpath(OUTPUT_DIR, "figure6_representative_city_metrics.csv"),
        update_city[(update_city.city .∈ Ref(APP_CITIES)) .& (update_city.update_frequency_hour .∈ Ref([1, 6, 12, 24])), :],
    )

    comparison_rows = NamedTuple[]
    for city in city_names
        u1 = update_city[(update_city.city .== city) .& (update_city.update_frequency_hour .== 1), :]
        u24 = update_city[(update_city.city .== city) .& (update_city.update_frequency_hour .== 24), :]
        (isempty(u1) || isempty(u24)) && continue
        rmse_reduction = isfinite(u1.rmse[1]) && isfinite(u24.rmse[1]) && u24.rmse[1] != 0 ?
            100 * (u24.rmse[1] - u1.rmse[1]) / u24.rmse[1] : NaN
        r_improvement = isfinite(u1.r[1]) && isfinite(u24.r[1]) && abs(u24.r[1]) > 1e-12 ?
            100 * (u1.r[1] - u24.r[1]) / abs(u24.r[1]) : NaN
        push!(comparison_rows, (city=city, rmse_reduction_percent=rmse_reduction, r_improvement_percent=r_improvement))
    end
    update_comparison = DataFrame(comparison_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure4_hourly_vs_daily_city_improvement.csv"), update_comparison)

    # Abstract: compare U1 and U24 only within the first 12 samples of each U24 block.
    first12_city_rows = NamedTuple[]
    first12_mask = [mod(position - 1, 24) < 12 for position in eachindex(main_times)]
    for (city_index, city) in enumerate(city_names)
        observation = [obs_city[time_lookup[time], city_index] for time in main_times]
        metric1 = paired_metrics(update_series[(city_index, 1)][first12_mask], observation[first12_mask])
        metric24 = paired_metrics(update_series[(city_index, 24)][first12_mask], observation[first12_mask])
        reduction = isfinite(metric1.rmse) && isfinite(metric24.rmse) && metric24.rmse != 0 ?
            100 * (metric24.rmse - metric1.rmse) / metric24.rmse : NaN
        push!(first12_city_rows, (
            city=city, u1_n=metric1.n, u1_rmse=metric1.rmse,
            u24_n=metric24.n, u24_rmse=metric24.rmse, rmse_reduction_percent=reduction,
        ))
    end
    first12_city = DataFrame(first12_city_rows)
    CSV.write(joinpath(OUTPUT_DIR, "abstract_initial12_city_metrics.csv"), first12_city)

    # Figure 4 arrival calculation, using the manuscript's 5-h smoother, C=200, G=0.1, R=3, tolerance=6 h.
    fig4_times = hours(FIG4_START, FIG4_END)
    arrival_rows = NamedTuple[]
    for (city_index, city) in enumerate(city_names)
        observation = Float64.([obs_city[time_lookup[time], city_index] for time in fig4_times])
        obs_arrival = definitive_arrival(observation)
        max_warning = 0
        lead24_arrival = 0
        lead24_error = NaN
        lead12_23_first_errors = Float64[]
        for lead in 1:72
            model = fixed_lead_city_series(city_fc, cycle_lookup, city_index, fig4_times, lead)
            fc_arrival = definitive_arrival(model)
            if obs_arrival.definitive > 0 && fc_arrival.definitive > 0 &&
                    abs(obs_arrival.definitive - fc_arrival.definitive) <= 6
                max_warning = lead
            end
            if lead == 24
                lead24_arrival = fc_arrival.definitive
                if obs_arrival.definitive > 0 && fc_arrival.definitive > 0
                    lead24_error = abs(obs_arrival.definitive - fc_arrival.definitive)
                end
            end
            # Original figure3_v2.jl arrival-error block used the first candidate at leads 12:23.
            if 12 <= lead <= 23 && obs_arrival.first_candidate > 0 && fc_arrival.first_candidate > 0
                push!(lead12_23_first_errors, abs(obs_arrival.first_candidate - fc_arrival.first_candidate))
            end
        end
        obs_time = obs_arrival.definitive == 0 ? "" : stamp(fig4_times[obs_arrival.definitive])
        fc24_time = lead24_arrival == 0 ? "" : stamp(fig4_times[lead24_arrival])
        push!(arrival_rows, (
            city=city, observed_event=obs_arrival.definitive > 0,
            observed_arrival_bjt=obs_time, max_warning_lead_hour=max_warning,
            lead24_arrival_bjt=fc24_time, lead24_absolute_arrival_error_hour=lead24_error,
            legacy_lead12_23_first_candidate_mae=finite_mean(lead12_23_first_errors),
        ))
    end
    arrivals = DataFrame(arrival_rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure4_city_arrival_metrics.csv"), arrivals)

    # Figure 5: exact original periods; forecast is all 0:72 leads from one 2025041115 cycle,
    # while observations cover 2025041100:2025041500. Only RUF3D is available in this run.
    fig5_rows = NamedTuple[]
    fig5_cycle_pos = get(cycle_lookup, FIG5_CYCLE, 0)
    fig5_cycle_pos == 0 && error("Figure 5 cycle $(stamp(FIG5_CYCLE)) is missing")
    fig5_obs_times = hours(FIG5_OBS_START, FIG5_OBS_END)
    for city in FIG5_CITIES
        ci = city_lookup[city]
        observation = Float64.([obs_city[time_lookup[time], ci] for time in fig5_obs_times])
        model = Float64.(city_fc[fig5_cycle_pos, :, ci])
        push!(fig5_rows, (
            city=city, ruf3d_cycle_bjt=stamp(FIG5_CYCLE), ruf3d_leads="0-72",
            ruf3d_mean=finite_mean(model), ruf3d_max=finite_max(model),
            observation_period="$(stamp(FIG5_OBS_START))-$(stamp(FIG5_OBS_END))",
            observation_mean=finite_mean(observation), observation_max=finite_max(observation),
        ))
    end
    CSV.write(joinpath(OUTPUT_DIR, "figure5_ruf3d_and_observation_stats.csv"), DataFrame(fig5_rows))

    # Appendix Table S3 and Figure S3 RUF 3-h values.
    total_hit = sum(b_hit)
    total_fa = sum(b_false_alarm)
    total_miss = sum(b_miss)
    total_cn = sum(b_correct_negative)
    far = total_fa / max(1, total_hit + total_fa)
    miss_rate = total_miss / max(1, total_hit + total_miss)
    ts = total_hit / max(1, total_hit + total_fa + total_miss)
    appendix_summary = DataFrame([(
        model="RUF_DryDeposition", lead_hour=3, period="$(stamp(APP_START))-$(stamp(APP_END))",
        station_type="B", station_count=b_count, valid_pairs=sum(b_valid),
        n_correct=total_hit, n_false_alarm=total_fa, n_miss=total_miss,
        n_correct_negative=total_cn, false_alarm_rate=far, miss_rate=miss_rate, ts=ts,
    )])
    CSV.write(joinpath(OUTPUT_DIR, "appendix_table_s3_ruf_summary.csv"), appendix_summary)
    station_appendix = DataFrame(
        station_id=forecast_ids[b_positions], station_name=clean.(template.staInfo.staName[b_positions]),
        city=clean.(template.staInfo.city[b_positions]), province=clean.(template.staInfo.province[b_positions]),
        longitude=Float64.(template.staInfo.lon[b_positions]), latitude=Float64.(template.staInfo.lat[b_positions]),
        valid_pairs=b_valid, n_correct=b_hit, n_false_alarm=b_false_alarm, n_miss=b_miss,
        n_correct_negative=b_correct_negative,
        false_alarm_rate=[(b_hit[i]+b_false_alarm[i]) == 0 ? NaN : b_false_alarm[i]/(b_hit[i]+b_false_alarm[i]) for i in eachindex(b_hit)],
        miss_rate=[(b_hit[i]+b_miss[i]) == 0 ? NaN : b_miss[i]/(b_hit[i]+b_miss[i]) for i in eachindex(b_hit)],
        ts=[(b_hit[i]+b_false_alarm[i]+b_miss[i]) == 0 ? NaN : b_hit[i]/(b_hit[i]+b_false_alarm[i]+b_miss[i]) for i in eachindex(b_hit)],
    )
    CSV.write(joinpath(OUTPUT_DIR, "appendix_table_s3_ruf_station_metrics.csv"), station_appendix)

    app_times = hours(APP_START, APP_END)
    app_city_rows = NamedTuple[]
    for city in APP_CITIES
        ci = city_lookup[city]
        observation = Float64.([obs_city[time_lookup[time], ci] for time in app_times])
        model = fixed_lead_city_series(city_fc, cycle_lookup, ci, app_times, 3)
        metric = paired_metrics(model, observation)
        push!(app_city_rows, (city=city, lead_hour=3, n=metric.n, rmse=metric.rmse, r=metric.r, bias=metric.bias))
    end
    CSV.write(joinpath(OUTPUT_DIR, "appendix_figure3_ruf_city_metrics.csv"), DataFrame(app_city_rows))

    # One compact key-value summary for downstream report generation.
    s6 = Dict(row.update_frequency_hour => row for row in eachrow(update_summary))
    rmse_reduction_12 = 100 * (s6[24].mean_city_rmse - s6[12].mean_city_rmse) / s6[24].mean_city_rmse
    rmse_reduction_6 = 100 * (s6[24].mean_city_rmse - s6[6].mean_city_rmse) / s6[24].mean_city_rmse
    diff12 = update_city[update_city.update_frequency_hour .== 12, :rmse] .-
        update_city[update_city.update_frequency_hour .== 24, :rmse]
    diff6 = update_city[update_city.update_frequency_hour .== 6, :rmse] .-
        update_city[update_city.update_frequency_hour .== 24, :rmse]
    valid_diff12 = isfinite.(diff12)
    valid_diff6 = isfinite.(diff6)
    impacted = arrivals.observed_event
    rmse_valid = isfinite.(update_comparison.rmse_reduction_percent)
    r_valid = isfinite.(update_comparison.r_improvement_percent)
    initial12_reduction = finite_mean(first12_city.rmse_reduction_percent)
    lead24_errors = arrivals.lead24_absolute_arrival_error_hour
    representative_mask = arrivals.city .∈ Ref(FIG4_REPRESENTATIVE_CITIES)
    legacy_errors = arrivals.legacy_lead12_23_first_candidate_mae[representative_mask]
    summary = DataFrame(
        key=String[
            "forecast_files_loaded", "forecast_cycles_missing", "observation_station_count",
            "forecast_station_count", "national_control_B_count", "county_control_C_count",
            "township_control_D_count", "abstract_initial12_mean_city_rmse_reduction_percent",
            "figure6_mean_city_rmse_U24", "figure6_mean_city_rmse_U12", "figure6_mean_city_rmse_U6",
            "figure6_U12_vs_U24_reduction_percent", "figure6_U6_vs_U24_reduction_percent",
            "figure6_city_percent_U12_rmse_drop_over_100", "figure6_city_percent_U6_rmse_drop_over_100",
            "figure4_impacted_city_count", "figure4_impacted_percent_warning_11_20",
            "figure4_impacted_percent_warning_over_20", "figure4_all_city_percent_warning_11_20",
            "figure4_all_city_percent_warning_over_20", "figure4_city_percent_rmse_reduction_over_50",
            "figure4_city_percent_r_improvement_over_30", "figure4_lead24_arrival_mae_impacted_cities",
            "figure4_legacy_representative_lead12_23_first_candidate_mae",
            "appendix_RUF3h_false_alarm_rate", "appendix_RUF3h_miss_rate", "appendix_RUF3h_TS",
        ],
        value=Float64[
            length(paths), length(missing_cycles), nrow(observations.staInfo), length(forecast_ids),
            count(==("B"), clean.(observations.staInfo.type)), count(==("C"), clean.(observations.staInfo.type)),
            count(==("D"), clean.(observations.staInfo.type)), initial12_reduction,
            s6[24].mean_city_rmse, s6[12].mean_city_rmse, s6[6].mean_city_rmse,
            rmse_reduction_12, rmse_reduction_6,
            safe_pct(count(diff12[valid_diff12] .< -100), count(valid_diff12)),
            safe_pct(count(diff6[valid_diff6] .< -100), count(valid_diff6)),
            count(impacted), safe_pct(count((arrivals.max_warning_lead_hour .>= 11) .& (arrivals.max_warning_lead_hour .<= 20) .& impacted), count(impacted)),
            safe_pct(count((arrivals.max_warning_lead_hour .> 20) .& impacted), count(impacted)),
            safe_pct(count((arrivals.max_warning_lead_hour .>= 11) .& (arrivals.max_warning_lead_hour .<= 20)), nrow(arrivals)),
            safe_pct(count(arrivals.max_warning_lead_hour .> 20), nrow(arrivals)),
            safe_pct(count(update_comparison.rmse_reduction_percent[rmse_valid] .> 50), count(rmse_valid)),
            safe_pct(count(update_comparison.r_improvement_percent[r_valid] .> 30), count(r_valid)),
            finite_mean(lead24_errors[impacted]), finite_mean(legacy_errors), far, miss_rate, ts,
        ],
    )
    CSV.write(joinpath(OUTPUT_DIR, "paper_text_summary.csv"), summary)

    println("Completed. Output: $OUTPUT_DIR")
    show(stdout, MIME("text/plain"), summary); println()
end

main()
