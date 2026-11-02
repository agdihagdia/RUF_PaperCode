#!/usr/bin/env julia

push!(LOAD_PATH, "/home/zsm/codeLib")

using CSV
using DataFrames
using Dates
using GradDataM
using JLD2
using PyCall
using StaDataM
using Statistics
using TimeZones
using utilities

include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle

const BJT = tz"Asia/Shanghai"
const EVAL_FIRST = ZonedDateTime(DateTime(2025, 4, 10, 0), BJT)
const EVAL_LAST = ZonedDateTime(DateTime(2025, 4, 16, 0), BJT)
const EVAL_TIMES = collect(EVAL_FIRST:Hour(1):EVAL_LAST)
const CITIES = ["呼和浩特", "太原", "北京", "洛阳", "南京", "海口"]
const CITY_ENGLISH = ["Hohhot", "Taiyuan", "Beijing", "Luoyang", "Nanjing", "Haikou"]
const SERIES_CITY_INDICES = [2, 3, 5, 6]

const OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const RUF_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const RUF_GRID = "/home/zsm/IETM/RUF_2d/datagrid.jld2"
const AURORA_ROOT = "/home/zsm/IETM/Aurora"
const AI_ROOT = "/home/zsm/IETM/AI-GAMFS/inference"
const LEAD_STATISTICS = "/home/zsm/IETM/figure/PapreFigure/20260818/PaperFigure_AppendixCompare3Modle_statistics.csv"
const OUTPUT_DIR = "/home/zsm/IETM/figure/PapreFigure/20260826"
const OUTPUT_BASE = "PaperFigure_AppendixCompare3ModelSpatial_PeakWindow"
const STAMP_FORMAT = dateformat"yyyymmddHH"
const PEAK_HALF_WINDOW_HOURS = 24
const DISPLAY_LEAD_HOURS = 12

# The spatial panels and displayed city time series compare Aurora and RUF at
# the same 12 h forecast lead. AI-GAMFS is retained only in legacy helper
# definitions below and is not read, plotted, or evaluated.
# Temporal interpolation is used only by the lead-time verification curves.
const RUF_MAP_TIMES = [
    DateTime(2025, 4, 11, 8), DateTime(2025, 4, 12, 8), DateTime(2025, 4, 13, 8),
]
const AURORA_MAP_SPECS = [
    ("aurora_20250410T12.nc", 12, DateTime(2025, 4, 11, 8)),
    ("aurora_20250411T12.nc", 12, DateTime(2025, 4, 12, 8)),
    ("aurora_20250412T12.nc", 12, DateTime(2025, 4, 13, 8)),
]
const AI_MAP_SPECS = [
    ("20250410_2230", "AI_GAMFS.20250410_2230+20250411_0130.V01.nc", DateTime(2025, 4, 11, 9)),
    ("20250411_2230", "AI_GAMFS.20250411_2230+20250412_0130.V01.nc", DateTime(2025, 4, 12, 9)),
    ("20250412_2230", "AI_GAMFS.20250412_2230+20250413_0130.V01.nc", DateTime(2025, 4, 13, 9)),
]

# MATLAB character arrays can arrive as vectors such as ["太", "原"]. Keep
# vector joining here; string(value) would make city matching fail.
clean_text(value) = value === missing ? "" : value isa AbstractVector ? join(string.(value)) : strip(string(value))
valid_obs(value) = isfinite(value) && 0.0 <= value < 2000.0
valid_forecast(value) = isfinite(value) && value >= 0.0
time_key(value) = DateTime(astimezone(value, BJT))
time_iso(value) = Dates.format(time_key(value), dateformat"yyyy-mm-ddTHH:MM:SS")

function finite_mean(values, validator)
    selected = Float64[value for value in values if validator(value)]
    return isempty(selected) ? NaN : mean(selected)
end

mutable struct PeakAccumulator
    n::Vector{Int}
    sum_error2::Vector{Float64}
    sum_obs::Vector{Float64}
    sum_fcst::Vector{Float64}
    sum_obs2::Vector{Float64}
    sum_fcst2::Vector{Float64}
    sum_cross::Vector{Float64}
end

PeakAccumulator(size::Int) = PeakAccumulator(
    zeros(Int, size), zeros(size), zeros(size), zeros(size),
    zeros(size), zeros(size), zeros(size),
)

function update_peak!(accumulator::PeakAccumulator, reference, predict, in_peak_window)
    valid = in_peak_window .& valid_obs.(reference) .& valid_forecast.(predict)
    indices = findall(valid)
    isempty(indices) && return
    obs = Float64.(reference[indices])
    fcst = Float64.(predict[indices])
    error = fcst .- obs
    accumulator.n[indices] .+= 1
    accumulator.sum_error2[indices] .+= error .^ 2
    accumulator.sum_obs[indices] .+= obs
    accumulator.sum_fcst[indices] .+= fcst
    accumulator.sum_obs2[indices] .+= obs .^ 2
    accumulator.sum_fcst2[indices] .+= fcst .^ 2
    accumulator.sum_cross[indices] .+= obs .* fcst
end

function finalize_peak(accumulator::PeakAccumulator)
    sample_count = Float64.(accumulator.n)
    station_rmse = fill(NaN, length(sample_count))
    valid_rmse = accumulator.n .> 0
    station_rmse[valid_rmse] .= sqrt.(accumulator.sum_error2[valid_rmse] ./ sample_count[valid_rmse])

    numerator = sample_count .* accumulator.sum_cross .- accumulator.sum_obs .* accumulator.sum_fcst
    denominator2 = (sample_count .* accumulator.sum_obs2 .- accumulator.sum_obs .^ 2) .*
        (sample_count .* accumulator.sum_fcst2 .- accumulator.sum_fcst .^ 2)
    station_r = fill(NaN, length(sample_count))
    valid_r = (accumulator.n .>= 3) .& (denominator2 .> 0.0)
    station_r[valid_r] .= numerator[valid_r] ./ sqrt.(denominator2[valid_r])

    finite_rmse = isfinite.(station_rmse)
    finite_r = isfinite.(station_r)
    return (
        rmse=any(finite_rmse) ? mean(station_rmse[finite_rmse]) : NaN,
        corr=any(finite_r) ? mean(station_r[finite_r]) : NaN,
        rmse_station_count=count(finite_rmse),
        r_station_count=count(finite_r),
        paired_values=sum(accumulator.n),
    )
end

function pooled_lead_rmse(table::DataFrame, model::String)
    selected = table[table.model .== model, :]
    leads = sort(unique(Int.(selected.lead_hour)))
    rmse = Float64[]
    counts = Int[]
    for lead in leads
        rows = selected[selected.lead_hour .== lead, :]
        nrow(rows) == length(CITIES) || error("Expected $(length(CITIES)) city rows for $model lead $lead, found $(nrow(rows))")
        valid = isfinite.(Float64.(rows.rmse_ug_m3)) .& (Int.(rows.valid_city_pairs) .> 0)
        n = Int.(rows.valid_city_pairs[valid])
        e = Float64.(rows.rmse_ug_m3[valid])
        isempty(n) && error("No valid RMSE samples for $model lead $lead")
        # RMSE_i^2 * n_i is the city sum of squared errors. Summing it and
        # dividing by the total count gives the pooled six-city RMSE.
        push!(rmse, sqrt(sum(n .* e .^ 2) / sum(n)))
        push!(counts, sum(n))
    end
    return leads, rmse, counts
end

function pooled_lead_correlation(table::DataFrame, model::String)
    selected = table[table.model .== model, :]
    leads = sort(unique(Int.(selected.lead_hour)))
    correlations = Float64[]
    counts = Int[]
    for lead in leads
        rows = selected[selected.lead_hour .== lead, :]
        nrow(rows) == length(CITIES) || error("Expected $(length(CITIES)) city rows for $model lead $lead, found $(nrow(rows))")
        sample_counts = Int.(rows.valid_city_pairs)
        city_r = Float64.(rows.correlation_r)
        valid = isfinite.(city_r) .& (sample_counts .> 3)
        any(valid) || error("No valid correlations for $model lead $lead")
        weights = sample_counts[valid] .- 3
        bounded_r = clamp.(city_r[valid], -0.999999, 0.999999)
        push!(correlations, tanh(sum(weights .* atanh.(bounded_r)) / sum(weights)))
        push!(counts, sum(sample_counts[valid]))
    end
    return leads, correlations, counts
end

function main()
    mkpath(OUTPUT_DIR)
    println("Evaluation window (BJT): $(time_iso(EVAL_FIRST)) -- $(time_iso(EVAL_LAST))")

    observations = staDATA(OBS_FILE)
    info = observations.staInfo
    cities_all = clean_text.(info.city)
    types_all = clean_text.(info.type)
    ids_all = clean_text.(info.staID)
    station_mask = (types_all .== "B") .& in.(cities_all, Ref(CITIES))
    selected_positions = findall(station_mask)
    isempty(selected_positions) && error("No national-control stations were selected")

    station_ids = ids_all[selected_positions]
    station_cities = cities_all[selected_positions]
    station_lons = Float64.(info.lon[selected_positions])
    station_lats = Float64.(info.lat[selected_positions])
    city_positions = [findall(==(city), station_cities) for city in CITIES]
    any(isempty, city_positions) && error("At least one requested city has no national-control station")
    println("National-control station counts: ", join(["$(CITIES[i])=$(length(city_positions[i]))" for i in eachindex(CITIES)], ", "))

    pm10_position = findfirst(==("PM10"), uppercase.(clean_text.(observations.varNames)))
    pm10_position === nothing && error("PM10 is absent from $OBS_FILE")
    obs_time_lookup = Dict(time_key(t) => i for (i, t) in enumerate(observations.times))
    obs_indices = [get(obs_time_lookup, time_key(t), 0) for t in EVAL_TIMES]
    any(==(0), obs_indices) && error("Observation hours are missing inside the evaluation window")
    obs_field = Float64.(observations.staDataAll[obs_indices, pm10_position, selected_positions, 1])
    obs_field[.!valid_obs.(obs_field)] .= NaN
    obs_city = fill(NaN, length(EVAL_TIMES), length(CITIES))
    for city_index in eachindex(CITIES), time_index in eachindex(EVAL_TIMES)
        obs_city[time_index, city_index] = finite_mean(@view(obs_field[time_index, city_positions[city_index]]), valid_obs)
    end
    obs_row_lookup = Dict(time_key(t) => i for (i, t) in enumerate(EVAL_TIMES))

    # Observation spatial panels retain every national-control station in the
    # source file. They are plotted as coloured station points, not interpolated.
    national_positions = findall(types_all .== "B")
    national_lons = Float64.(info.lon[national_positions])
    national_lats = Float64.(info.lat[national_positions])
    coordinate_valid = isfinite.(national_lons) .& isfinite.(national_lats)
    national_positions = national_positions[coordinate_valid]
    national_lons = national_lons[coordinate_valid]
    national_lats = national_lats[coordinate_valid]
    national_ids = ids_all[national_positions]
    national_obs_field = Float64.(observations.staDataAll[obs_indices, pm10_position, national_positions, 1])
    national_obs_field[.!valid_obs.(national_obs_field)] .= NaN

    # For every station, locate its maximum observed PM10 hour and build the
    # inclusive observation-centred window, clipped to the evaluation period at
    # either boundary.
    peak_rows = zeros(Int, length(national_ids))
    peak_window_mask = falses(length(EVAL_TIMES), length(national_ids))
    for station_index in eachindex(national_ids)
        valid_rows = findall(valid_obs, @view national_obs_field[:, station_index])
        isempty(valid_rows) && continue
        local_peak = argmax(@view national_obs_field[valid_rows, station_index])
        peak_row = valid_rows[local_peak]
        peak_rows[station_index] = peak_row
        first_row = max(1, peak_row - PEAK_HALF_WINDOW_HOURS)
        last_row = min(length(EVAL_TIMES), peak_row + PEAK_HALF_WINDOW_HOURS)
        peak_window_mask[first_row:last_row, station_index] .= true
    end
    peak_station_count = count(>(0), peak_rows)
    println("Stations with a valid observed peak: $peak_station_count/$(length(national_ids))")

    # Mean-observation benchmark for the same station-level peak-window sample
    # used by the model lead curves. The constant is the mean of all valid
    # observations inside the observation-centred windows; its plotted RMSE is
    # then computed per station and averaged across stations, exactly as for the
    # three models. Pearson R is undefined for this constant forecast because
    # its variance is zero; the correlation panel therefore shows the neutral
    # R=0 reference level as the baseline instead of inventing a Pearson value.
    baseline_mask = peak_window_mask .& valid_obs.(national_obs_field)
    baseline_reference = Float64.(national_obs_field[baseline_mask])
    isempty(baseline_reference) && error("No valid observations for the peak-window mean-observation benchmark")
    baseline_value = mean(baseline_reference)
    baseline_forecast = fill(baseline_value, length(national_ids))
    baseline_accumulator = PeakAccumulator(length(national_ids))
    for time_index in eachindex(EVAL_TIMES)
        update_peak!(baseline_accumulator, @view(national_obs_field[time_index, :]),
            baseline_forecast, @view(peak_window_mask[time_index, :]))
    end
    baseline_metrics = finalize_peak(baseline_accumulator)
    baseline_rmse = baseline_metrics.rmse
    baseline_r_reference = 0.0
    println("Peak-window mean-observation benchmark: mean=$baseline_value, RMSE=$baseline_rmse, " *
        "stations=$(baseline_metrics.rmse_station_count), paired_values=$(baseline_metrics.paired_values), " *
        "Pearson R=undefined, plotted R reference=$baseline_r_reference")

    peak_accumulators = Dict{Tuple{String, Int}, PeakAccumulator}()
    function consume_peak!(model::String, valid_time::DateTime, lead::Int, values::AbstractVector{<:Real})
        length(values) == length(national_ids) || error("$model national-station vector length changed at $valid_time")
        obs_index = get(obs_row_lookup, valid_time, 0)
        obs_index == 0 && return
        accumulator = get!(peak_accumulators, (model, lead), PeakAccumulator(length(national_ids)))
        update_peak!(accumulator, @view(national_obs_field[obs_index, :]), values,
            @view(peak_window_mask[obs_index, :]))
    end

    obs_map_fields = Vector{Vector{Float64}}()
    for valid_time in RUF_MAP_TIMES
        time_position = get(obs_row_lookup, valid_time, 0)
        time_position == 0 && error("Missing observation: valid=$(valid_time) BJT")
        push!(obs_map_fields, copy(vec(@view national_obs_field[time_position, :])))
    end

    plot_times = Dict((model, city) => DateTime[] for model in ("Aurora", "RUF") for city in CITIES)
    plot_values = Dict((model, city) => Float64[] for model in ("Aurora", "RUF") for city in CITIES)
    aurora_hourly_metric_times = Dict(city => DateTime[] for city in CITIES)
    aurora_hourly_metric_values = Dict(city => Float64[] for city in CITIES)

    national_lookup = Dict(id => index for (index, id) in enumerate(national_ids))
    all(haskey(national_lookup, id) for id in station_ids) || error("Displayed city stations are absent from the national-control set")
    selected_from_national = [national_lookup[id] for id in station_ids]

    function consume_stitched!(model::String, valid_time::DateTime, values::AbstractVector{<:Real})
        length(values) == length(station_ids) || error("$model station-vector length changed at $valid_time")
        haskey(obs_row_lookup, valid_time) || return
        for city_index in eachindex(CITIES)
            forecast_city = finite_mean(@view(values[city_positions[city_index]]), valid_forecast)
            push!(plot_times[(model, CITIES[city_index])], valid_time)
            push!(plot_values[(model, CITIES[city_index])], forecast_city)
        end
    end

    function consume_aurora_hourly_metric!(valid_time::DateTime, values::AbstractVector{<:Real})
        length(values) == length(station_ids) || error("Aurora station-vector length changed at $valid_time")
        haskey(obs_row_lookup, valid_time) || return
        for city_index in eachindex(CITIES)
            forecast_city = finite_mean(@view(values[city_positions[city_index]]), valid_forecast)
            push!(aurora_hourly_metric_times[CITIES[city_index]], valid_time)
            push!(aurora_hourly_metric_values[CITIES[city_index]], forecast_city)
        end
    end

    # RUF peak-window statistics use every native lead. The displayed rolling
    # series is taken from lead 12 of the same pass, avoiding a second read.
    println("Reading RUF cycles for peak-window lead statistics and lead-12 series")
    ruf_template_stamp = Dates.format(first(EVAL_TIMES), STAMP_FORMAT)
    ruf_template_path = joinpath(RUF_ROOT, ruf_template_stamp, "staData.$ruf_template_stamp.mat")
    isfile(ruf_template_path) || error("Missing RUF template: $ruf_template_path")
    ruf_template = staDATA(ruf_template_path)
    ruf_lookup = Dict(id => i for (i, id) in enumerate(clean_text.(ruf_template.staInfo.staID)))
    all(haskey(ruf_lookup, id) for id in national_ids) || error("RUF staData misses one or more national-control stations")
    ruf_national_positions = [ruf_lookup[id] for id in national_ids]
    ruf_cycle_first = max(EVAL_FIRST - Hour(72), ZonedDateTime(DateTime(2025, 4, 8, 0), BJT))
    ruf_cycles = collect(ruf_cycle_first:Hour(1):EVAL_LAST)
    for (cycle_index, cycle) in enumerate(ruf_cycles)
        stamp = Dates.format(cycle, STAMP_FORMAT)
        path = joinpath(RUF_ROOT, stamp, "staData.$stamp.mat")
        isfile(path) || error("Missing RUF cycle for peak-window evaluation: init=$stamp, file=$path")
        object = staDATA(path)
        pm_position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
        pm_position === nothing && error("PM10 absent in $path")
        lead_lookup = Dict(Int(lead) => i for (i, lead) in enumerate(object.leadTimes))
        for lead in 0:72
            valid_time = time_key(cycle + Hour(lead))
            time_key(EVAL_FIRST) <= valid_time <= time_key(EVAL_LAST) || continue
            lead_position = get(lead_lookup, lead, 0)
            lead_position == 0 && error("Missing RUF: init=$stamp, lead=$lead h, valid=$valid_time")
            values = Float64.(object.staDataAll[1, pm_position, ruf_national_positions, lead_position])
            consume_peak!("RUF", valid_time, lead, values)
            lead == DISPLAY_LEAD_HOURS &&
                consume_stitched!("RUF", valid_time, @view(values[selected_from_national]))
        end
        (cycle_index == 1 || cycle_index % 24 == 0 || cycle_index == length(ruf_cycles)) &&
            println("  RUF cycle $cycle_index/$(length(ruf_cycles)): $stamp")
    end

    # Reuse the exact NetCDF definitions, component sum, unit conversion and
    # BilinearSampler implementation from ThreeModelEva.
    py"""
import sys, re
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
sys.path.insert(0, "/home/zsm/IETM/code/ThreeModelEva")
from three_model_eval_netcdf_plot import NetCDF, BilinearSampler, AI_COMPONENTS

def spatial_series_aurora(root, first_iso, last_iso, station_lats, station_lons):
    first = datetime.fromisoformat(first_iso)
    last = datetime.fromisoformat(last_iso)
    series_first = first - timedelta(hours=24)
    series_last = last + timedelta(hours=24)
    lats = np.asarray(station_lats, dtype=np.float64)
    lons = np.asarray(station_lons, dtype=np.float64)
    rows_t, rows_l, rows_v, missing = [], [], [], []
    sampler = None
    first_init = (series_first - timedelta(hours=72)).date()
    last_init = series_last.date()
    day = first_init
    while day <= last_init:
        init_bjt = datetime.combine(day, datetime.min.time()).replace(hour=20)
        path = Path(root) / f"aurora_{day:%Y%m%d}T12.nc"
        selected = [lead for lead in range(0, 73, 12) if series_first <= init_bjt + timedelta(hours=lead) <= series_last]
        if selected and not path.is_file():
            for lead in selected:
                valid = init_bjt + timedelta(hours=lead)
                missing.append(f"Aurora init={init_bjt.isoformat()} lead={lead} h valid={valid.isoformat()} file={path}")
            day += timedelta(days=1)
            continue
        if selected:
            with NetCDF(path) as ds:
                if sampler is None:
                    sampler = BilinearSampler(ds.read_double("lat"), ds.read_double("lon"), lats, lons)
                lookup = {int(v): i for i, v in enumerate(ds.read_int("lead_time").tolist())}
                pm10 = ds.read_float("PM10")[0]
                for lead in selected:
                    valid = init_bjt + timedelta(hours=lead)
                    if lead not in lookup:
                        missing.append(f"Aurora init={init_bjt.isoformat()} lead={lead} h valid={valid.isoformat()} file={path}")
                        continue
                    rows_t.append(valid.isoformat())
                    rows_l.append(lead)
                    rows_v.append(sampler.sample(pm10[lookup[lead]]) * 1.0e9)
        day += timedelta(days=1)
    values = np.vstack(rows_v) if rows_v else np.empty((0, len(lats)), dtype=np.float64)
    return {"times": rows_t, "leads": np.asarray(rows_l, dtype=np.int64), "values": values, "missing": missing}

def spatial_series_ai(root, first_iso, last_iso, station_lats, station_lons):
    first = datetime.fromisoformat(first_iso)
    last = datetime.fromisoformat(last_iso)
    series_first = first - timedelta(hours=3)
    series_last = last + timedelta(hours=3)
    lats = np.asarray(station_lats, dtype=np.float64)
    lons = np.asarray(station_lons, dtype=np.float64)
    rows_t, rows_l, rows_v, missing = [], [], [], []
    sampler = None
    for directory in sorted(Path(root).glob("????????_????")):
        try:
            init_label_utc = datetime.strptime(directory.name, "%Y%m%d_%H%M")
        except ValueError:
            continue
        init_bjt = init_label_utc - timedelta(minutes=30) + timedelta(hours=8)
        selected = [lead for lead in range(3, 73, 3) if series_first <= init_bjt + timedelta(hours=lead) <= series_last]
        for lead in selected:
            valid = init_bjt + timedelta(hours=lead)
            valid_label_utc = init_label_utc + timedelta(hours=lead)
            path = directory / f"AI_GAMFS.{init_label_utc:%Y%m%d_%H%M}+{valid_label_utc:%Y%m%d_%H%M}.V01.nc"
            if not path.is_file():
                missing.append(f"AI-GAMFS init={init_bjt.isoformat()} lead={lead} h valid={valid.isoformat()} file={path}")
                continue
            with NetCDF(path) as ds:
                if sampler is None:
                    sampler = BilinearSampler(ds.read_double("latitude"), ds.read_double("longitude"), lats, lons)
                forecast = np.zeros(len(lats), dtype=np.float64)
                for component in AI_COMPONENTS:
                    forecast += sampler.sample(ds.read_float(component)[0])
                forecast *= 1.0e9
            rows_t.append(valid.isoformat())
            rows_l.append(lead)
            rows_v.append(forecast)
    values = np.vstack(rows_v) if rows_v else np.empty((0, len(lats)), dtype=np.float64)
    return {"times": rows_t, "leads": np.asarray(rows_l, dtype=np.int64), "values": values, "missing": missing}

def spatial_extract_maps(aurora_root, ai_root, aurora_specs, ai_specs):
    aurora_fields, ai_fields = [], []
    aurora_lat = aurora_lon = ai_lat = ai_lon = None
    for filename, lead in aurora_specs:
        path = Path(aurora_root) / filename
        if not path.is_file():
            raise FileNotFoundError(f"Aurora map init={filename} lead={lead} h file={path}")
        with NetCDF(path) as ds:
            if aurora_lat is None:
                aurora_lat = ds.read_double("lat")
                aurora_lon = ds.read_double("lon")
            lookup = {int(v): i for i, v in enumerate(ds.read_int("lead_time").tolist())}
            if int(lead) not in lookup:
                raise RuntimeError(f"Aurora map missing lead={lead} h in {path}")
            aurora_fields.append(ds.read_float("PM10")[0, lookup[int(lead)]] * 1.0e9)
    for dirname, filename in ai_specs:
        path = Path(ai_root) / dirname / filename
        if not path.is_file():
            raise FileNotFoundError(f"AI-GAMFS map missing lead=3 h file={path}")
        with NetCDF(path) as ds:
            if ai_lat is None:
                ai_lat = ds.read_double("latitude")
                ai_lon = ds.read_double("longitude")
            field = np.zeros((len(ai_lat), len(ai_lon)), dtype=np.float64)
            for component in AI_COMPONENTS:
                field += ds.read_float(component)[0]
            ai_fields.append(field * 1.0e9)
    return {
        "aurora_lat": aurora_lat, "aurora_lon": aurora_lon, "aurora_fields": aurora_fields,
        "ai_lat": ai_lat, "ai_lon": ai_lon, "ai_fields": ai_fields,
    }
"""

    function consume_python_series!(model::String, result, native_interval_hours::Int)
        missing = convert(Vector{String}, result["missing"])
        for item in missing
            println("MISSING: $item")
        end
        times = convert(Vector{String}, result["times"])
        valid_times = [DateTime(value, dateformat"yyyy-mm-ddTHH:MM:SS") for value in times]
        leads = Int.(vec(Array(result["leads"])))
        values = Float64.(Array(result["values"]))
        size(values, 1) == length(times) == length(leads) || error("$model extraction arrays are inconsistent")

        # For each fixed lead and station, verify at hourly valid times. Exact
        # native values are retained and interior hours are linearly interpolated.
        # Only consecutive normal updates are joined, so missing-file gaps are
        # neither bridged nor extrapolated.
        interpolated_time_count = 0
        display_interpolated_time_count = 0
        skipped_gap_count = 0
        for lead in sort(unique(leads))
            rows = findall(==(lead), leads)
            sort!(rows; by=row -> valid_times[row])
            for row in rows
                consume_peak!(model, valid_times[row], lead, @view(values[row, :]))
                if model == "Aurora" && lead == DISPLAY_LEAD_HOURS
                    displayed_values = @view values[row, selected_from_national]
                    consume_stitched!(model, valid_times[row], displayed_values)
                    consume_aurora_hourly_metric!(valid_times[row], displayed_values)
                end
            end
            for pair_index in 1:(length(rows) - 1)
                first_row = rows[pair_index]
                last_row = rows[pair_index + 1]
                first_time = valid_times[first_row]
                last_time = valid_times[last_row]
                gap_milliseconds = Dates.value(last_time - first_time)
                gap_milliseconds % 3_600_000 == 0 || error("$model has a non-hourly valid-time gap: $first_time -- $last_time")
                gap_hours = gap_milliseconds ÷ 3_600_000
                if gap_hours != native_interval_hours
                    skipped_gap_count += 1
                    continue
                end
                first_values = @view values[first_row, :]
                last_values = @view values[last_row, :]
                valid_endpoints = valid_forecast.(first_values) .& valid_forecast.(last_values)
                for offset_hour in 1:(gap_hours - 1)
                    valid_time = first_time + Hour(offset_hour)
                    haskey(obs_row_lookup, valid_time) || continue
                    weight = offset_hour / gap_hours
                    hourly_values = fill(NaN, length(national_ids))
                    hourly_values[valid_endpoints] .=
                        (1.0 - weight) .* first_values[valid_endpoints] .+
                        weight .* last_values[valid_endpoints]
                    consume_peak!(model, valid_time, lead, hourly_values)
                    if model == "Aurora" && lead == DISPLAY_LEAD_HOURS
                        consume_aurora_hourly_metric!(
                            valid_time, @view(hourly_values[selected_from_national]))
                        display_interpolated_time_count += 1
                    end
                    interpolated_time_count += 1
                end
            end
        end
        println("$model hourly verification: native interval=$(native_interval_hours) h, interpolated valid-time rows=$interpolated_time_count, skipped abnormal gaps=$skipped_gap_count")
        model == "Aurora" && println(
            "Aurora displayed lead $(DISPLAY_LEAD_HOURS) h city metrics: " *
            "interpolated hourly rows=$display_interpolated_time_count")
    end

    first_iso = time_iso(EVAL_FIRST)
    last_iso = time_iso(EVAL_LAST)
    println("Reading Aurora all leads for peak-window statistics and lead-12 series")
    consume_python_series!("Aurora", py"spatial_series_aurora"(AURORA_ROOT, first_iso, last_iso, national_lats, national_lons), 24)

    stitched_metrics = Dict{Tuple{String, String}, Any}()
    for city_index in SERIES_CITY_INDICES, model in ("Aurora", "RUF")
        city = CITIES[city_index]
        reference = Float64[]
        predict = Float64[]
        metric_times = model == "Aurora" ? aurora_hourly_metric_times[city] : plot_times[(model, city)]
        metric_values = model == "Aurora" ? aurora_hourly_metric_values[city] : plot_values[(model, city)]
        for (valid_time, forecast_city) in zip(metric_times, metric_values)
            obs_index = get(obs_row_lookup, valid_time, 0)
            obs_index == 0 && continue
            reference_city = obs_city[obs_index, city_index]
            if valid_obs(reference_city) && valid_forecast(forecast_city)
                push!(reference, reference_city)
                push!(predict, forecast_city)
            end
        end
        isempty(reference) && error("No valid stitched pairs for $model in $city")
        stitched_metrics[(model, city)] = getStatistics(reference, predict)
        metrics = stitched_metrics[(model, city)]
        println("City-series $model $(DISPLAY_LEAD_HOURS) h $city: paired_values=$(length(reference)), " *
            "RMSE=$(metrics.rmse), R=$(metrics.corr)")
    end

    model_leads = Dict(
        "RUF" => collect(0:72),
        "Aurora" => collect(0:12:72),
    )
    lead_curves = Dict{String, Tuple{Vector{Int}, Vector{Float64}, Vector{Int}}}()
    lead_correlations = Dict{String, Tuple{Vector{Int}, Vector{Float64}, Vector{Int}}}()
    peak_statistics = DataFrame(
        model=String[], lead_hour=Int[], mean_station_rmse_ug_m3=Float64[],
        mean_station_correlation_r=Float64[], rmse_station_count=Int[],
        r_station_count=Int[], paired_values=Int[], evaluation_time_step_hour=Int[],
        temporal_interpolation=String[], native_valid_interval_hour=Int[], peak_window_first_offset_hour=Int[],
        peak_window_last_offset_hour=Int[], evaluation_start_bjt=String[], evaluation_end_bjt=String[],
    )
    for model in ("RUF", "Aurora")
        leads = model_leads[model]
        rmse_values = Float64[]
        r_values = Float64[]
        rmse_counts = Int[]
        r_counts = Int[]
        for lead in leads
            accumulator = get(peak_accumulators, (model, lead), PeakAccumulator(length(national_ids)))
            metrics = finalize_peak(accumulator)
            push!(rmse_values, metrics.rmse)
            push!(r_values, metrics.corr)
            push!(rmse_counts, metrics.rmse_station_count)
            push!(r_counts, metrics.r_station_count)
            push!(peak_statistics, (
                model, lead, metrics.rmse, metrics.corr, metrics.rmse_station_count,
                metrics.r_station_count, metrics.paired_values, 1,
                model == "RUF" ? "native hourly" : "linear; no extrapolation or missing-gap bridging",
                model == "RUF" ? 1 : model == "Aurora" ? 24 : 3, -PEAK_HALF_WINDOW_HOURS,
                PEAK_HALF_WINDOW_HOURS, first_iso, last_iso,
            ))
        end
        lead_curves[model] = (leads, rmse_values, rmse_counts)
        lead_correlations[model] = (leads, r_values, r_counts)
        println("Peak-window $model: RMSE stations $(minimum(rmse_counts))-$(maximum(rmse_counts)); R stations $(minimum(r_counts))-$(maximum(r_counts))")
    end
    peak_statistics_path = joinpath(OUTPUT_DIR, "$(OUTPUT_BASE)_lead_statistics.csv")
    CSV.write(peak_statistics_path, peak_statistics)
    println("Wrote: $peak_statistics_path")

    println("Reading three RUF lead-12 concentration grids")
    ruf_lon, ruf_lat = load(RUF_GRID, "lonM", "latM")
    ruf_fields = Matrix{Float64}[]
    for valid_time in RUF_MAP_TIMES
        cycle_time = valid_time - Hour(DISPLAY_LEAD_HOURS)
        stamp = Dates.format(cycle_time, STAMP_FORMAT)
        ctl = joinpath(RUF_ROOT, stamp, "ruf.ctl")
        isfile(ctl) || error("Missing RUF map: init=$stamp, lead=$(DISPLAY_LEAD_HOURS) h, valid=$valid_time, ctl=$ctl")
        grid = gradData(ctl)
        setTime!(grid, ZonedDateTime(valid_time, BJT))
        field = Float64.(read(grid, "PM10", 1))
        field[.!isfinite.(field) .| (field .< 0)] .= NaN
        push!(ruf_fields, field)
    end

    aurora_specs_python = [(item[1], item[2]) for item in AURORA_MAP_SPECS]
    println("Reading Aurora direct concentration grids")
    map_result = py"spatial_extract_maps"(AURORA_ROOT, "", aurora_specs_python, Tuple{String, String}[])
    aurora_lat = Float64.(vec(Array(map_result["aurora_lat"])))
    aurora_lon = Float64.(vec(Array(map_result["aurora_lon"])))
    aurora_fields = [Float64.(Array(field)) for field in map_result["aurora_fields"]]
    for field in aurora_fields
        field[.!isfinite.(field) .| (field .< 0)] .= NaN
    end

    println("Drawing composite spatial and city-series figure")
    matplotlib = pyimport("matplotlib")
    configure_matplotlib!(matplotlib)
    plt = pyimport("matplotlib.pyplot")
    mdates = pyimport("matplotlib.dates")
    ticker = matplotlib.ticker
    colors_module = matplotlib.colors
    ccrs = pyimport("cartopy.crs")
    shpreader = pyimport("cartopy.io.shapereader")
    fm = matplotlib.font_manager

    font_regular_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_bold_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    figure_fonts = font_properties(fm, font_regular_path, font_bold_path)
    font_regular = figure_fonts.axis_label
    font_small = figure_fonts.annotation
    font_series_legend = figure_fonts.figure_legend
    font_tick = figure_fonts.tick
    font_map_tick = figure_fonts.tick
    font_metric = figure_fonts.annotation
    font_axis_label = figure_fonts.axis_label
    font_series_ylabel = figure_fonts.axis_label
    font_panel_label = fm.FontProperties(fname=font_bold_path, size=8.5)
    font_row_label = figure_fonts.row_label

    colors = Dict("Observation" => OBS_COLOR, "Aurora" => AURORA_COLOR,
        "RUF" => RUF_COLOR)
    figure = plt.figure(figsize=(FIGURE_WIDTH, 3.85), dpi=600)

    # The left half is a 3 x 3 spatial matrix. Peak-window lead-time panels are
    # intentionally omitted from this figure; their statistics remain in CSV.
    left_x0, left_total_width = 0.055, 0.555
    # (a)-(i) Three direct-concentration fields for RUF and Aurora plus
    # national-control station observations. All panels use
    # one color normalization so spatial intensities remain directly comparable.
    plate = ccrs.PlateCarree()
    map_rows = [
        ("RUF", ruf_lon, ruf_lat, ruf_fields, RUF_MAP_TIMES),
        ("Aurora", aurora_lon, aurora_lat, aurora_fields, [item[3] for item in AURORA_MAP_SPECS]),
        ("Observation", national_lons, national_lats, obs_map_fields, RUF_MAP_TIMES),
    ]
    map_x0, map_total_width, map_gap = left_x0, left_total_width, 0.0
    map_width = (map_total_width - 2 * map_gap) / 3
    map_bottom, map_top = 0.07, 0.975
    map_height = (map_top - map_bottom) / 3
    map_y = [map_bottom + 2 * map_height, map_bottom + map_height, map_bottom]
    boundary_file = "/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf"
    cmap = plt.get_cmap("afmhot_r")
    norm = colors_module.Normalize(vmin=0.0, vmax=1000.0)
    map_handle = nothing
    for (row_index, (model, lon, lat, fields, times)) in enumerate(map_rows)
        for column_index in 1:3
            x = map_x0 + (column_index - 1) * (map_width + map_gap)
            axis = figure.add_axes([x, map_y[row_index], map_width, map_height], projection=plate)
            field = fields[column_index]
            if model == "Observation"
                valid = isfinite.(field)
                map_handle = axis.scatter(lon[valid], lat[valid], c=field[valid], transform=plate,
                    cmap=cmap, norm=norm, s=2.2, marker="o", edgecolors="none", zorder=2)
            else
                map_handle = axis.pcolormesh(lon, lat, field, transform=plate, cmap=cmap, norm=norm,
                    shading="auto", rasterized=true, zorder=1)
            end
            axis.add_geometries(shpreader.Reader(boundary_file).geometries(), plate,
                facecolor="none", edgecolor="#555555", linewidth=0.6, zorder=3)
            axis.set_extent([73, 136, 16, 55], crs=plate)
            axis.set_aspect("auto")
            gl = axis.gridlines(crs=plate, draw_labels=true, linewidth=GRID_LINEWIDTH,
                color="#777777", alpha=0.55, linestyle="--", zorder=2)
            gl.top_labels = false
            gl.right_labels = false
            gl.left_labels = column_index == 1
            gl.bottom_labels = row_index == 3
            gl.xlocator = ticker.FixedLocator([80, 100, 120])
            gl.ylocator = ticker.FixedLocator([20, 30, 40, 50])
            gl.xlabel_style = Dict("fontproperties" => font_map_tick)
            gl.ylabel_style = Dict("fontproperties" => font_map_tick)
            gl.xpadding = 1.0
            gl.ypadding = 1.0
            time_label = Dates.format(times[column_index], dateformat"mm-dd HH:MM")
            panel_index = 3 * (row_index - 1) + column_index
            panel_letter = Char('a' + panel_index - 1)
            title = "($panel_letter) $time_label BJT"
            axis.text(0.015, 0.975, title, transform=axis.transAxes, va="top", ha="left",
                fontproperties=font_panel_label, color="black", zorder=5,
                bbox=Dict("facecolor" => "white", "edgecolor" => "none", "alpha" => 0.72, "pad" => 0.5))
            if column_index == 1
                row_label = model == "Observation" ? "obs" : model
                axis.text(-0.22, 0.5, row_label, transform=axis.transAxes, rotation=90,
                    va="center", ha="center", color="black", fontproperties=font_row_label)
            end

            # South China Sea inset: its lower-right corner coincides exactly
            # with the parent map's lower-right corner (no positional padding).
            inset_width = map_width * 0.17
            inset_height = map_height * 0.26
            inset_axis = figure.add_axes([
                x + map_width - inset_width,
                map_y[row_index],
                inset_width,
                inset_height,
            ], projection=plate)
            if model == "Observation"
                valid = isfinite.(field)
                inset_axis.scatter(lon[valid], lat[valid], c=field[valid], transform=plate,
                    cmap=cmap, norm=norm, s=1.5, marker="o", edgecolors="none", zorder=2)
            else
                inset_axis.pcolormesh(lon, lat, field, transform=plate, cmap=cmap, norm=norm,
                    shading="auto", rasterized=true, zorder=1)
            end
            inset_axis.add_geometries(shpreader.Reader(boundary_file).geometries(), plate,
                facecolor="none", edgecolor="#555555", linewidth=0.6, zorder=3)
            inset_axis.set_extent([106, 120, 2, 24], crs=plate)
            inset_axis.set_aspect("auto")
            inset_axis.set_xticks([])
            inset_axis.set_yticks([])
            inset_axis.tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
            for spine in inset_axis.spines.values()
                spine.set_linewidth(0.6)
            end
            inset_axis.set_zorder(6)
        end
    end
    colorbar_axis = figure.add_axes([map_x0 + 0.095, 0.012, map_total_width - 0.19, 0.010])
    colorbar = figure.colorbar(map_handle, cax=colorbar_axis, orientation="horizontal", extend="max")
    colorbar.set_label(raw"PM$_{10}$ ($\mathrm{\mu}$g m$^{-3}$)", fontproperties=font_regular, labelpad=1.0)
    colorbar.set_ticks([0, 200, 400, 600, 800, 1000])
    colorbar.ax.tick_params(labelsize=TICK_SIZE, width=AXIS_LINEWIDTH, length=2.8, pad=1.0)
    for label in colorbar.ax.get_xticklabels()
        label.set_fontproperties(font_map_tick)
    end

    # (j)-(m) City-average time series, preserving the original all-national-
    # station averaging and fastest-update stitching rules.
    series_axes = PyObject[]
    right_x, right_width = 0.680, 0.305
    right_bottom, right_top = 0.045, 0.970
    right_height = (right_top - right_bottom) / 4
    right_y = [right_bottom + 3 * right_height, right_bottom + 2 * right_height,
        right_bottom + right_height, right_bottom]
    obs_dates = DateTime.(EVAL_TIMES)
    series_handles = nothing
    series_labels = nothing
    for (panel_index, city_index) in enumerate(SERIES_CITY_INDICES)
        city = CITIES[city_index]
        axis = figure.add_axes([right_x, right_y[panel_index], right_width, right_height])
        push!(series_axes, axis)
        axis.bar(obs_dates, obs_city[:, city_index], width=0.02, color=colors["Observation"],
            edgecolor="none", label="Obs", zorder=1)
        for (model, marker) in (("Aurora", "s"), ("RUF", nothing))
            times = plot_times[(model, city)]
            values = plot_values[(model, city)]
            order = sortperm(times)
            kwargs = Dict{Symbol, Any}(:color => colors[model], :linewidth => DATA_LINEWIDTH,
                :linestyle => "-", :label => model, :zorder => 3)
            if marker !== nothing
                kwargs[:marker] = marker
                kwargs[:markersize] = 2.8
                kwargs[:markeredgewidth] = 0.6
            end
            axis.plot(times[order], values[order]; kwargs...)
        end
        letter = Char('j' + panel_index - 1)
        axis.text(0.015, 0.965, "($letter) $(CITY_ENGLISH[city_index])", transform=axis.transAxes,
            va="top", ha="left", fontproperties=font_panel_label)
        metric_label_x = city == "海口" ? 0.04 : city == "太原" ? 0.47 : 0.57
        metric_value_x = city == "海口" ? (0.25, 0.43) :
            city == "太原" ? (0.74, 0.89) : (0.78, 0.91)
        rmse_label_y = city == "海口" ? 0.46 : 0.90
        rmse_value_y = rmse_label_y
        corr_label_y = city == "海口" ? 0.30 : 0.74
        corr_value_y = corr_label_y
        axis.text(metric_label_x, rmse_label_y, "RMSE:", transform=axis.transAxes, va="center", ha="left",
            color="black", fontproperties=font_metric, zorder=10)
        axis.text(metric_label_x, corr_label_y, "R:", transform=axis.transAxes, va="center", ha="left",
            color="black", fontproperties=font_metric, zorder=10)
        for (model_index, model) in enumerate(("Aurora", "RUF"))
            metrics = stitched_metrics[(model, city)]
            axis.text(metric_value_x[model_index], rmse_value_y, string(round(metrics.rmse; digits=1)),
                transform=axis.transAxes, va="center", ha="center", color=colors[model],
                fontproperties=font_metric, zorder=10)
            axis.text(metric_value_x[model_index], corr_value_y, string(round(metrics.corr; digits=2)),
                transform=axis.transAxes, va="center", ha="center", color=colors[model],
                fontproperties=font_metric, zorder=10)
        end
        axis.set_xlim(time_key(EVAL_FIRST), time_key(EVAL_LAST))
        axis.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        axis.grid(true, color="#D9D9D9", linewidth=GRID_LINEWIDTH, linestyle="--", alpha=0.75)
        axis.tick_params(axis="both", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH, length=2.8, pad=1.6)
        for label in vcat(axis.get_xticklabels(), axis.get_yticklabels())
            label.set_fontproperties(font_map_tick)
        end
        panel_index < length(SERIES_CITY_INDICES) && axis.tick_params(labelbottom=false)
        if panel_index == 1
            series_handles, series_labels = axis.get_legend_handles_labels()
        end
    end
    figure.text(0.635, (right_bottom + right_top) / 2, raw"PM$_{10}$ ($\mathrm{\mu}$g m$^{-3}$)",
        rotation=90, va="center", ha="center", fontproperties=font_series_ylabel)
    original_labels = convert(Vector{String}, series_labels)
    legend_order = [findfirst(==(label), original_labels) for label in ("Obs", "Aurora", "RUF")]
    any(isnothing, legend_order) && error("Time-series legend entries are incomplete")
    ordered_handles = [series_handles[index] for index in Int.(legend_order)]
    ordered_labels = [original_labels[index] for index in Int.(legend_order)]
    figure.legend(ordered_handles, ordered_labels, loc="lower center", bbox_to_anchor=(0.835, -0.065),
        ncol=3, frameon=false, prop=font_series_legend, handlelength=1.8, columnspacing=0.8)

    png_path = joinpath(OUTPUT_DIR, "$OUTPUT_BASE.png")
    pdf_path = joinpath(OUTPUT_DIR, "$OUTPUT_BASE.pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=0.025)
    figure.savefig(pdf_path, dpi=600, bbox_inches="tight", pad_inches=0.025)
    plt.close(figure)
    println("Wrote: $png_path")
    println("Wrote: $pdf_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
