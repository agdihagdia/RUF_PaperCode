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
const OUTPUT_BASE = "PaperFigure_AppendixCompare3ModelSpatial"
const STAMP_FORMAT = dateformat"yyyymmddHH"

# Aurora uses lead 12 h. RUF lead 1 and observations use the exact Aurora valid
# times. AI-GAMFS uses its direct lead-3 output one hour later; no temporal or
# spatial interpolation is performed and every panel states its actual time.
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
    national_city_names = cities_all[national_positions]
    coordinate_valid = isfinite.(national_lons) .& isfinite.(national_lats)
    national_positions = national_positions[coordinate_valid]
    national_lons = national_lons[coordinate_valid]
    national_lats = national_lats[coordinate_valid]
    national_city_names = national_city_names[coordinate_valid]
    national_cities = sort(unique(filter(value -> !isempty(value), national_city_names)))
    national_city_groups = [findall(==(city), national_city_names) for city in national_cities]
    national_obs_field = Float64.(observations.staDataAll[obs_indices, pm10_position, national_positions, 1])
    national_obs_field[.!valid_obs.(national_obs_field)] .= NaN
    national_obs_city = fill(NaN, length(EVAL_TIMES), length(national_cities))
    for city_index in eachindex(national_cities), time_index in eachindex(EVAL_TIMES)
        national_obs_city[time_index, city_index] = finite_mean(
            @view(national_obs_field[time_index, national_city_groups[city_index]]), valid_obs,
        )
    end
    println("National-control cities available for the observation benchmark: $(length(national_cities))")
    obs_map_fields = Vector{Vector{Float64}}()
    for valid_time in RUF_MAP_TIMES
        time_position = get(obs_row_lookup, valid_time, 0)
        time_position == 0 && error("Missing observation: valid=$(valid_time) BJT")
        push!(obs_map_fields, copy(vec(@view national_obs_field[time_position, :])))
    end

    plot_times = Dict((model, city) => DateTime[] for model in ("AI-GAMFS", "Aurora", "RUF") for city in CITIES)
    plot_values = Dict((model, city) => Float64[] for model in ("AI-GAMFS", "Aurora", "RUF") for city in CITIES)

    function consume_stitched!(model::String, valid_time::DateTime, values::AbstractVector{<:Real})
        length(values) == length(station_ids) || error("$model station-vector length changed at $valid_time")
        haskey(obs_row_lookup, valid_time) || return
        for city_index in eachindex(CITIES)
            forecast_city = finite_mean(@view(values[city_positions[city_index]]), valid_forecast)
            push!(plot_times[(model, CITIES[city_index])], valid_time)
            push!(plot_values[(model, CITIES[city_index])], forecast_city)
        end
    end

    # RUF freshest-update series: the cycle valid at each hour, lead 0.
    println("Reading RUF hourly lead-0 station forecasts")
    ruf_template_stamp = Dates.format(first(EVAL_TIMES), STAMP_FORMAT)
    ruf_template_path = joinpath(RUF_ROOT, ruf_template_stamp, "staData.$ruf_template_stamp.mat")
    isfile(ruf_template_path) || error("Missing RUF template: $ruf_template_path")
    ruf_template = staDATA(ruf_template_path)
    ruf_lookup = Dict(id => i for (i, id) in enumerate(clean_text.(ruf_template.staInfo.staID)))
    all(haskey(ruf_lookup, id) for id in station_ids) || error("RUF staData misses selected national-control stations")
    ruf_positions = [ruf_lookup[id] for id in station_ids]
    for (time_index, valid_zdt) in enumerate(EVAL_TIMES)
        stamp = Dates.format(valid_zdt, STAMP_FORMAT)
        path = joinpath(RUF_ROOT, stamp, "staData.$stamp.mat")
        isfile(path) || error("Missing RUF: init=$stamp, lead=0 h, valid=$(time_iso(valid_zdt)), file=$path")
        object = staDATA(path)
        pm_position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
        pm_position === nothing && error("PM10 absent in $path")
        lead_position = findfirst(==(0), Int.(object.leadTimes))
        lead_position === nothing && error("Missing RUF: init=$stamp, lead=0 h, valid=$(time_iso(valid_zdt))")
        values = Float64.(object.staDataAll[1, pm_position, ruf_positions, lead_position])
        consume_stitched!("RUF", time_key(valid_zdt), values)
        (time_index == 1 || time_index % 24 == 0 || time_index == length(EVAL_TIMES)) &&
            println("  RUF $time_index/$(length(EVAL_TIMES)): $stamp")
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
    lats = np.asarray(station_lats, dtype=np.float64)
    lons = np.asarray(station_lons, dtype=np.float64)
    rows_t, rows_v, missing = [], [], []
    sampler = None
    first_init = (first - timedelta(hours=12)).date()
    last_init = last.date()
    day = first_init
    while day <= last_init:
        init_bjt = datetime.combine(day, datetime.min.time()).replace(hour=20)
        path = Path(root) / f"aurora_{day:%Y%m%d}T12.nc"
        for lead in (0, 12):
            valid = init_bjt + timedelta(hours=lead)
            if not (first <= valid <= last):
                continue
            if not path.is_file():
                missing.append(f"Aurora init={init_bjt.isoformat()} lead={lead} h valid={valid.isoformat()} file={path}")
                continue
            with NetCDF(path) as ds:
                if sampler is None:
                    sampler = BilinearSampler(ds.read_double("lat"), ds.read_double("lon"), lats, lons)
                lookup = {int(v): i for i, v in enumerate(ds.read_int("lead_time").tolist())}
                if lead not in lookup:
                    missing.append(f"Aurora init={init_bjt.isoformat()} lead={lead} h valid={valid.isoformat()} file={path}")
                    continue
                field = ds.read_float("PM10")[0, lookup[lead]]
                rows_t.append(valid.isoformat())
                rows_v.append(sampler.sample(field) * 1.0e9)
        day += timedelta(days=1)
    values = np.vstack(rows_v) if rows_v else np.empty((0, len(lats)), dtype=np.float64)
    return {"times": rows_t, "values": values, "missing": missing}

def spatial_series_ai(root, first_iso, last_iso, station_lats, station_lons):
    first = datetime.fromisoformat(first_iso)
    last = datetime.fromisoformat(last_iso)
    lats = np.asarray(station_lats, dtype=np.float64)
    lons = np.asarray(station_lons, dtype=np.float64)
    rows_t, rows_v, missing = [], [], []
    sampler = None
    for directory in sorted(Path(root).glob("????????_????")):
        try:
            init_label_utc = datetime.strptime(directory.name, "%Y%m%d_%H%M")
        except ValueError:
            continue
        init_bjt = init_label_utc - timedelta(minutes=30) + timedelta(hours=8)
        lead = 3
        valid = init_bjt + timedelta(hours=lead)
        if not (first <= valid <= last):
            continue
        valid_label_utc = init_label_utc + timedelta(hours=lead)
        path = directory / f"AI_GAMFS.{init_label_utc:%Y%m%d_%H%M}+{valid_label_utc:%Y%m%d_%H%M}.V01.nc"
        if not path.is_file():
            missing.append(f"AI-GAMFS init={init_bjt.isoformat()} lead=3 h valid={valid.isoformat()} file={path}")
            continue
        with NetCDF(path) as ds:
            if sampler is None:
                sampler = BilinearSampler(ds.read_double("latitude"), ds.read_double("longitude"), lats, lons)
            forecast = np.zeros(len(lats), dtype=np.float64)
            for component in AI_COMPONENTS:
                forecast += sampler.sample(ds.read_float(component)[0])
            forecast *= 1.0e9
        rows_t.append(valid.isoformat())
        rows_v.append(forecast)
    values = np.vstack(rows_v) if rows_v else np.empty((0, len(lats)), dtype=np.float64)
    return {"times": rows_t, "values": values, "missing": missing}

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

    function consume_python_series!(model::String, result)
        missing = convert(Vector{String}, result["missing"])
        for item in missing
            println("MISSING: $item")
        end
        times = convert(Vector{String}, result["times"])
        values = Float64.(Array(result["values"]))
        size(values, 1) == length(times) || error("$model extraction arrays are inconsistent")
        for row in eachindex(times)
            consume_stitched!(model, DateTime(times[row], dateformat"yyyy-mm-ddTHH:MM:SS"), @view(values[row, :]))
        end
    end

    first_iso = time_iso(EVAL_FIRST)
    last_iso = time_iso(EVAL_LAST)
    println("Reading Aurora lead-0/12 stitched station forecasts")
    consume_python_series!("Aurora", py"spatial_series_aurora"(AURORA_ROOT, first_iso, last_iso, station_lats, station_lons))
    println("Reading AI-GAMFS lead-3 stitched station forecasts")
    consume_python_series!("AI-GAMFS", py"spatial_series_ai"(AI_ROOT, first_iso, last_iso, station_lats, station_lons))

    stitched_metrics = Dict{Tuple{String, String}, Any}()
    for city_index in SERIES_CITY_INDICES, model in ("AI-GAMFS", "Aurora", "RUF")
        city = CITIES[city_index]
        reference = Float64[]
        predict = Float64[]
        for (valid_time, forecast_city) in zip(plot_times[(model, city)], plot_values[(model, city)])
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
    end

    # Mean-observation benchmark over every valid national-control city-hour in
    # the evaluation window. Stations are averaged within each city first, so
    # cities (rather than individual stations) are the statistical units.
    baseline_reference = Float64.(national_obs_city[valid_obs.(national_obs_city)])
    isempty(baseline_reference) && error("No valid observations for the mean-observation benchmark")
    baseline_value = mean(baseline_reference)
    baseline_predict = fill(baseline_value, length(baseline_reference))
    baseline_metrics = getStatistics(baseline_reference, baseline_predict)
    baseline_rmse = baseline_metrics.rmse
    baseline_r = NaN  # Pearson R is undefined for any constant forecast.
    println("Mean-observation benchmark: cities=$(length(national_cities)), mean=$baseline_value, n=$(length(baseline_reference)), RMSE=$baseline_rmse, R=$baseline_r (undefined for zero forecast variance)")

    isfile(LEAD_STATISTICS) || error("Missing audited lead statistics: $LEAD_STATISTICS")
    lead_table = CSV.read(LEAD_STATISTICS, DataFrame)
    lead_curves = Dict(model => pooled_lead_rmse(lead_table, model) for model in ("RUF", "Aurora", "AI-GAMFS"))
    lead_correlations = Dict(model => pooled_lead_correlation(lead_table, model) for model in ("RUF", "Aurora", "AI-GAMFS"))

    println("Reading three RUF lead-1 concentration grids")
    ruf_lon, ruf_lat = load(RUF_GRID, "lonM", "latM")
    ruf_fields = Matrix{Float64}[]
    for valid_time in RUF_MAP_TIMES
        cycle_time = valid_time - Hour(1)
        stamp = Dates.format(cycle_time, STAMP_FORMAT)
        ctl = joinpath(RUF_ROOT, stamp, "ruf.ctl")
        isfile(ctl) || error("Missing RUF map: init=$stamp, lead=1 h, valid=$valid_time, ctl=$ctl")
        grid = gradData(ctl)
        setTime!(grid, ZonedDateTime(valid_time, BJT))
        field = Float64.(read(grid, "PM10", 1))
        field[.!isfinite.(field) .| (field .< 0)] .= NaN
        push!(ruf_fields, field)
    end

    aurora_specs_python = [(item[1], item[2]) for item in AURORA_MAP_SPECS]
    ai_specs_python = [(item[1], item[2]) for item in AI_MAP_SPECS]
    println("Reading Aurora and AI-GAMFS direct concentration grids")
    map_result = py"spatial_extract_maps"(AURORA_ROOT, AI_ROOT, aurora_specs_python, ai_specs_python)
    aurora_lat = Float64.(vec(Array(map_result["aurora_lat"])))
    aurora_lon = Float64.(vec(Array(map_result["aurora_lon"])))
    aurora_fields = [Float64.(Array(field)) for field in map_result["aurora_fields"]]
    ai_lat = Float64.(vec(Array(map_result["ai_lat"])))
    ai_lon = Float64.(vec(Array(map_result["ai_lon"])))
    ai_fields = [Float64.(Array(field)) for field in map_result["ai_fields"]]
    for field in vcat(aurora_fields, ai_fields)
        field[.!isfinite.(field) .| (field .< 0)] .= NaN
    end

    println("Drawing composite lead-time, spatial and city-series figure")
    matplotlib = pyimport("matplotlib")
    plt = pyimport("matplotlib.pyplot")
    mdates = pyimport("matplotlib.dates")
    ticker = matplotlib.ticker
    colors_module = matplotlib.colors
    ccrs = pyimport("cartopy.crs")
    shpreader = pyimport("cartopy.io.shapereader")
    fm = matplotlib.font_manager

    font_regular_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_bold_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    font_regular = fm.FontProperties(fname=font_regular_path, size=6)
    font_small = fm.FontProperties(fname=font_regular_path, size=5)
    font_tick = fm.FontProperties(fname=font_regular_path, size=7)
    font_map_tick = fm.FontProperties(fname=font_regular_path, size=6)
    font_metric = fm.FontProperties(fname=font_regular_path, size=5.9)
    font_axis_label = fm.FontProperties(fname=font_regular_path, size=7.5)
    font_series_ylabel = fm.FontProperties(fname=font_regular_path, size=6.5)
    font_panel_label = fm.FontProperties(fname=font_bold_path, size=6.5)
    font_row_label = fm.FontProperties(fname=font_bold_path, size=8.0)
    matplotlib.rcParams["axes.linewidth"] = 0.6
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42

    colors = Dict("Observation" => "deepskyblue", "AI-GAMFS" => "#2ca02c", "Aurora" => "#F58518", "RUF" => "#C44E52")
    figure = plt.figure(figsize=(12.0, 8.0), dpi=600)

    # The lower 1 x 2 lead-time row and upper 4 x 3 spatial matrix use exactly
    # the same narrowed left and right boundaries.
    left_x0, left_total_width = 0.055, 0.555
    lead_gap = 0.032
    lead_width = (left_total_width - lead_gap) / 2
    lead_y, lead_height = 0.045, 0.20
    rmse_axis = figure.add_axes([left_x0, lead_y, lead_width, lead_height])
    for (model, marker) in (("RUF", nothing), ("Aurora", "s"), ("AI-GAMFS", "o"))
        leads, values, _ = lead_curves[model]
        kwargs = Dict{Symbol, Any}(:color => colors[model], :linewidth => 1.15, :label => model, :zorder => 3)
        if marker !== nothing
            kwargs[:marker] = marker
            kwargs[:markersize] = 2.4
            kwargs[:markeredgewidth] = 0.0
        end
        rmse_axis.plot(leads, values; kwargs...)
    end
    baseline_label = "Mean Obs ($(round(baseline_value; digits=1)) " * raw"$\mu$g m$^{-3}$" * ")"
    rmse_axis.plot([0, 72], [baseline_rmse, baseline_rmse], color="#4D4D4D",
        linewidth=1.1, linestyle="--", label=baseline_label, zorder=2)
    rmse_axis.text(0.015, 0.965, "(a) Lead-time RMSE", transform=rmse_axis.transAxes,
        va="top", ha="left", fontproperties=font_panel_label)
    rmse_axis.set_xlim(0, 72)
    rmse_axis.set_xticks(collect(0:12:72))
    rmse_axis.set_xlabel("Lead time (h)", fontproperties=font_regular, labelpad=1.5)
    rmse_axis.set_ylabel(raw"RMSE ($\mu$g m$^{-3}$)", fontproperties=font_axis_label, labelpad=1.5)
    rmse_axis.grid(true, color="#D9D9D9", linewidth=0.4, linestyle="--", alpha=0.8)
    rmse_axis.tick_params(axis="both", labelsize=6.5, width=0.5, length=2.4, pad=1.5)
    rmse_axis.legend(loc="best", frameon=false, prop=font_small, ncol=1,
        handlelength=1.8, labelspacing=0.25)
    for label in vcat(rmse_axis.get_xticklabels(), rmse_axis.get_yticklabels())
        label.set_fontproperties(font_tick)
    end

    correlation_axis = figure.add_axes([left_x0 + lead_width + lead_gap, lead_y, lead_width, lead_height])
    for (model, marker) in (("RUF", nothing), ("Aurora", "s"), ("AI-GAMFS", "o"))
        leads, values, _ = lead_correlations[model]
        kwargs = Dict{Symbol, Any}(:color => colors[model], :linewidth => 1.15, :label => model, :zorder => 3)
        if marker !== nothing
            kwargs[:marker] = marker
            kwargs[:markersize] = 2.4
            kwargs[:markeredgewidth] = 0.0
        end
        correlation_axis.plot(leads, values; kwargs...)
    end
    correlation_axis.axhline(0.0, color="#777777", linewidth=0.45, linestyle="--", zorder=1)
    correlation_axis.text(0.015, 0.965, "(b) Lead-time R", transform=correlation_axis.transAxes,
        va="top", ha="left", fontproperties=font_panel_label)
    baseline_r_note = "Mean Obs ($(round(baseline_value; digits=1)) " * raw"$\mu$g m$^{-3}$" * "): R undefined"
    correlation_axis.text(0.98, 0.08, baseline_r_note,
        transform=correlation_axis.transAxes, va="bottom", ha="right",
        color="#4D4D4D", fontproperties=font_small,
        bbox=Dict("facecolor" => "white", "edgecolor" => "none", "alpha" => 0.8, "pad" => 0.6), zorder=5)
    correlation_axis.set_xlim(0, 72)
    correlation_axis.set_xticks(collect(0:12:72))
    correlation_axis.set_ylim(-0.4, 1.0)
    correlation_axis.set_yticks(collect(-0.4:0.2:1.0))
    correlation_axis.set_xlabel("Lead time (h)", fontproperties=font_regular, labelpad=1.5)
    correlation_axis.set_ylabel("Correlation coefficient (R)", fontproperties=font_axis_label, labelpad=0.5)
    correlation_axis.grid(true, color="#D9D9D9", linewidth=0.4, linestyle="--", alpha=0.8)
    correlation_axis.tick_params(axis="both", labelsize=6.5, width=0.5, length=2.4, pad=1.5)
    for label in vcat(correlation_axis.get_xticklabels(), correlation_axis.get_yticklabels())
        label.set_fontproperties(font_tick)
    end

    # (b)-(m) Three direct-concentration fields for each model plus national-
    # control station observations. All panels use
    # one color normalization so spatial intensities remain directly comparable.
    plate = ccrs.PlateCarree()
    map_rows = [
        ("RUF", ruf_lon, ruf_lat, ruf_fields, RUF_MAP_TIMES),
        ("Aurora", aurora_lon, aurora_lat, aurora_fields, [item[3] for item in AURORA_MAP_SPECS]),
        ("AI-GAMFS", ai_lon, ai_lat, ai_fields, [item[3] for item in AI_MAP_SPECS]),
        ("Observation", national_lons, national_lats, obs_map_fields, RUF_MAP_TIMES),
    ]
    map_x0, map_total_width, map_gap = left_x0, left_total_width, 0.0
    map_width = (map_total_width - 2 * map_gap) / 3
    map_y = [0.805, 0.645, 0.485, 0.325]
    map_height = 0.16
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
                facecolor="none", edgecolor="#555555", linewidth=0.35, zorder=3)
            axis.set_extent([73, 136, 16, 55], crs=plate)
            axis.set_aspect("auto")
            gl = axis.gridlines(crs=plate, draw_labels=true, linewidth=0.25,
                color="#777777", alpha=0.55, linestyle="--", zorder=2)
            gl.top_labels = false
            gl.right_labels = false
            gl.left_labels = column_index == 1
            gl.bottom_labels = row_index == 4
            gl.xlocator = ticker.FixedLocator([80, 100, 120])
            gl.ylocator = ticker.FixedLocator([20, 30, 40, 50])
            gl.xlabel_style = Dict("fontproperties" => font_map_tick)
            gl.ylabel_style = Dict("fontproperties" => font_map_tick)
            gl.xpadding = 1.0
            gl.ypadding = 1.0
            time_label = Dates.format(times[column_index], dateformat"mm-dd HH:MM")
            title = "$time_label BJT"
            axis.text(0.015, 0.975, title, transform=axis.transAxes, va="top", ha="left",
                fontproperties=font_panel_label, color="black", zorder=5,
                bbox=Dict("facecolor" => "white", "edgecolor" => "none", "alpha" => 0.72, "pad" => 0.5))
            if column_index == 1
                axis.text(-0.105, 0.5, model, transform=axis.transAxes, rotation=90,
                    va="center", ha="center", color="black", fontproperties=font_row_label)
            end

            # South China Sea inset: its lower-right corner coincides exactly
            # with the parent map's lower-right corner (no positional padding).
            inset_width = map_width * 0.17
            inset_height = map_height * 0.42
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
                facecolor="none", edgecolor="#555555", linewidth=0.3, zorder=3)
            inset_axis.set_extent([106, 120, 2, 24], crs=plate)
            inset_axis.set_aspect("auto")
            inset_axis.set_xticks([])
            inset_axis.set_yticks([])
            inset_axis.tick_params(left=false, bottom=false, labelleft=false, labelbottom=false)
            for spine in inset_axis.spines.values()
                spine.set_linewidth(0.45)
            end
            inset_axis.set_zorder(6)
        end
    end
    colorbar_axis = figure.add_axes([map_x0 + 0.095, 0.295, map_total_width - 0.19, 0.012])
    colorbar = figure.colorbar(map_handle, cax=colorbar_axis, orientation="horizontal", extend="max")
    colorbar.set_label(raw"PM$_{10}$ ($\mu$g m$^{-3}$)", fontproperties=font_regular, labelpad=1.0)
    colorbar.set_ticks([0, 200, 400, 600, 800, 1000])
    colorbar.ax.tick_params(labelsize=6, width=0.4, length=2.0, pad=1.0)
    for label in colorbar.ax.get_xticklabels()
        label.set_fontproperties(font_map_tick)
    end

    # (n)-(q) City-average time series, preserving the original all-national-
    # station averaging and fastest-update stitching rules.
    series_axes = PyObject[]
    right_x, right_width, right_height = 0.64, 0.345, 0.205
    right_y = [0.765, 0.525, 0.285, 0.045]
    obs_dates = DateTime.(EVAL_TIMES)
    series_handles = nothing
    series_labels = nothing
    for (panel_index, city_index) in enumerate(SERIES_CITY_INDICES)
        city = CITIES[city_index]
        axis = figure.add_axes([right_x, right_y[panel_index], right_width, right_height])
        push!(series_axes, axis)
        axis.bar(obs_dates, obs_city[:, city_index], width=0.02, color=colors["Observation"],
            edgecolor="none", label="Observation", zorder=1)
        for (model, marker) in (("AI-GAMFS", "o"), ("Aurora", "s"), ("RUF", nothing))
            times = plot_times[(model, city)]
            values = plot_values[(model, city)]
            order = sortperm(times)
            kwargs = Dict{Symbol, Any}(:color => colors[model], :linewidth => 0.9, :label => model, :zorder => 3)
            if marker !== nothing
                kwargs[:marker] = marker
                kwargs[:markersize] = 1.7
                kwargs[:markeredgewidth] = 0.0
            end
            axis.plot(times[order], values[order]; kwargs...)
        end
        letter = Char('c' + panel_index - 1)
        axis.text(0.015, 0.965, "($letter) $(CITY_ENGLISH[city_index])", transform=axis.transAxes,
            va="top", ha="left", fontproperties=font_panel_label)
        axis.text(0.53, 0.92, "RMSE:", transform=axis.transAxes, va="center", ha="left",
            color="black", fontproperties=font_metric, zorder=10)
        axis.text(0.53, 0.83, "R:", transform=axis.transAxes, va="center", ha="left",
            color="black", fontproperties=font_metric, zorder=10)
        value_x = (0.68, 0.79, 0.90)
        for (model_index, model) in enumerate(("AI-GAMFS", "Aurora", "RUF"))
            metrics = stitched_metrics[(model, city)]
            axis.text(value_x[model_index], 0.92, string(round(metrics.rmse; digits=1)),
                transform=axis.transAxes, va="center", ha="left", color=colors[model],
                fontproperties=font_metric, zorder=10)
            axis.text(value_x[model_index], 0.83, string(round(metrics.corr; digits=2)),
                transform=axis.transAxes, va="center", ha="left", color=colors[model],
                fontproperties=font_metric, zorder=10)
        end
        axis.set_xlim(time_key(EVAL_FIRST), time_key(EVAL_LAST))
        axis.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        axis.grid(true, color="#D9D9D9", linewidth=0.35, linestyle="--", alpha=0.75)
        axis.tick_params(axis="both", labelsize=6.2, width=0.5, length=2.2, pad=1.4)
        for label in vcat(axis.get_xticklabels(), axis.get_yticklabels())
            label.set_fontproperties(font_map_tick)
        end
        panel_index < length(SERIES_CITY_INDICES) && axis.tick_params(labelbottom=false)
        axis.set_ylabel(raw"PM$_{10}$ ($\mu$g m$^{-3}$)", fontproperties=font_series_ylabel, labelpad=1.0)
        if panel_index == 1
            series_handles, series_labels = axis.get_legend_handles_labels()
        end
    end
    original_labels = convert(Vector{String}, series_labels)
    legend_order = [findfirst(==(label), original_labels) for label in ("Observation", "AI-GAMFS", "Aurora", "RUF")]
    any(isnothing, legend_order) && error("Time-series legend entries are incomplete")
    ordered_handles = [series_handles[index] for index in Int.(legend_order)]
    ordered_labels = [original_labels[index] for index in Int.(legend_order)]
    figure.legend(ordered_handles, ordered_labels, loc="lower center", bbox_to_anchor=(0.835, 0.005),
        ncol=4, frameon=false, prop=font_small, handlelength=1.8, columnspacing=0.8)

    png_path = joinpath(OUTPUT_DIR, "$OUTPUT_BASE.png")
    pdf_path = joinpath(OUTPUT_DIR, "$OUTPUT_BASE.pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=0.025)
    figure.savefig(pdf_path, dpi=600, bbox_inches="tight", pad_inches=0.025)
    plt.close(figure)
    println("Wrote: $png_path")
    println("Wrote: $pdf_path")
end

main()
