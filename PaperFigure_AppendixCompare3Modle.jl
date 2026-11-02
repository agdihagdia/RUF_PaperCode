#!/usr/bin/env julia

push!(LOAD_PATH, "/home/zsm/codeLib")

using CSV
using DataFrames
using Dates
using NaNStatistics
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
const OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const RUF_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const AURORA_ROOT = "/home/zsm/IETM/Aurora"
const AI_ROOT = "/home/zsm/IETM/AI-GAMFS/inference"
const OUTPUT_DIR = "/home/zsm/IETM/figure/PapreFigure/20260826"
const OUTPUT_BASE = "PaperFigure_AppendixCompare3Modle"
const STAMP_FORMAT = dateformat"yyyymmddHH"

clean_text(value) = value === missing ? "" : value isa AbstractVector ? join(string.(value)) : strip(string(value))
valid_obs(value) = isfinite(value) && 0.0 <= value < 2000.0
valid_forecast(value) = isfinite(value) && value >= 0.0

function finite_mean(values, validator)
    selected = Float64[value for value in values if validator(value)]
    return isempty(selected) ? NaN : mean(selected)
end

time_key(value) = DateTime(astimezone(value, BJT))
time_iso(value) = Dates.format(time_key(value), dateformat"yyyy-mm-ddTHH:MM:SS")

function empty_statistics()
    return (
        rmse=NaN, mae=NaN, corr=NaN, bias=NaN, pstd=NaN,
        rstd=NaN, crmsd=NaN, nmb=NaN, nmae=NaN,
    )
end

function main()
    mkpath(OUTPUT_DIR)
    println("Evaluation window (BJT): $(time_iso(EVAL_FIRST)) -- $(time_iso(EVAL_LAST))")

    # The observation file is the authoritative source of national-control station
    # IDs and coordinates. MATLAB character arrays are normalized with clean_text.
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
    unique_ids = unique(station_ids)
    length(unique_ids) == length(station_ids) || error("Duplicated national-control station IDs in observation metadata")
    city_positions = [findall(==(city), station_cities) for city in CITIES]
    any(isempty, city_positions) && error("At least one requested city has no national-control station")
    println("National-control station counts: ", join(["$(CITIES[i])=$(length(city_positions[i]))" for i in eachindex(CITIES)], ", "))

    pm10_position = findfirst(==("PM10"), uppercase.(clean_text.(observations.varNames)))
    pm10_position === nothing && error("PM10 is absent from $OBS_FILE")
    obs_time_lookup = Dict(time_key(t) => i for (i, t) in enumerate(observations.times))
    obs_indices = [get(obs_time_lookup, time_key(t), 0) for t in EVAL_TIMES]
    any(==(0), obs_indices) && error("Observation hours are missing inside the requested evaluation window")
    obs_field = Float64.(observations.staDataAll[obs_indices, pm10_position, selected_positions, 1])
    obs_field[.!valid_obs.(obs_field)] .= NaN
    obs_city = fill(NaN, length(EVAL_TIMES), length(CITIES))
    for city_index in eachindex(CITIES), time_index in eachindex(EVAL_TIMES)
        obs_city[time_index, city_index] = finite_mean(@view(obs_field[time_index, city_positions[city_index]]), valid_obs)
    end
    obs_row_lookup = Dict(time_key(t) => i for (i, t) in enumerate(EVAL_TIMES))

    pair_store = Dict{Tuple{String, String, Int}, Tuple{Vector{Float64}, Vector{Float64}}}()
    plot_times = Dict((model, city) => DateTime[] for model in ("AI-GAMFS", "Aurora", "RUF") for city in CITIES)
    plot_values = Dict((model, city) => Float64[] for model in ("AI-GAMFS", "Aurora", "RUF") for city in CITIES)

    function consume_forecast!(model::String, valid_time::DateTime, lead::Int, values::AbstractVector{<:Real}; plot_point::Bool=false)
        haskey(obs_row_lookup, valid_time) || return
        length(values) == length(station_ids) || error("$model station-vector length changed at $valid_time")
        obs_index = obs_row_lookup[valid_time]
        for city_index in eachindex(CITIES)
            city = CITIES[city_index]
            forecast_city = finite_mean(@view(values[city_positions[city_index]]), valid_forecast)
            reference_city = obs_city[obs_index, city_index]
            if valid_obs(reference_city) && valid_forecast(forecast_city)
                reference, predict = get!(pair_store, (model, city, lead), (Float64[], Float64[]))
                push!(reference, reference_city)
                push!(predict, forecast_city)
            end
            if plot_point
                push!(plot_times[(model, city)], valid_time)
                push!(plot_values[(model, city)], forecast_city)
            end
        end
    end

    # RUF: one cycle every hour; lead 0 gives the freshest rolling time series.
    println("Reading hourly RUF cycles and all native lead hours")
    ruf_cycle_first = max(EVAL_FIRST - Hour(72), ZonedDateTime(DateTime(2025, 4, 8, 0), BJT))
    ruf_cycles = collect(ruf_cycle_first:Hour(1):EVAL_LAST)
    template_stamp = Dates.format(first(ruf_cycles), STAMP_FORMAT)
    template_path = joinpath(RUF_ROOT, template_stamp, "staData.$template_stamp.mat")
    isfile(template_path) || error("Missing RUF template: $template_path")
    ruf_template = staDATA(template_path)
    ruf_ids = clean_text.(ruf_template.staInfo.staID)
    ruf_lookup = Dict(id => i for (i, id) in enumerate(ruf_ids))
    all(haskey(ruf_lookup, id) for id in station_ids) || error("RUF staData is missing one or more selected national-control stations")
    ruf_positions = [ruf_lookup[id] for id in station_ids]
    for (cycle_index, cycle) in enumerate(ruf_cycles)
        stamp = Dates.format(cycle, STAMP_FORMAT)
        path = joinpath(RUF_ROOT, stamp, "staData.$stamp.mat")
        isfile(path) || error("Missing RUF cycle required by the shortened analysis: $path")
        object = staDATA(path)
        pm_position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
        pm_position === nothing && error("PM10 absent in $path")
        lead_lookup = Dict(Int(lead) => i for (i, lead) in enumerate(object.leadTimes))
        for lead in 0:72
            valid_time = time_key(cycle + Hour(lead))
            time_key(EVAL_FIRST) <= valid_time <= time_key(EVAL_LAST) || continue
            haskey(lead_lookup, lead) || error("RUF lead $lead absent in $path")
            field = Float64.(object.staDataAll[1, pm_position, ruf_positions, lead_lookup[lead]])
            consume_forecast!("RUF", valid_time, lead, field; plot_point=(lead == 0))
        end
        (cycle_index == 1 || cycle_index % 24 == 0 || cycle_index == length(ruf_cycles)) &&
            println("  RUF cycle $cycle_index/$(length(ruf_cycles)): $stamp")
    end

    # Reuse the tested NetCDF reader and bilinear interpolation implementation in
    # ThreeModelEva. Python only extracts station vectors; all city averaging and
    # accuracy statistics remain in this Julia script.
    py"""
import sys, re
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
sys.path.insert(0, "/home/zsm/IETM/code/ThreeModelEva")
from three_model_eval_netcdf_plot import NetCDF, BilinearSampler, AI_COMPONENTS

def paper_extract_aurora(root, first_iso, last_iso, station_lats, station_lons):
    first = datetime.fromisoformat(first_iso)
    last = datetime.fromisoformat(last_iso)
    lats = np.asarray(station_lats, dtype=np.float64)
    lons = np.asarray(station_lons, dtype=np.float64)
    rows_t, rows_l, rows_v = [], [], []
    sampler = None
    for path in sorted(Path(root).glob("aurora_*T12.nc")):
        match = re.fullmatch(r"aurora_(\d{8})T12\.nc", path.name)
        if not match:
            continue
        init_bjt = datetime.strptime(match.group(1) + "12", "%Y%m%d%H") + timedelta(hours=8)
        selected = [lead for lead in range(0, 73, 12) if first <= init_bjt + timedelta(hours=lead) <= last]
        if not selected:
            continue
        with NetCDF(path) as ds:
            if sampler is None:
                sampler = BilinearSampler(ds.read_double("lat"), ds.read_double("lon"), lats, lons)
            lookup = {int(v): i for i, v in enumerate(ds.read_int("lead_time").tolist())}
            pm10 = ds.read_float("PM10")[0]
            for lead in selected:
                rows_t.append((init_bjt + timedelta(hours=lead)).isoformat())
                rows_l.append(lead)
                rows_v.append(sampler.sample(pm10[lookup[lead]]) * 1.0e9)
        print("Aurora", path.name, "done", flush=True)
    return {"times": rows_t, "leads": np.asarray(rows_l, dtype=np.int64), "values": np.vstack(rows_v)}

def paper_extract_ai(root, first_iso, last_iso, station_lats, station_lons):
    first = datetime.fromisoformat(first_iso)
    last = datetime.fromisoformat(last_iso)
    lats = np.asarray(station_lats, dtype=np.float64)
    lons = np.asarray(station_lons, dtype=np.float64)
    rows_t, rows_l, rows_v = [], [], []
    sampler = None
    for directory in sorted(Path(root).glob("????????_????")):
        try:
            init_label_utc = datetime.strptime(directory.name, "%Y%m%d_%H%M")
        except ValueError:
            continue
        # Empty placeholder directories exist before the first archived AI-GAMFS
        # cycle. Follow ThreeModelEva's established rule and skip unavailable
        # cycle/lead pairs; valid sample counts are retained in the statistics CSV.
        if not any(directory.glob("*.nc")):
            continue
        init_bjt = init_label_utc - timedelta(minutes=30) + timedelta(hours=8)
        selected = [lead for lead in range(3, 73, 3) if first <= init_bjt + timedelta(hours=lead) <= last]
        if not selected:
            continue
        for lead in selected:
            valid_label_utc = init_label_utc + timedelta(hours=lead)
            filename = f"AI_GAMFS.{init_label_utc:%Y%m%d_%H%M}+{valid_label_utc:%Y%m%d_%H%M}.V01.nc"
            path = directory / filename
            if not path.is_file():
                continue
            with NetCDF(path) as ds:
                if sampler is None:
                    sampler = BilinearSampler(ds.read_double("latitude"), ds.read_double("longitude"), lats, lons)
                forecast = np.zeros(len(lats), dtype=np.float64)
                for component in AI_COMPONENTS:
                    forecast += sampler.sample(ds.read_float(component)[0])
                forecast *= 1.0e9
            rows_t.append((init_bjt + timedelta(hours=lead)).isoformat())
            rows_l.append(lead)
            rows_v.append(forecast)
        print("AI-GAMFS", directory.name, "done", flush=True)
    return {"times": rows_t, "leads": np.asarray(rows_l, dtype=np.int64), "values": np.vstack(rows_v)}
"""

    function consume_python_rows!(model::String, result; plot_rule)
        times = convert(Vector{String}, result["times"])
        leads = Int.(vec(Array(result["leads"])))
        values = Float64.(Array(result["values"]))
        size(values, 1) == length(times) == length(leads) || error("$model extraction arrays are inconsistent")
        for row in eachindex(times)
            valid_time = DateTime(times[row], dateformat"yyyy-mm-ddTHH:MM:SS")
            consume_forecast!(model, valid_time, leads[row], @view(values[row, :]); plot_point=plot_rule(valid_time, leads[row]))
        end
    end

    first_iso = time_iso(EVAL_FIRST)
    last_iso = time_iso(EVAL_LAST)
    println("Reading Aurora daily cycles (12-hour native output)")
    aurora_result = py"paper_extract_aurora"(AURORA_ROOT, first_iso, last_iso, station_lats, station_lons)
    consume_python_rows!("Aurora", aurora_result; plot_rule=(time, lead) -> lead in (0, 12))

    println("Reading AI-GAMFS three-hour cycles")
    ai_result = py"paper_extract_ai"(AI_ROOT, first_iso, last_iso, station_lats, station_lons)
    consume_python_rows!("AI-GAMFS", ai_result; plot_rule=(time, lead) -> lead == 3)

    # Statistics shown inside each panel use the same freshest-update stitched
    # points that are drawn in the figure, paired with city-mean observations at
    # the corresponding valid times.
    stitched_metrics = Dict{Tuple{String, String}, Any}()
    for (city_index, city) in enumerate(CITIES), model in ("AI-GAMFS", "Aurora", "RUF")
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
        println("  Stitched metrics $city $model: n=$(length(reference)), RMSE=$(stitched_metrics[(model, city)].rmse), R=$(stitched_metrics[(model, city)].corr)")
    end

    # Each row is based on city-mean observation and city-mean forecast pairs.
    println("Calculating city-by-lead statistics with utilities.getStatistics")
    model_leads = Dict("RUF" => collect(0:72), "AI-GAMFS" => collect(3:3:72), "Aurora" => collect(0:12:72))
    model_update = Dict("RUF" => 1, "AI-GAMFS" => 3, "Aurora" => 24)
    native_output = Dict("RUF" => 1, "AI-GAMFS" => 3, "Aurora" => 12)
    statistics_table = DataFrame(
        city=String[], city_english=String[], national_station_count=Int[], model=String[],
        update_interval_hour=Int[], native_output_interval_hour=Int[], lead_hour=Int[], valid_city_pairs=Int[],
        evaluation_start_bjt=String[], evaluation_end_bjt=String[], rmse_ug_m3=Float64[], mae_ug_m3=Float64[],
        correlation_r=Float64[], bias_ug_m3=Float64[], forecast_std_ug_m3=Float64[], observation_std_ug_m3=Float64[],
        centered_rmse_ug_m3=Float64[], normalized_mean_bias=Float64[], normalized_mean_absolute_error=Float64[],
    )
    for (city_index, city) in enumerate(CITIES), model in ("AI-GAMFS", "Aurora", "RUF"), lead in model_leads[model]
        reference, predict = get(pair_store, (model, city, lead), (Float64[], Float64[]))
        metrics = isempty(reference) ? empty_statistics() : getStatistics(reference, predict)
        push!(statistics_table, (
            city, CITY_ENGLISH[city_index], length(city_positions[city_index]), model,
            model_update[model], native_output[model], lead, length(reference), first_iso, last_iso,
            metrics.rmse, metrics.mae, metrics.corr, metrics.bias, metrics.pstd, metrics.rstd,
            metrics.crmsd, metrics.nmb, metrics.nmae,
        ))
    end
    statistics_path = joinpath(OUTPUT_DIR, "$(OUTPUT_BASE)_statistics.csv")
    CSV.write(statistics_path, statistics_table)

    println("Drawing the 2 x 3 city-average comparison figure")
    matplotlib = pyimport("matplotlib")
    configure_matplotlib!(matplotlib)
    plt = pyimport("matplotlib.pyplot")
    mdates = pyimport("matplotlib.dates")
    fm = matplotlib.font_manager
    font_regular_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_bold_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    figure_fonts = font_properties(fm, font_regular_path, font_bold_path)
    font_regular = figure_fonts.tick
    font_bold = figure_fonts.panel_title
    font_metric = figure_fonts.annotation
    font_axis = figure_fonts.axis_label
    font_figure_legend = figure_fonts.figure_legend

    colors = Dict("Observation" => OBS_COLOR, "AI-GAMFS" => AI_COLOR,
        "Aurora" => AURORA_COLOR, "RUF" => RUF_COLOR)
    line_styles = Dict("AI-GAMFS" => "--", "Aurora" => "-.", "RUF" => "-")
    figure, axes = plt.subplots(2, 3, figsize=(FIGURE_WIDTH, 4.8), dpi=600, sharex=true)
    axes_flat = vec(permutedims(Array(axes)))
    obs_dates = DateTime.(EVAL_TIMES)
    for city_index in eachindex(CITIES)
        city = CITIES[city_index]
        axis = axes_flat[city_index]
        axis.bar(obs_dates, obs_city[:, city_index], width=0.02, color=colors["Observation"],
            edgecolor="none", label="Observation", zorder=1)
        for (model, marker) in (("AI-GAMFS", "o"), ("Aurora", "s"), ("RUF", nothing))
            times = plot_times[(model, city)]
            values = plot_values[(model, city)]
            order = sortperm(times)
            kwargs = Dict{Symbol, Any}(:color => colors[model], :linewidth => DATA_LINEWIDTH,
                :linestyle => line_styles[model], :label => model, :zorder => 3)
            if marker !== nothing
                kwargs[:marker] = marker
                kwargs[:markersize] = 1.7
                kwargs[:markeredgewidth] = 0.0
            end
            axis.plot(times[order], values[order]; kwargs...)
        end
        letter = Char('a' + city_index - 1)
        axis.text(0.02, 0.96, "($letter) $(CITY_ENGLISH[city_index])", transform=axis.transAxes,
            va="top", ha="left", fontproperties=font_bold)
        # Panels (a)-(e): upper-right blank area. Panel (f): upper-left blank
        # area below the city label. Values follow the legend order and use the
        # corresponding model colors.
        if city_index < 6
            label_x, value_x = 0.58, (0.62, 0.80, 0.98)
            rmse_label_y, rmse_value_y, corr_label_y, corr_value_y = 0.94, 0.865, 0.79, 0.715
        else
            label_x, value_x = 0.02, (0.16, 0.33, 0.50)
            rmse_label_y, rmse_value_y, corr_label_y, corr_value_y = 0.82, 0.745, 0.67, 0.595
        end
        axis.text(label_x, rmse_label_y, "RMSE:", transform=axis.transAxes,
            va="center", ha="left", color="black", fontproperties=font_metric, zorder=10)
        axis.text(label_x, corr_label_y, "R:", transform=axis.transAxes,
            va="center", ha="left", color="black", fontproperties=font_metric, zorder=10)
        for (model_index, model) in enumerate(("AI-GAMFS", "Aurora", "RUF"))
            metrics = stitched_metrics[(model, city)]
            axis.text(value_x[model_index], rmse_value_y, string(round(metrics.rmse; digits=1)),
                transform=axis.transAxes, va="center", ha="center", color=colors[model],
                fontproperties=font_metric, zorder=10)
            axis.text(value_x[model_index], corr_value_y, string(round(metrics.corr; digits=2)),
                transform=axis.transAxes, va="center", ha="center", color=colors[model],
                fontproperties=font_metric, zorder=10)
        end
        axis.set_xlim(time_key(EVAL_FIRST), time_key(EVAL_LAST))
        axis.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        axis.grid(true, color="#D9D9D9", linewidth=GRID_LINEWIDTH, linestyle="--", alpha=0.75)
        axis.tick_params(axis="both", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH, length=2.8, pad=1.5)
        for label in axis.get_xticklabels()
            label.set_fontproperties(font_regular)
            label.set_rotation(0)
        end
        for label in axis.get_yticklabels()
            label.set_fontproperties(font_regular)
        end
        city_index in (1, 4) && axis.set_ylabel(raw"PM$_{10}$ ($\mu$g m$^{-3}$)", fontproperties=font_axis)
    end
    handles, labels = axes_flat[1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.005), ncol=4,
        frameon=false, prop=font_figure_legend, handlelength=2.5, columnspacing=1.4)
    figure.subplots_adjust(left=0.09, right=0.99, bottom=0.15, top=0.97, wspace=0.20, hspace=0.18)
    png_path = joinpath(OUTPUT_DIR, "$(OUTPUT_BASE).png")
    pdf_path = joinpath(OUTPUT_DIR, "$(OUTPUT_BASE).pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    figure.savefig(pdf_path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(figure)

    println("Wrote: $png_path")
    println("Wrote: $pdf_path")
    println("Wrote: $statistics_path")
end

main()
