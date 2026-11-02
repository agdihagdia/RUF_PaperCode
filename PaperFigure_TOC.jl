#!/usr/bin/env julia

push!(LOAD_PATH, "/home/zsm/codeLib")

using CSV
using DataFrames
using Dates
using PyCall
using StaDataM
using Statistics
using TimeZones

include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle

const BJT = tz"Asia/Shanghai"
const OBSERVATION_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const BASELINE_ROOT = "/media/dataLL/zsm/paper_data/baseFolder"
const RUF3D_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const OUTPUT_DIR = "/home/zsm/IETM/figure/PapreFigure/20260903"
const OUTPUT_BASE = "PaperFigure_TOC"
const FORECAST_LEADS = collect(0:72)
const STAMP_FORMAT = dateformat"yyyymmddHH"

clean_text(value) = value === missing ? "" :
    value isa AbstractVector ? join(string.(value)) : strip(string(value))
time_key(value) = DateTime(astimezone(value, BJT))
valid_observation(value) = isfinite(value) && 0.0 <= value < 2000.0
valid_forecast(value) = isfinite(value) && value >= 0.0

function discover_cycles(root, model)
    cycles = NamedTuple{(:stamp, :initialization, :path), Tuple{String, DateTime, String}}[]
    for directory in sort(readdir(root; join=true))
        isdir(directory) || continue
        stamp = basename(directory)
        occursin(r"^\d{10}$", stamp) || continue
        path = joinpath(directory, "staData.$stamp.mat")
        isfile(path) || continue
        initialization = DateTime(stamp, STAMP_FORMAT)
        model == "Baseline" && (initialization += Hour(8))
        push!(cycles, (stamp=stamp, initialization=initialization, path=path))
    end
    isempty(cycles) && error("No forecast cycles found for $model in $root")
    return cycles
end

function pm10_position(object, label)
    position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
    position === nothing && error("PM10 is absent from $label")
    return position
end

function station_positions(object_ids, station_ids, label)
    length(unique(object_ids)) == length(object_ids) || error("Duplicate station IDs in $label")
    lookup = Dict(id => index for (index, id) in enumerate(object_ids))
    missing_ids = [id for id in station_ids if !haskey(lookup, id)]
    isempty(missing_ids) || error("$label misses $(length(missing_ids)) required stations")
    return [lookup[id] for id in station_ids]
end

mutable struct MetricAccumulator
    paired_values::Int
    sum_error2::Float64
    sum_observation::Float64
    sum_forecast::Float64
    sum_observation2::Float64
    sum_forecast2::Float64
    sum_cross::Float64
    cycle_count::Int
    station_seen::BitVector
end

MetricAccumulator(station_count) = MetricAccumulator(
    0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, falses(station_count),
)

function update_metrics!(accumulator, reference, prediction)
    paired = valid_observation.(reference) .& valid_forecast.(prediction)
    indices = findall(paired)
    isempty(indices) && return
    observed = Float64.(reference[indices])
    forecast = Float64.(prediction[indices])
    errors = forecast .- observed
    accumulator.paired_values += length(indices)
    accumulator.sum_error2 += sum(errors .^ 2)
    accumulator.sum_observation += sum(observed)
    accumulator.sum_forecast += sum(forecast)
    accumulator.sum_observation2 += sum(observed .^ 2)
    accumulator.sum_forecast2 += sum(forecast .^ 2)
    accumulator.sum_cross += sum(observed .* forecast)
    accumulator.cycle_count += 1
    accumulator.station_seen[indices] .= true
end

function finalize_metrics(accumulator)
    n = accumulator.paired_values
    n >= 3 || error("Fewer than three paired station-cycle values are available")
    numerator = n * accumulator.sum_cross - accumulator.sum_observation * accumulator.sum_forecast
    denominator2 =
        (n * accumulator.sum_observation2 - accumulator.sum_observation^2) *
        (n * accumulator.sum_forecast2 - accumulator.sum_forecast^2)
    denominator2 > 0.0 || error("Pearson correlation is undefined")
    return (
        rmse=sqrt(accumulator.sum_error2 / n),
        correlation=numerator / sqrt(denominator2),
        paired_values=n,
        cycle_count=accumulator.cycle_count,
        station_count=count(accumulator.station_seen),
    )
end

function main()
    mkpath(OUTPUT_DIR)
    isfile(OBSERVATION_FILE) || error("Missing observation file: $OBSERVATION_FILE")
    baseline_cycles = discover_cycles(BASELINE_ROOT, "Baseline")
    ruf3d_cycles = discover_cycles(RUF3D_ROOT, "RUF3d")
    println("Discovered Baseline cycles: $(length(baseline_cycles)), " *
        "$(first(baseline_cycles).stamp)-$(last(baseline_cycles).stamp)")
    println("Discovered RUF3d cycles: $(length(ruf3d_cycles)), " *
        "$(first(ruf3d_cycles).stamp)-$(last(ruf3d_cycles).stamp)")

    println("Loading observations and template forecast metadata")
    observations = staDATA(OBSERVATION_FILE)
    baseline_template = staDATA(first(baseline_cycles).path)
    ruf3d_template = staDATA(first(ruf3d_cycles).path)
    baseline_ids = clean_text.(baseline_template.staInfo.staID)
    baseline_types = clean_text.(baseline_template.staInfo.type)
    all(==("B"), baseline_types) || error("Baseline station set is not exclusively national-control stations")
    observation_ids = Set(clean_text.(observations.staInfo.staID))
    ruf3d_ids = Set(clean_text.(ruf3d_template.staInfo.staID))
    station_ids = [id for id in baseline_ids if id in observation_ids && id in ruf3d_ids]
    isempty(station_ids) && error("No common national-control stations")
    println("Common national-control stations: $(length(station_ids))/$(length(baseline_ids)) baseline stations")

    observation_object_ids = clean_text.(observations.staInfo.staID)
    ruf3d_template_ids = clean_text.(ruf3d_template.staInfo.staID)
    observation_positions = station_positions(observation_object_ids, station_ids, "observations")
    baseline_positions = station_positions(baseline_ids, station_ids, "Baseline template")
    ruf3d_positions = station_positions(ruf3d_template_ids, station_ids, "RUF3d template")
    observation_pm10 = pm10_position(observations, "observations")
    observation_time_lookup = Dict(time_key(time) => index for (index, time) in enumerate(observations.times))
    observation_cache = Dict{DateTime, Vector{Float64}}()

    function observation_values(valid_time)
        return get!(observation_cache, valid_time) do
            time_position = get(observation_time_lookup, valid_time, 0)
            time_position > 0 || error("Observation is missing at $valid_time BJT")
            Float64.(observations.staDataAll[
                time_position, observation_pm10, observation_positions, 1,
            ])
        end
    end

    accumulators = Dict(
        "Baseline" => [MetricAccumulator(length(station_ids)) for _ in FORECAST_LEADS],
        "RUF3d" => [MetricAccumulator(length(station_ids)) for _ in FORECAST_LEADS],
    )

    function consume_cycles!(model, cycles, template_ids, positions)
        println("Reading all $model cycles")
        for (cycle_index, cycle) in enumerate(cycles)
            object = staDATA(cycle.path)
            object_ids = clean_text.(object.staInfo.staID)
            object_ids == template_ids || error("$model station ordering differs in $(cycle.path)")
            time_key(first(object.times)) == cycle.initialization ||
                error("$model initialization metadata differs in $(cycle.path)")
            pm10_index = pm10_position(object, cycle.path)
            lead_lookup = Dict(Int(lead) => index for (index, lead) in enumerate(object.leadTimes))
            for lead in FORECAST_LEADS
                lead_position = get(lead_lookup, lead, 0)
                lead_position > 0 || error("$model lead $lead h is missing in $(cycle.path)")
                valid_time = cycle.initialization + Hour(lead)
                reference = observation_values(valid_time)
                prediction = Float64.(object.staDataAll[1, pm10_index, positions, lead_position])
                update_metrics!(accumulators[model][lead + 1], reference, prediction)
            end
            (cycle_index == 1 || cycle_index % 24 == 0 || cycle_index == length(cycles)) &&
                println("  $model cycle $cycle_index/$(length(cycles)): $(cycle.stamp)")
        end
    end

    consume_cycles!("Baseline", baseline_cycles, baseline_ids, baseline_positions)
    baseline_template = nothing
    GC.gc()
    consume_cycles!("RUF3d", ruf3d_cycles, ruf3d_template_ids, ruf3d_positions)

    results = Dict(
        model => [finalize_metrics(accumulator) for accumulator in model_accumulators]
        for (model, model_accumulators) in accumulators
    )
    statistics_table = DataFrame(
        model=String[], lead_hour=Int[], rmse_ug_m3=Float64[], correlation_r=Float64[],
        initialization_count=Int[], paired_values=Int[], valid_station_count=Int[],
        first_initialization_bjt=String[], last_initialization_bjt=String[],
        station_scope=String[], aggregation=String[],
    )
    for (model, cycles) in (("Baseline", baseline_cycles), ("RUF3d", ruf3d_cycles))
        for lead in FORECAST_LEADS
            metric = results[model][lead + 1]
            push!(statistics_table, (
                model, lead, metric.rmse, metric.correlation, metric.cycle_count,
                metric.paired_values, metric.station_count,
                Dates.format(first(cycles).initialization, dateformat"yyyy-mm-ddTHH:MM:SS"),
                Dates.format(last(cycles).initialization, dateformat"yyyy-mm-ddTHH:MM:SS"),
                "common national-control stations",
                "pooled valid station-initialization pairs; model cycles evaluated independently",
            ))
        end
    end
    statistics_path = joinpath(OUTPUT_DIR, "$(OUTPUT_BASE)_statistics.csv")
    CSV.write(statistics_path, statistics_table)

    baseline_rmse = [metric.rmse for metric in results["Baseline"]]
    ruf3d_rmse = [metric.rmse for metric in results["RUF3d"]]
    baseline_correlation = [metric.correlation for metric in results["Baseline"]]
    ruf3d_correlation = [metric.correlation for metric in results["RUF3d"]]
    for (label, values) in (("Baseline", baseline_rmse), ("RUF3d", ruf3d_rmse))
        println("$label RMSE range: $(minimum(values))-$(maximum(values)) μg m^-3")
    end
    for (label, values) in (("Baseline", baseline_correlation), ("RUF3d", ruf3d_correlation))
        println("$label R range: $(minimum(values))-$(maximum(values))")
    end

    matplotlib = pyimport("matplotlib")
    configure_matplotlib!(matplotlib)
    plt = pyimport("matplotlib.pyplot")
    font_manager = matplotlib.font_manager
    regular_font = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    toc_font_increment = 2.0
    fonts = (
        axis_label=font_manager.FontProperties(
            fname=regular_font, size=AXIS_LABEL_SIZE + toc_font_increment,
        ),
        tick=font_manager.FontProperties(
            fname=regular_font, size=TICK_SIZE + toc_font_increment,
        ),
        figure_legend=font_manager.FontProperties(
            fname=regular_font, size=FIGURE_LEGEND_SIZE + toc_font_increment,
        ),
    )

    function draw_metric_figure(series, ylabel, output_base, legend_location; ylimits=nothing)
        figure, axis = plt.subplots(figsize=(3.6, 3.6), dpi=600)
        figure.subplots_adjust(left=0.18, right=0.97, bottom=0.16, top=0.96)
        axis.plot(FORECAST_LEADS, series[1], color="#4D4D4D", linestyle=":",
            linewidth=1.3, label="Baseline", zorder=3)
        axis.plot(FORECAST_LEADS, series[2], color=RUF_COLOR, linestyle="-",
            linewidth=1.4, label="RUF", zorder=4)
        axis.set_xlim(0, 72)
        ylimits === nothing ? axis.set_ylim(bottom=0) : axis.set_ylim(ylimits[1], ylimits[2])
        axis.set_xticks(0:12:72)
        axis.set_xlabel("Lead Time (h)", fontproperties=fonts.axis_label, labelpad=3.0)
        axis.set_ylabel(ylabel, fontproperties=fonts.axis_label, labelpad=3.0)
        axis.tick_params(axis="both", which="major", width=AXIS_LINEWIDTH, length=3.0, pad=2.0)
        for label in axis.get_xticklabels()
            label.set_fontproperties(fonts.tick)
        end
        for label in axis.get_yticklabels()
            label.set_fontproperties(fonts.tick)
        end
        axis.grid(true, color="#D9D9D9", linewidth=GRID_LINEWIDTH, linestyle="--", alpha=0.8)
        axis.set_axisbelow(true)
        legend = axis.legend(loc=legend_location, frameon=false, handlelength=2.6,
            handletextpad=0.6, borderaxespad=0.6, prop=fonts.figure_legend)
        legend.set_zorder(10)
        png_path = joinpath(OUTPUT_DIR, "$output_base.png")
        pdf_path = joinpath(OUTPUT_DIR, "$output_base.pdf")
        figure.savefig(png_path, dpi=600)
        figure.savefig(pdf_path)
        plt.close(figure)
        return png_path, pdf_path
    end

    rmse_png_path, rmse_pdf_path = draw_metric_figure(
        (baseline_rmse, ruf3d_rmse), "RMSE", OUTPUT_BASE, "lower right",
    )
    all_correlations = vcat(baseline_correlation, ruf3d_correlation)
    correlation_lower = max(-1.0, min(0.0, floor(10 * (minimum(all_correlations) - 0.05)) / 10))
    correlation_upper = min(1.0, max(0.5, ceil(10 * (maximum(all_correlations) + 0.05)) / 10))
    correlation_png_path, correlation_pdf_path = draw_metric_figure(
        (baseline_correlation, ruf3d_correlation), "R", "$(OUTPUT_BASE)_R", "upper right",
        ylimits=(correlation_lower, correlation_upper),
    )
    println("Wrote: $statistics_path")
    println("Wrote: $rmse_png_path")
    println("Wrote: $rmse_pdf_path")
    println("Wrote: $correlation_png_path")
    println("Wrote: $correlation_pdf_path")
end

main()
