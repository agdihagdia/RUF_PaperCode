#!/usr/bin/env julia

push!(LOAD_PATH, "/home/zsm/codeLib")

using CSV
using DataFrames
using Dates
using PyCall
using StaDataM
using TimeZones

include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle

const OBSERVATION_FILE = "/home/zsm/IETM/obs/qc202504.mat"
const FORECAST_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const OUTPUT_DIR = "/home/zsm/IETM/figure/PapreFigure/20260910"
const FORECAST_LEADS = collect(0:24)
const STATION_SCOPE = "national"
const BJT = tz"Asia/Shanghai"
const STAMP_FORMAT = dateformat"yyyymmddHH"
# The reference panel's outer bounds are approximately 212 by 115 pixels.
const REFERENCE_PANEL_RATIO = 212 / 115
const BAR_COLOR = "#BED8F3"

clean_text(value) = value === missing ? "" :
    value isa AbstractVector ? join(string.(value)) : strip(string(value))
time_key(value) = DateTime(astimezone(value, BJT))
valid_observation(value) = isfinite(value) && 0.0 <= value < 2000.0
valid_forecast(value) = isfinite(value) && value >= 0.0

function pm10_position(object, label)
    position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
    position === nothing && error("PM10 is absent from $label")
    return position
end

function discover_cycles()
    cycles = NamedTuple{(:stamp, :initialization, :path), Tuple{String, DateTime, String}}[]
    for directory in sort(readdir(FORECAST_ROOT; join=true))
        isdir(directory) || continue
        stamp = basename(directory)
        occursin(r"^\d{10}$", stamp) || continue
        path = joinpath(directory, "staData.$stamp.mat")
        isfile(path) || continue
        push!(cycles, (stamp=stamp, initialization=DateTime(stamp, STAMP_FORMAT), path=path))
    end
    isempty(cycles) && error("No complete RUF3d station forecast files found")
    return cycles
end

function draw_bars(plt, fm, statistics, basename, square, upper_limit)
    regular_font = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    axis_font = fm.FontProperties(fname=regular_font, size=AXIS_LABEL_SIZE)
    tick_font = fm.FontProperties(fname=regular_font, size=TICK_SIZE)
    width = 3.6
    height = square ? width : width / REFERENCE_PANEL_RATIO
    figure, axis = plt.subplots(figsize=(width, height), dpi=600)
    # Fixed canvas bounds preserve both requested output aspect ratios.
    figure.subplots_adjust(left=0.18, right=0.98, bottom=square ? 0.15 : 0.25, top=0.96)
    axis.bar(statistics.lead_hour, statistics.rmse_ug_m3,
        width=0.78, color=BAR_COLOR, edgecolor="none", linewidth=0)
    axis.set_xlim(first(FORECAST_LEADS) - 0.6, last(FORECAST_LEADS) + 0.6)
    axis.set_ylim(0, upper_limit)
    axis.set_xticks([first(FORECAST_LEADS), div(last(FORECAST_LEADS), 2), last(FORECAST_LEADS)])
    axis.set_xlabel("Lead time (h)", fontproperties=axis_font, labelpad=2.5)
    axis.set_ylabel(raw"RMSE (μg m$^{-3}$)", fontproperties=axis_font, labelpad=2.5)
    axis.tick_params(axis="both", which="major", width=AXIS_LINEWIDTH,
        length=2.5, pad=2.0, direction="out", colors="#333333")
    axis.spines["top"].set_visible(false)
    axis.spines["right"].set_visible(false)
    for side in ("left", "bottom")
        axis.spines[side].set_color("#555555")
    end
    for label in vcat(collect(axis.get_xticklabels()), collect(axis.get_yticklabels()))
        label.set_fontproperties(tick_font)
    end
    for extension in ("png", "pdf")
        path = joinpath(OUTPUT_DIR, "$basename.$extension")
        figure.savefig(path, dpi=600, facecolor="white")
        println("Wrote: $path")
    end
    plt.close(figure)
end

function main()
    cycles = discover_cycles()
    println("RUF3d cycles: $(length(cycles)), $(first(cycles).stamp)-$(last(cycles).stamp)")
    println("Loading specified observations: $OBSERVATION_FILE")
    flush(stdout)
    observations = staDATA(OBSERVATION_FILE)
    template = staDATA(first(cycles).path)
    observation_ids = clean_text.(observations.staInfo.staID)
    template_ids = clean_text.(template.staInfo.staID)
    length(unique(observation_ids)) == length(observation_ids) || error("Duplicate observation station IDs")
    length(unique(template_ids)) == length(template_ids) || error("Duplicate RUF3d station IDs")
    observation_lookup = Dict(id => index for (index, id) in enumerate(observation_ids))
    station_types = clean_text.(template.staInfo.type)
    STATION_SCOPE in ("national", "all") || error("Unsupported station scope")
    positions = [index for (index, id) in enumerate(template_ids)
        if haskey(observation_lookup, id) && (STATION_SCOPE == "all" || station_types[index] == "B")]
    isempty(positions) && error("No shared stations in the selected station scope")
    station_ids = template_ids[positions]
    observation_positions = [observation_lookup[id] for id in station_ids]
    println("Shared stations ($STATION_SCOPE): $(length(station_ids)); leads: $(first(FORECAST_LEADS))-$(last(FORECAST_LEADS)) h")
    flush(stdout)

    obs_pm10 = pm10_position(observations, OBSERVATION_FILE)
    observation_time_lookup = Dict(time_key(time) => index for (index, time) in enumerate(observations.times))
    length(observation_time_lookup) == length(observations.times) || error("Duplicate observation timestamps")
    observation_cache = Dict{DateTime, Vector{Float64}}()
    function observation_values(valid_time)
        return get!(observation_cache, valid_time) do
            index = get(observation_time_lookup, valid_time, 0)
            index > 0 || error("Missing observation timestamp: $valid_time BJT")
            Float64.(observations.staDataAll[index, obs_pm10, observation_positions, 1])
        end
    end

    squared_errors = zeros(Float64, length(FORECAST_LEADS))
    paired_counts = zeros(Int, length(FORECAST_LEADS))
    cycle_counts = zeros(Int, length(FORECAST_LEADS))
    station_seen = [falses(length(station_ids)) for _ in FORECAST_LEADS]
    cycle_statistics = DataFrame(initialization_bjt=String[], lead_hour=Int[],
        paired_values=Int[], sum_squared_error=Float64[])
    for (cycle_index, cycle) in enumerate(cycles)
        forecast = staDATA(cycle.path)
        clean_text.(forecast.staInfo.staID) == template_ids || error("Station ordering changed: $(cycle.path)")
        length(forecast.times) == 1 || error("Expected a single initialization: $(cycle.path)")
        time_key(first(forecast.times)) == cycle.initialization || error("Initialization metadata mismatch: $(cycle.path)")
        model_pm10 = pm10_position(forecast, cycle.path)
        lead_lookup = Dict(Int(lead) => index for (index, lead) in enumerate(forecast.leadTimes))
        length(lead_lookup) == length(forecast.leadTimes) || error("Duplicate forecast leads")
        for (lead_index, lead) in enumerate(FORECAST_LEADS)
            model_lead = get(lead_lookup, lead, 0)
            model_lead > 0 || error("Missing lead $lead h: $(cycle.path)")
            reference = observation_values(cycle.initialization + Hour(lead))
            prediction = Float64.(forecast.staDataAll[1, model_pm10, positions, model_lead])
            paired = valid_observation.(reference) .& valid_forecast.(prediction)
            count_paired = count(paired)
            error_sum = sum(abs2, prediction[paired] .- reference[paired])
            squared_errors[lead_index] += error_sum
            paired_counts[lead_index] += count_paired
            cycle_counts[lead_index] += count_paired > 0
            station_seen[lead_index] .|= paired
            push!(cycle_statistics, (string(cycle.initialization), lead, count_paired, error_sum))
        end
        if cycle_index == 1 || cycle_index % 24 == 0 || cycle_index == length(cycles)
            println("Processed RUF3d cycle $cycle_index/$(length(cycles)): $(cycle.stamp)")
            flush(stdout)
        end
    end
    all(>(0), paired_counts) || error("No valid pairs at one or more leads")
    rmse = sqrt.(squared_errors ./ paired_counts)
    all(isfinite, rmse) || error("Non-finite RMSE")
    statistics = DataFrame(model=fill("RUF3d", length(FORECAST_LEADS)), lead_hour=FORECAST_LEADS,
        rmse_ug_m3=rmse, initialization_count=cycle_counts, paired_values=paired_counts,
        valid_station_count=count.(station_seen), eligible_station_count=fill(length(station_ids), length(FORECAST_LEADS)),
        first_initialization_bjt=fill(string(first(cycles).initialization), length(FORECAST_LEADS)),
        last_initialization_bjt=fill(string(last(cycles).initialization), length(FORECAST_LEADS)),
        station_scope=fill(STATION_SCOPE, length(FORECAST_LEADS)),
        aggregation=fill("pooled valid station-initialization pairs", length(FORECAST_LEADS)))
    mkpath(OUTPUT_DIR)
    CSV.write(joinpath(OUTPUT_DIR, "PaperFigure_MethodFig_statistics.csv"), statistics)
    CSV.write(joinpath(OUTPUT_DIR, "PaperFigure_MethodFig_cycle_statistics.csv"), cycle_statistics)
    CSV.write(joinpath(OUTPUT_DIR, "PaperFigure_MethodFig_stations.csv"),
        DataFrame(station_id=station_ids, station_type=station_types[positions]))
    println("RMSE range: $(minimum(rmse))-$(maximum(rmse)) μg m^-3")

    matplotlib = pyimport("matplotlib")
    configure_matplotlib!(matplotlib)
    plt = pyimport("matplotlib.pyplot")
    # Both variants use identical bars and axis limits.
    upper_limit = max(20.0, ceil(maximum(rmse) * 1.10 / 20) * 20)
    draw_bars(plt, matplotlib.font_manager, statistics, "PaperFigure_MethodFig", false, upper_limit)
    draw_bars(plt, matplotlib.font_manager, statistics, "PaperFigure_MethodFig_square", true, upper_limit)
end

main()
