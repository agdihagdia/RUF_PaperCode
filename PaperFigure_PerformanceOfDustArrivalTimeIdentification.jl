#!/usr/bin/env julia

push!(LOAD_PATH, "/home/zsm/codeLib")

using Dates
using NaNStatistics
using PyCall
using StaDataM
using Statistics
using TimeZones

include(joinpath(@__DIR__, "FigureStyle.jl"))
using .FigureStyle

const BJT = tz"Asia/Shanghai"
const OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const RUF_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const OUTPUT_DIR = "/home/zsm/IETM/figure/PapreFigure/20260826"
const OUTPUT_BASE = "PaperFigure_PerformanceOfDustArrivalTimeIdentification"
const STAMP_FORMAT = dateformat"yyyymmddHH"
const FIRST_VALID = DateTime(2025, 4, 11, 0)
const LAST_VALID = DateTime(2025, 4, 14, 0)
const VALID_TIMES = collect(FIRST_VALID:Hour(1):LAST_VALID)
const LEAD_HOUR = 2
const DERIVATIVE_THRESHOLD = 0.1
const CONCENTRATION_THRESHOLD = 200.0

const CITY_PANEL_GROUPS = (
    top = ["太原", "呼和浩特", "北京", "长治", "晋城"],
    left = ["邵阳", "永州", "桂林"],
    right = ["洛阳", "南阳", "襄阳"],
    bottom = ["清远", "茂名", "海口", "常德", "宜昌"],
)
const CITY_ENGLISH = Dict(
    "太原" => "Taiyuan", "呼和浩特" => "Hohhot", "北京" => "Beijing",
    "长治" => "Changzhi", "晋城" => "Jincheng", "邵阳" => "Shaoyang",
    "永州" => "Yongzhou", "桂林" => "Guilin", "洛阳" => "Luoyang",
    "南阳" => "Nanyang", "襄阳" => "Xiangyang", "清远" => "Qingyuan",
    "茂名" => "Maoming", "海口" => "Haikou", "常德" => "Changde",
    "宜昌" => "Yichang",
)
const CITIES = vcat(
    CITY_PANEL_GROUPS.top, CITY_PANEL_GROUPS.left,
    CITY_PANEL_GROUPS.right, CITY_PANEL_GROUPS.bottom,
)

clean_text(value) = value === missing ? "" : value isa AbstractVector ? join(string.(value)) : strip(string(value))
stamp(value::DateTime) = Dates.format(value, STAMP_FORMAT)
valid_obs(value) = isfinite(value) && 0.0 <= value < 2000.0
valid_forecast(value) = isfinite(value) && value >= 0.0

function finite_mean(values, validator)
    selected = Float64[value for value in values if validator(value)]
    return isempty(selected) ? NaN : mean(selected)
end

include(joinpath(@__DIR__, "Dependency", "change_point_detect.jl"))

function main()
    mkpath(OUTPUT_DIR)
    println("Valid-time window (BJT): $FIRST_VALID -- $LAST_VALID")
    println("RUF fixed lead: $(LEAD_HOUR) h")

    template_cycle = FIRST_VALID - Hour(LEAD_HOUR)
    template_text = stamp(template_cycle)
    template_path = joinpath(RUF_ROOT, template_text, "staData.$template_text.mat")
    isfile(template_path) || error("Missing RUF template: $template_path")
    template = staDATA(template_path)
    template_ids = clean_text.(template.staInfo.staID)
    template_cities = clean_text.(template.staInfo.city)
    selected_positions = findall(in(CITIES), template_cities)
    station_ids = template_ids[selected_positions]
    station_cities = template_cities[selected_positions]
    station_lons = Float64.(template.staInfo.lon[selected_positions])
    station_lats = Float64.(template.staInfo.lat[selected_positions])
    city_groups = Dict(city => findall(==(city), station_cities) for city in CITIES)
    any(isempty(city_groups[city]) for city in CITIES) && error("At least one requested city has no RUF station")
    println("Selected city-station counts: ", join(["$city=$(length(city_groups[city]))" for city in CITIES], ", "))

    city_lon = Dict(city => mean(station_lons[city_groups[city]]) for city in CITIES)
    city_lat = Dict(city => mean(station_lats[city_groups[city]]) for city in CITIES)

    observations = staDATA(OBS_FILE)
    observation_ids = clean_text.(observations.staInfo.staID)
    observation_lookup = Dict(id => index for (index, id) in enumerate(observation_ids))
    all(haskey(observation_lookup, id) for id in station_ids) || error("Observation file is missing selected RUF stations")
    observation_positions = [observation_lookup[id] for id in station_ids]
    pm10_obs_position = findfirst(==("PM10"), uppercase.(clean_text.(observations.varNames)))
    pm10_obs_position === nothing && error("PM10 is absent from $OBS_FILE")
    observation_time_lookup = Dict(
        DateTime(astimezone(time, BJT)) => index for (index, time) in enumerate(observations.times)
    )
    observation_indices = [get(observation_time_lookup, time, 0) for time in VALID_TIMES]
    any(==(0), observation_indices) && error("Observation hours are missing inside the plotting window")
    observation_field = Float64.(
        observations.staDataAll[observation_indices, pm10_obs_position, observation_positions, 1]
    )

    observation_city = Dict(city => fill(NaN, length(VALID_TIMES)) for city in CITIES)
    for city in CITIES, time_index in eachindex(VALID_TIMES)
        observation_city[city][time_index] = finite_mean(
            @view(observation_field[time_index, city_groups[city]]), valid_obs
        )
    end

    ruf_city = Dict(city => fill(NaN, length(VALID_TIMES)) for city in CITIES)
    missing_cycles = String[]
    for (time_index, valid_time) in enumerate(VALID_TIMES)
        cycle = valid_time - Hour(LEAD_HOUR)
        cycle_text = stamp(cycle)
        path = joinpath(RUF_ROOT, cycle_text, "staData.$cycle_text.mat")
        if !isfile(path)
            push!(missing_cycles, cycle_text)
            println("MISSING|model=RUF|init_bjt=$cycle|lead=$LEAD_HOUR|valid_bjt=$valid_time|file=$path")
            continue
        end
        object = staDATA(path)
        clean_text.(object.staInfo.staID) == template_ids || error("RUF station ordering changed at $cycle_text")
        pm10_position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
        pm10_position === nothing && error("PM10 absent in $path")
        lead_position = findfirst(==(LEAD_HOUR), Int.(object.leadTimes))
        lead_position === nothing && error("RUF lead $LEAD_HOUR absent in $path")
        field = Float64.(object.staDataAll[1, pm10_position, selected_positions, lead_position])
        for city in CITIES
            ruf_city[city][time_index] = finite_mean(@view(field[city_groups[city]]), valid_forecast)
        end
        (time_index == 1 || time_index % 24 == 0 || time_index == length(VALID_TIMES)) &&
            println("  RUF valid time $time_index/$(length(VALID_TIMES)): $valid_time")
    end
    isempty(missing_cycles) || println("Missing RUF cycles: ", join(missing_cycles, ", "))

    matplotlib = pyimport("matplotlib")
    configure_matplotlib!(matplotlib)
    plt = pyimport("matplotlib.pyplot")
    mdates = pyimport("matplotlib.dates")
    ticker = pyimport("matplotlib.ticker")
    ccrs = pyimport("cartopy.crs")
    shpreader = pyimport("cartopy.io.shapereader")
    patches = pyimport("matplotlib.patches")
    lines = pyimport("matplotlib.lines")
    fm = matplotlib.font_manager

    font_regular_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
    font_bold_path = "/home/zsm/.julia/conda/3/x86_64/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
    figure_fonts = font_properties(fm, font_regular_path, font_bold_path)
    font_regular = figure_fonts.tick
    font_bold = figure_fonts.panel_title
    font_axis = figure_fonts.axis_label
    font_legend = figure_fonts.figure_legend

    ruf_color = RUF_COLOR
    obs_color = OBS_COLOR
    change_color = AURORA_COLOR
    city_color = "red"

    figure = plt.figure(figsize=(FIGURE_WIDTH, 3.85), dpi=600)
    top_x = [0.035, 0.23, 0.425, 0.62, 0.815]
    panel_width, panel_height = 0.145, 0.15
    axes_by_city = Dict{String, Any}()
    panel_side = Dict{String, Symbol}()

    for (index, city) in enumerate(CITY_PANEL_GROUPS.top)
        axes_by_city[city] = figure.add_axes([top_x[index], 0.815, panel_width, panel_height])
        panel_side[city] = :top
    end
    for (index, city) in enumerate(CITY_PANEL_GROUPS.left)
        axes_by_city[city] = figure.add_axes([0.035, 0.625 - (index - 1) * 0.19, panel_width, panel_height])
        panel_side[city] = :left
    end
    for (index, city) in enumerate(CITY_PANEL_GROUPS.right)
        axes_by_city[city] = figure.add_axes([0.82, 0.625 - (index - 1) * 0.19, panel_width, panel_height])
        panel_side[city] = :right
    end
    for (index, city) in enumerate(CITY_PANEL_GROUPS.bottom)
        axes_by_city[city] = figure.add_axes([top_x[index], 0.025, panel_width, panel_height])
        panel_side[city] = :bottom
    end

    plate_carree = ccrs.PlateCarree()
    projection = ccrs.LambertConformal(
        central_longitude=105, central_latitude=83, standard_parallels=[25, 45], globe=nothing
    )
    map_axis = figure.add_axes([0.27, 0.245, 0.46, 0.51], projection=projection)
    boundary_path = "/home/zsm/codeLib/MyPlotM/china_basic_map/bou2_4l.dbf"
    map_axis.add_geometries(
        shpreader.Reader(boundary_path).geometries(), plate_carree,
        facecolor="none", edgecolor="#777777", zorder=1, linewidth=0.6,
    )
    map_axis.set_extent([72.3, 136.3, 16.9, 54.2], crs=plate_carree)
    map_axis.gridlines(
        crs=plate_carree, draw_labels=false, linewidth=GRID_LINEWIDTH,
        color="#888888", alpha=0.55, linestyle="--", zorder=0,
    )
    map_axis.text(
        0.02, 0.97, "(a) Cities position", transform=map_axis.transAxes,
        va="top", ha="left", fontproperties=font_bold, zorder=20,
    )
    map_axis.scatter(
        [city_lon[city] for city in CITIES], [city_lat[city] for city in CITIES],
        s=6, color=city_color, edgecolors="none", transform=plate_carree,
        zorder=20, label="city",
    )

    # South China Sea inset, matching the paper-map style used by the other scripts.
    map_position = map_axis.get_position()
    figure_aspect = figure.figbbox.height / figure.figbbox.width
    inset_width = map_position.width / 5.5
    inset_height = inset_width * 1.75 / figure_aspect
    inset_axis = figure.add_axes(
        [map_position.xmax - inset_width, map_position.ymin, inset_width, inset_height],
        projection=projection,
    )
    inset_axis.set_anchor("SE")
    inset_axis.set_extent([106, 120, 2, 24], crs=plate_carree)
    inset_axis.add_geometries(
        shpreader.Reader(boundary_path).geometries(), plate_carree,
        facecolor="none", edgecolor="#888888", zorder=1, linewidth=0.6,
    )
    inset_cities = [city for city in CITIES if 106 <= city_lon[city] <= 120 && 2 <= city_lat[city] <= 24]
    inset_axis.scatter(
        [city_lon[city] for city in inset_cities], [city_lat[city] for city in inset_cities],
        s=4, color=city_color, edgecolors="none", transform=plate_carree, zorder=10,
    )

    for (city_index, city) in enumerate(CITIES)
        axis = axes_by_city[city]
        ruf_values = ruf_city[city]
        obs_values = observation_city[city]
        axis.plot(VALID_TIMES, ruf_values, color=ruf_color, linewidth=DATA_LINEWIDTH,
            linestyle="-", label="RUF", zorder=3)
        axis.plot(VALID_TIMES, obs_values, color=obs_color, linewidth=DATA_LINEWIDTH,
            linestyle="--", label="obs", zorder=2)

        ruf_change_points, _ = detect_starting_points(
            copy(ruf_values), DERIVATIVE_THRESHOLD, CONCENTRATION_THRESHOLD
        )
        if !isempty(ruf_change_points)
            change_index = first(ruf_change_points)
            if isfinite(ruf_values[change_index])
                axis.scatter(
                    [VALID_TIMES[change_index]], [ruf_values[change_index]], marker="*", s=48,
                    color=change_color, edgecolors="none", label="first change point", zorder=9,
                )
            end
        end

        obs_change_points, _ = detect_starting_points(
            copy(obs_values), DERIVATIVE_THRESHOLD, CONCENTRATION_THRESHOLD
        )
        if !isempty(obs_change_points)
            change_index = first(obs_change_points)
            if isfinite(obs_values[change_index])
                axis.scatter(
                    [VALID_TIMES[change_index]], [obs_values[change_index]], marker="*", s=48,
                    color=change_color, edgecolors="none", label="first change point", zorder=10,
                )
            end
        end

        paired = vcat(
            Float64[value for value in ruf_values if isfinite(value)],
            Float64[value for value in obs_values if isfinite(value)],
        )
        ymax = isempty(paired) ? 1.0 : maximum(paired)
        axis.set_ylim(0, ymax <= 0 ? 1.0 : ymax * 1.12)
        axis.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, integer=true))
        axis.set_xlim(FIRST_VALID, LAST_VALID)
        axis.set_xticks([FIRST_VALID, LAST_VALID])
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        axis.tick_params(axis="both", labelsize=TICK_SIZE, width=AXIS_LINEWIDTH,
            length=2.6, pad=1.2)
        for label in axis.get_xticklabels()
            label.set_fontproperties(font_regular)
        end
        for label in axis.get_yticklabels()
            label.set_fontproperties(font_regular)
        end
        xticklabels = axis.get_xticklabels()
        !isempty(xticklabels) && xticklabels[1].set_ha("left")
        length(xticklabels) > 1 && xticklabels[end].set_ha("right")
        axis.text(
            0.02, 0.96, "($(Char('b' + city_index - 1))) $(CITY_ENGLISH[city])",
            transform=axis.transAxes,
            va="top", ha="left", fontproperties=font_bold, zorder=10,
        )
        side = panel_side[city]
        anchor = side == :top ? (0.5, 0.0) : side == :bottom ? (0.5, 1.0) :
            side == :left ? (1.0, 0.5) : (0.0, 0.5)
        transform_a = plate_carree["_as_mpl_transform"](map_axis)
        connector = patches.ConnectionPatch(
            (city_lon[city], city_lat[city]), anchor,
            transform_a, "axes fraction", axesA=map_axis, axesB=axis,
            arrowstyle="-", linestyle="--", mutation_scale=6,
            linewidth=0.6, color="#777777", clip_on=false, zorder=30,
        )
        figure.add_artist(connector)
    end

    figure.text(-0.018, 0.50, raw"PM$_{10}$ (μg m$^{-3}$)", rotation=90,
        va="center", ha="center", fontproperties=font_axis)

    legend_handles = [
        lines.Line2D([], [], color=ruf_color, linewidth=DATA_LINEWIDTH),
        lines.Line2D([], [], color=obs_color, linewidth=DATA_LINEWIDTH, linestyle="--"),
        lines.Line2D([], [], color=change_color, marker="*", linestyle="none", markersize=6.0),
        lines.Line2D([], [], color=city_color, marker="o", linestyle="none", markersize=2.7),
    ]
    figure.legend(
        legend_handles, ["RUF", "obs", "first change point", "city"],
        loc="lower center", bbox_to_anchor=(0.5, 0.19), ncol=4,
        frameon=true, fancybox=false, edgecolor="#777777", framealpha=1.0,
        prop=font_legend, handlelength=2.4, columnspacing=1.4,
    )

    png_path = joinpath(OUTPUT_DIR, "$OUTPUT_BASE.png")
    pdf_path = joinpath(OUTPUT_DIR, "$OUTPUT_BASE.pdf")
    figure.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    figure.savefig(pdf_path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(figure)
    println("Wrote: $png_path")
    println("Wrote: $pdf_path")
end

main()
