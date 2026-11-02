using Dates
using GradDataM
using JLD2
using NaNStatistics
using Serialization
using StaDataM
using Statistics
using TimeZones

const FIG5_OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const FIG5_RUF2D_FILE = "/home/zsm/IETM/RUF_2d/DATA/RUF/2025041118/staData.2025041118.mat"
const FIG5_RUF3D_FILE = "/home/zsm/IETM/RUF_DryDeposition/2025041118/staData.2025041118.mat"
const FIG5_RUF3D_CTL = "/home/zsm/IETM/RUF_DryDeposition/2025041118/ruf.ctl"
const FIG5_DATAGRID = "/home/zsm/IETM/RUF_2d/datagrid.jld2"
const FIG5_HEIGHT_CTL = "/home/zsm/IETM/AuxData/grid.d2.ctl"
const FIG5_CACHE = "/home/zsm/IETM/figure/PapreFigure/20260817/figure5_zsm_prepared_cache_20260817.bin"
const FIG5_BJT = tz"Asia/Shanghai"

figure5_clean(value) = value === missing ? String() : value isa AbstractVector ? join(string.(value)) : strip(string(value))
figure5_city(value) = figure5_clean(value)

function figure5_finite_mean(values)
    finite_values = Float64.(values[isfinite.(values)])
    return isempty(finite_values) ? NaN : mean(finite_values)
end

function figure5_finite_max(values)
    finite_values = Float64.(values[isfinite.(values)])
    return isempty(finite_values) ? NaN : maximum(finite_values)
end

function get_truncated_cmap(base_cmap_name, min_val=0.0, max_val=1.0, n=256)
    base_cmap = plt.get_cmap(base_cmap_name)
    colors_resampled = base_cmap(np.linspace(min_val, max_val, n))
    return mcolors.LinearSegmentedColormap.from_list("trunc_afmhot", colors_resampled)
end

println("Loading Figure 5 observation and station metadata")
staObs = staDATA(FIG5_OBS_FILE)
ruf3d_object = staDATA(FIG5_RUF3D_FILE)
full_station_ids = figure5_clean.(ruf3d_object.staInfo.staID)
full_city_labels = figure5_city.(ruf3d_object.staInfo.city)
valid_station_indices = findall(!isempty, full_city_labels)
staInfo = ruf3d_object.staInfo[valid_station_indices, :]
city_labels = full_city_labels[valid_station_indices]
PlaCar = ccrs.PlateCarree()
lonM, latM = load(FIG5_DATAGRID, "lonM", "latM")

cities = [
    "呼和浩特", "太原", "长治", "晋城", "洛阳", "南阳", "襄阳", "宜昌",
    "常德", "邵阳", "永州", "桂林", "梧州", "茂名", "海口",
]
cities_error = ["太原", "洛阳", "常德", "海口"]
cities_all = unique(city_labels)
all_city_groups = [findall(==(city), city_labels) for city in cities_all]
all_city_lookup = Dict(city => index for (index, city) in enumerate(cities_all))
route_groups = [all_city_groups[all_city_lookup[city]] for city in cities]
selected_groups = [all_city_groups[all_city_lookup[city]] for city in cities_error]

cycle_time = DateTime(2025, 4, 11, 18)
valid_times = collect(cycle_time:Hour(1):(cycle_time + Hour(72)))
t_obs = ZonedDateTime.(valid_times, FIG5_BJT)
t2plot = valid_times
t2plot_label = Dates.format.(valid_times, "mm-dd")

observation_ids = figure5_clean.(staObs.staInfo.staID)
observation_position = Dict(id => index for (index, id) in enumerate(observation_ids))
selected_station_ids = full_station_ids[valid_station_indices]
all(haskey(observation_position, id) for id in selected_station_ids) ||
    error("forecast station IDs are absent from the observation file")
selected_observation_positions = [observation_position[id] for id in selected_station_ids]
pm10_position = findfirst(==("PM10"), uppercase.(figure5_clean.(staObs.varNames)))
pm10_position === nothing && error("PM10 is absent from observation file")
observation_time_lookup = Dict(
    DateTime(astimezone(time, FIG5_BJT)) => index for (index, time) in enumerate(staObs.times)
)
observation_indices = [observation_time_lookup[time] for time in valid_times]
observation_field = Float64.(
    staObs.staDataAll[observation_indices, pm10_position, selected_observation_positions, 1]
)
observation_field[
    .!isfinite.(observation_field) .| (observation_field .< 0) .| (observation_field .>= 2000)
] .= NaN

obs = Vector{Vector{Float64}}(undef, length(cities_error))
for (city_index, positions) in enumerate(selected_groups)
    obs[city_index] = vec(nanmean(@view(observation_field[:, positions]), dims=2))
end
mean_obs = figure5_finite_mean.(obs)
max_obs_selected = figure5_finite_max.(obs)

max_obs = Vector{Float64}(undef, length(cities_all))
cities_lon = Vector{Float64}(undef, length(cities_all))
cities_lat = Vector{Float64}(undef, length(cities_all))
for (city_index, positions) in enumerate(all_city_groups)
    series = vec(nanmean(@view(observation_field[:, positions]), dims=2))
    max_obs[city_index] = figure5_finite_max(series)
    cities_lon[city_index] = figure5_finite_mean(Float64.(staInfo.lon[positions]))
    cities_lat[city_index] = figure5_finite_mean(Float64.(staInfo.lat[positions]))
end

if isfile(FIG5_CACHE)
    println("Loading cached Figure 5 single-cycle series and vertical profiles")
    cached = deserialize(FIG5_CACHE)
    ietm2d_Data = cached.ietm2d_Data
    ietm3d_Data = cached.ietm3d_Data
    height = cached.height
    vertical_profiles = cached.vertical_profiles
else
    println("Preparing matched RUF2D/RUF3D cycle 2025041118")
    ruf2d_object = staDATA(FIG5_RUF2D_FILE)
    figure5_clean.(ruf2d_object.staInfo.staID) == full_station_ids ||
        error("RUF2D and RUF3D station ordering differs")
    Int.(ruf2d_object.leadTimes) == collect(0:72) || error("unexpected RUF2D lead times")
    Int.(ruf3d_object.leadTimes) == collect(0:72) || error("unexpected RUF3D lead times")
    field2d = Float64.(ruf2d_object.staDataAll[1, 1, valid_station_indices, :])
    field3d = Float64.(ruf3d_object.staDataAll[1, 1, valid_station_indices, :])
    field2d[.!isfinite.(field2d) .| (field2d .< 0)] .= NaN
    field3d[.!isfinite.(field3d) .| (field3d .< 0)] .= NaN
    ietm2d_Data = [
        vec(nanmean(@view(field2d[positions, :]), dims=1)) for positions in selected_groups
    ]
    ietm3d_Data = [
        vec(nanmean(@view(field3d[positions, :]), dims=1)) for positions in selected_groups
    ]

    route_locations = [
        (
            figure5_finite_mean(Float64.(staInfo.lon[positions])),
            figure5_finite_mean(Float64.(staInfo.lat[positions])),
        ) for positions in route_groups
    ]
    nearest = [
        argmin((Float64.(lonM) .- lon).^2 .+ (Float64.(latM) .- lat).^2)
        for (lon, lat) in route_locations
    ]

    height_reader = gradData(FIG5_HEIGHT_CTL)
    height = fill(1000.0, 13, length(cities) + 1)
    surface_height = read(height_reader, "heiz", 1)
    for level in 1:13
        level_height = read(height_reader, "heiz", level)
        for city_index in eachindex(cities)
            height[level, city_index] =
                level_height[nearest[city_index]] - surface_height[nearest[city_index]]
        end
    end
    height = reverse(height, dims=1)

    vertical_times = [
        DateTime(2025, 4, 11, 22), DateTime(2025, 4, 12, 2),
        DateTime(2025, 4, 12, 6), DateTime(2025, 4, 12, 10),
        DateTime(2025, 4, 12, 14), DateTime(2025, 4, 12, 18),
    ]
    vertical_profiles = Matrix{Float64}[]
    grid_reader = gradData(FIG5_RUF3D_CTL)
    for valid_time in vertical_times
        mesh = read(setTime!(grid_reader, ZonedDateTime(valid_time, FIG5_BJT)), "PM10")
        profile = fill(NaN, size(mesh, 3), length(cities))
        for city_index in eachindex(cities), level in axes(mesh, 3)
            profile[level, city_index] = mesh[nearest[city_index], level]
        end
        push!(vertical_profiles, profile)
    end
    serialize(FIG5_CACHE, (
        ietm2d_Data=ietm2d_Data,
        ietm3d_Data=ietm3d_Data,
        height=height,
        vertical_profiles=vertical_profiles,
    ))
    println("Saved Figure 5 cache: $FIG5_CACHE")
end

mean_2d = figure5_finite_mean.(ietm2d_Data)
mean_3d = figure5_finite_mean.(ietm3d_Data)
max_2d = figure5_finite_max.(ietm2d_Data)
max_3d = figure5_finite_max.(ietm3d_Data)

ax5data, ax6data, ax7data, ax8data, ax9data, ax10data = vertical_profiles
xtick = collect(1.5:2:15.5)
allpos = [
    [figure5_finite_mean(Float64.(staInfo.lon[group])) for group in route_groups],
    [figure5_finite_mean(Float64.(staInfo.lat[group])) for group in route_groups],
]
city2plot = cities[1:2:15]
citypos = [
    [figure5_finite_mean(Float64.(staInfo.lon[all_city_groups[all_city_lookup[city]]])) for city in city2plot],
    [figure5_finite_mean(Float64.(staInfo.lat[all_city_groups[all_city_lookup[city]]])) for city in city2plot],
]
ty_group = all_city_groups[all_city_lookup["太原"]]
mm_group = all_city_groups[all_city_lookup["茂名"]]
ty_pos = [
    figure5_finite_mean(Float64.(staInfo.lon[ty_group])),
    figure5_finite_mean(Float64.(staInfo.lat[ty_group])),
]
mm_pos = [
    figure5_finite_mean(Float64.(staInfo.lon[mm_group])),
    figure5_finite_mean(Float64.(staInfo.lat[mm_group])),
]

println("Figure 5 new-data preparation complete")
