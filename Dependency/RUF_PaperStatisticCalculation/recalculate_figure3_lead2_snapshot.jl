#!/usr/bin/env julia

pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "codeLib")))

using CSV
using DataFrames
using Dates
using NaNStatistics
using StaDataM
using TimeZones

const BJT = tz"Asia/Shanghai"
const STAMP = dateformat"yyyymmddHH"
const OBS_FILE = get(ENV, "OBS_FILE", "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat")
const FC_ROOT = get(ENV, "FC_ROOT", "/home/zsm/IETM/RUF_DryDeposition")
const OUTPUT_DIR = get(ENV, "OUTPUT_DIR", joinpath(@__DIR__, "output"))
const THRESHOLD = 200.0

clean(value) = value === missing ? "" : strip(string(value))
stamp(value) = Dates.format(value, STAMP)

function main()
    observations = staDATA(OBS_FILE)
    pm10_position = findfirst(==("PM10"), uppercase.(observations.varNames))
    valid_times = [DateTime(2025, 4, 11, hour) for hour in (10, 12, 14, 16)]
    first_cycle = first(valid_times) - Hour(2)
    template = staDATA(joinpath(FC_ROOT, stamp(first_cycle), "staData.$(stamp(first_cycle)).mat"))
    forecast_ids = clean.(template.staInfo.staID)
    obs_position = Dict(clean(id) => index for (index, id) in enumerate(observations.staInfo.staID))
    matched_obs = [obs_position[id] for id in forecast_ids]
    obs_lookup = Dict(DateTime(astimezone(time, BJT)) => index for (index, time) in enumerate(observations.times))
    city_labels = clean.(template.staInfo.city)
    cities = sort(filter(!isempty, unique(city_labels)))
    city_groups = [findall(==(city), city_labels) for city in cities]
    rows = NamedTuple[]
    for valid_time in valid_times
        cycle = valid_time - Hour(2)
        path = joinpath(FC_ROOT, stamp(cycle), "staData.$(stamp(cycle)).mat")
        object = staDATA(path)
        forecast = Float64.(object.staDataAll[1, 1, :, 3])
        observation = Float64.(observations.staDataAll[obs_lookup[valid_time], pm10_position, matched_obs, 1])
        city_forecast = [nanmean(forecast[positions]) for positions in city_groups]
        city_observation = [nanmean(observation[positions]) for positions in city_groups]
        valid = isfinite.(city_forecast) .& isfinite.(city_observation)
        obs_event = valid .& (city_observation .> THRESHOLD)
        fc_event = valid .& (city_forecast .> THRESHOLD)
        hit = count(obs_event .& fc_event)
        false_alarm = count(valid .& .!obs_event .& fc_event)
        miss = count(valid .& obs_event .& .!fc_event)
        ts = hit / (hit + false_alarm + miss)
        push!(rows, (
            valid_time_bjt=stamp(valid_time), lead_hour=2, valid_cities=count(valid),
            hit=hit, false_alarm=false_alarm, miss=miss, ts=ts,
        ))
    end
    result = DataFrame(rows)
    CSV.write(joinpath(OUTPUT_DIR, "figure3_allstations_snapshot_ts_lead2.csv"), result)
    show(stdout, MIME("text/plain"), result); println()
end

main()
