#!/usr/bin/env julia

push!(LOAD_PATH, "/home/zsm/codeLib")

using CSV
using DataFrames
using Dates
using PyCall
using StaDataM
using TimeZones

const BJT = tz"Asia/Shanghai"
const EVAL_FIRST = ZonedDateTime(DateTime(2025, 4, 10, 0), BJT)
const EVAL_LAST = ZonedDateTime(DateTime(2025, 4, 16, 0), BJT)
const EVAL_TIMES = collect(EVAL_FIRST:Hour(1):EVAL_LAST)
const CITIES = ["呼和浩特", "太原", "北京", "洛阳", "南京", "海口"]
const OBS_FILE = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const RUF_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const AURORA_ROOT = "/home/zsm/IETM/Aurora"
const AI_ROOT = "/home/zsm/IETM/AI-GAMFS/inference"
const OUTPUT_DIR = "/home/zsm/IETM/figure/PapreFigure/20260818"
const OUTPUT_FILE = joinpath(OUTPUT_DIR, "ParperTable_ThreeModelCompare.csv")
const EVENT_THRESHOLD = 200.0
const STAMP_FORMAT = dateformat"yyyymmddHH"

# MATLAB character metadata can arrive as vectors of one-character strings.
clean_text(value) = value === missing ? "" : value isa AbstractVector ? join(string.(value)) : strip(string(value))
valid_obs(value) = isfinite(value) && 0.0 <= value < 2000.0
valid_forecast(value) = isfinite(value) && value >= 0.0
time_key(value) = DateTime(astimezone(value, BJT))
time_iso(value) = Dates.format(time_key(value), dateformat"yyyy-mm-ddTHH:MM:SS")

function update_contingency!(counts::Vector{Int}, reference::AbstractVector, forecast::AbstractVector)
    length(reference) == length(forecast) || error("Observation and forecast station vectors have different lengths")
    for (observed, predicted) in zip(reference, forecast)
        valid_obs(observed) && valid_forecast(predicted) || continue
        observed_event = observed > EVENT_THRESHOLD
        forecast_event = predicted > EVENT_THRESHOLD
        if forecast_event && observed_event
            counts[1] += 1  # hit
        elseif forecast_event
            counts[2] += 1  # false alarm
        elseif observed_event
            counts[3] += 1  # miss
        else
            counts[4] += 1  # correct rejection
        end
    end
    return nothing
end

safe_ratio(numerator::Int, denominator::Int) = denominator == 0 ? NaN : numerator / denominator

function main()
    mkpath(OUTPUT_DIR)
    println("Evaluation window (BJT): $(time_iso(EVAL_FIRST)) -- $(time_iso(EVAL_LAST))")
    println("Dust-event threshold: station PM10 > $(EVENT_THRESHOLD) ug m^-3")

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
    length(unique(station_ids)) == length(station_ids) || error("Duplicated national-control station IDs")
    city_positions = [findall(==(city), station_cities) for city in CITIES]
    any(isempty, city_positions) && error("At least one requested city has no national-control station")
    println("National-control station counts: ",
        join(["$(CITIES[i])=$(length(city_positions[i]))" for i in eachindex(CITIES)], ", "))
    println("Total selected national-control stations: $(length(station_ids))")

    pm10_position = findfirst(==("PM10"), uppercase.(clean_text.(observations.varNames)))
    pm10_position === nothing && error("PM10 is absent from $OBS_FILE")
    obs_time_lookup_all = Dict(time_key(t) => i for (i, t) in enumerate(observations.times))
    obs_indices = [get(obs_time_lookup_all, time_key(t), 0) for t in EVAL_TIMES]
    any(==(0), obs_indices) && error("Observation hours are missing inside the requested evaluation window")
    obs_field = Float64.(observations.staDataAll[obs_indices, pm10_position, selected_positions, 1])
    obs_row_lookup = Dict(time_key(t) => i for (i, t) in enumerate(EVAL_TIMES))

    counts = Dict(
        "AI-GAMFS" => zeros(Int, 4),
        "Aurora" => zeros(Int, 4),
        "RUF" => zeros(Int, 4),
    )

    function consume_forecast!(model::String, valid_time::DateTime, values::AbstractVector{<:Real})
        obs_index = get(obs_row_lookup, valid_time, 0)
        obs_index == 0 && return
        length(values) == length(station_ids) || error("$model station-vector length changed at $valid_time")
        update_contingency!(counts[model], @view(obs_field[obs_index, :]), values)
    end

    # RUF: hourly cycles, evaluated at lead 3 h.
    println("Reading RUF lead 3 h")
    ruf_cycles = collect(EVAL_FIRST - Hour(3):Hour(1):EVAL_LAST - Hour(3))
    template_stamp = Dates.format(first(ruf_cycles), STAMP_FORMAT)
    template_path = joinpath(RUF_ROOT, template_stamp, "staData.$template_stamp.mat")
    isfile(template_path) || error("Missing RUF template: $template_path")
    ruf_template = staDATA(template_path)
    ruf_ids = clean_text.(ruf_template.staInfo.staID)
    ruf_lookup = Dict(id => i for (i, id) in enumerate(ruf_ids))
    all(haskey(ruf_lookup, id) for id in station_ids) || error("RUF is missing selected national-control stations")
    ruf_positions = [ruf_lookup[id] for id in station_ids]
    for (cycle_index, cycle) in enumerate(ruf_cycles)
        stamp = Dates.format(cycle, STAMP_FORMAT)
        path = joinpath(RUF_ROOT, stamp, "staData.$stamp.mat")
        if !isfile(path)
            @warn "Missing RUF cycle" model="RUF" init_bjt=time_iso(cycle) lead_hour=3 valid_bjt=time_iso(cycle + Hour(3)) path
            continue
        end
        object = staDATA(path)
        pm_position = findfirst(==("PM10"), uppercase.(clean_text.(object.varNames)))
        pm_position === nothing && error("PM10 absent in $path")
        lead_position = findfirst(==(3), Int.(object.leadTimes))
        if lead_position === nothing
            @warn "Missing RUF lead" model="RUF" init_bjt=time_iso(cycle) lead_hour=3 valid_bjt=time_iso(cycle + Hour(3)) path
            continue
        end
        field = Float64.(object.staDataAll[1, pm_position, ruf_positions, lead_position])
        consume_forecast!("RUF", time_key(cycle + Hour(3)), field)
        (cycle_index == 1 || cycle_index % 24 == 0 || cycle_index == length(ruf_cycles)) &&
            println("  RUF cycle $cycle_index/$(length(ruf_cycles)): $stamp")
    end

    # Reuse the established ThreeModelEva NetCDF reader, PM10 definitions,
    # bilinear interpolation, and kg m^-3 to ug m^-3 conversion.
    py"""
import sys, re
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
sys.path.insert(0, "/home/zsm/IETM/code/ThreeModelEva")
from three_model_eval_netcdf_plot import NetCDF, BilinearSampler, AI_COMPONENTS

def _stack_or_empty(rows, station_count):
    return np.vstack(rows) if rows else np.empty((0, station_count), dtype=np.float64)

def paper_table_extract_aurora(root, first_iso, last_iso, station_lats, station_lons):
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
        selected = [lead for lead in (12,) if first <= init_bjt + timedelta(hours=lead) <= last]
        if not selected:
            continue
        with NetCDF(path) as ds:
            if sampler is None:
                sampler = BilinearSampler(ds.read_double("lat"), ds.read_double("lon"), lats, lons)
            lookup = {int(v): i for i, v in enumerate(ds.read_int("lead_time").tolist())}
            pm10 = ds.read_float("PM10")[0]
            for lead in selected:
                valid_bjt = init_bjt + timedelta(hours=lead)
                if lead not in lookup:
                    print(f"MISSING|model=Aurora|init_bjt={init_bjt.isoformat()}|lead={lead}|valid_bjt={valid_bjt.isoformat()}|file={path}", flush=True)
                    continue
                rows_t.append(valid_bjt.isoformat())
                rows_l.append(lead)
                rows_v.append(sampler.sample(pm10[lookup[lead]]) * 1.0e9)
        print("Aurora", path.name, "done", flush=True)
    return {
        "times": rows_t,
        "leads": np.asarray(rows_l, dtype=np.int64),
        "values": _stack_or_empty(rows_v, len(lats)),
    }

def paper_table_extract_ai(root, first_iso, last_iso, station_lats, station_lons):
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
        init_bjt = init_label_utc - timedelta(minutes=30) + timedelta(hours=8)
        lead = 3
        valid_bjt = init_bjt + timedelta(hours=lead)
        if not first <= valid_bjt <= last:
            continue
        valid_label_utc = init_label_utc + timedelta(hours=lead)
        filename = f"AI_GAMFS.{init_label_utc:%Y%m%d_%H%M}+{valid_label_utc:%Y%m%d_%H%M}.V01.nc"
        path = directory / filename
        if not path.is_file():
            print(f"MISSING|model=AI-GAMFS|init_bjt={init_bjt.isoformat()}|lead={lead}|valid_bjt={valid_bjt.isoformat()}|file={path}", flush=True)
            continue
        with NetCDF(path) as ds:
            if sampler is None:
                sampler = BilinearSampler(ds.read_double("latitude"), ds.read_double("longitude"), lats, lons)
            forecast = np.zeros(len(lats), dtype=np.float64)
            for component in AI_COMPONENTS:
                forecast += sampler.sample(ds.read_float(component)[0])
            forecast *= 1.0e9
        rows_t.append(valid_bjt.isoformat())
        rows_l.append(lead)
        rows_v.append(forecast)
        print("AI-GAMFS", directory.name, "done", flush=True)
    return {
        "times": rows_t,
        "leads": np.asarray(rows_l, dtype=np.int64),
        "values": _stack_or_empty(rows_v, len(lats)),
    }
"""

    function consume_python_rows!(model::String, result, expected_leads::Set{Int})
        times = convert(Vector{String}, result["times"])
        leads = Int.(vec(Array(result["leads"])))
        values = Float64.(Array(result["values"]))
        size(values, 1) == length(times) == length(leads) || error("$model extraction arrays are inconsistent")
        for row in eachindex(times)
            leads[row] in expected_leads || error("Unexpected $model lead $(leads[row])")
            valid_time = DateTime(times[row], dateformat"yyyy-mm-ddTHH:MM:SS")
            consume_forecast!(model, valid_time, @view(values[row, :]))
        end
        println("$model forecast records used: $(length(times))")
    end

    first_iso = time_iso(EVAL_FIRST)
    last_iso = time_iso(EVAL_LAST)
    println("Reading Aurora lead 12 h")
    aurora_result = py"paper_table_extract_aurora"(AURORA_ROOT, first_iso, last_iso, station_lats, station_lons)
    consume_python_rows!("Aurora", aurora_result, Set([12]))

    println("Reading AI-GAMFS lead 3 h")
    ai_result = py"paper_table_extract_ai"(AI_ROOT, first_iso, last_iso, station_lats, station_lons)
    consume_python_rows!("AI-GAMFS", ai_result, Set([3]))

    models = ["AI-GAMFS", "Aurora", "RUF"]
    evaluated_leads = ["3 h", "12 h", "3 h"]
    false_alarm_rates = Float64[]
    miss_rates = Float64[]
    threat_scores = Float64[]
    for model in models
        hit, false_alarm, miss, correct_rejection = counts[model]
        false_alarm_rate = safe_ratio(false_alarm, hit + false_alarm)
        miss_rate = safe_ratio(miss, hit + miss)
        threat_score = safe_ratio(hit, hit + false_alarm + miss)
        println("$model contingency counts: hit=$hit, false_alarm=$false_alarm, miss=$miss, correct_rejection=$correct_rejection, valid_pairs=$(sum(counts[model]))")
        println("$model unrounded scores: false_alarm_rate=$false_alarm_rate, miss_rate=$miss_rate, TS=$threat_score")
        push!(false_alarm_rates, round(false_alarm_rate; digits=2))
        push!(miss_rates, round(miss_rate; digits=2))
        push!(threat_scores, round(threat_score; digits=2))
    end

    result = DataFrame(
        "Model" => models,
        "Evaluated lead time" => evaluated_leads,
        "False alarm rate" => false_alarm_rates,
        "Miss rate" => miss_rates,
        "TS" => threat_scores,
    )
    CSV.write(OUTPUT_FILE, result)
    println(result)
    println("Wrote: $OUTPUT_FILE")
end

main()
