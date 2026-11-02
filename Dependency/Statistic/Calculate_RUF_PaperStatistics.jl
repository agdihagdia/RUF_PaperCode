#!/usr/bin/env julia
# Run with Julia 1.11 and --project=/home/zsm/IETM/code/RUF_Module.
# Reads station MAT files only. No interpolation, forecast, or observation files are modified.
# The first 12 forecast hours mean leads 1:12, excluding the assimilated lead 0.
using Dates, MAT, NaNStatistics, Printf, SHA, Statistics, TOML, TimeZones
using Test

const BJT = tz"Asia/Shanghai"
const STAMP = dateformat"yyyymmddHH"
const DETECTOR = joinpath(@__DIR__, "..", "change_point_detect.jl")
const OBS_DEFAULT = "/home/zsm/数据库运维-Julia/data/verified/qc202504.mat"
const ALL_ROOT = "/home/zsm/IETM/RUF_DryDeposition"
const NATION_ROOT = "/home/zsm/IETM/RUF_DryDeposition_Nation/Ruf_forcast"
const BASE_ROOTS = ["/media/dataLL/zsm/paper_data/baseFolder", "/media/dataLL/zsm/BaseFolder_EPICC"]
const TOLERANCE_H = 6
const GRADIENT_THRESHOLD = 0.1
const CONCENTRATION_THRESHOLD = 200.0
const OBS_UPPER = 2000.0

# Load only function definitions from the manuscript detector, without importing
# its unrelated station IO dependencies or executing a plotting/data-preparation script.
module Arrival
using NaNStatistics
end
for expr in Meta.parseall(read(DETECTOR, String)).args
    if expr isa Expr && expr.head == :function
        Core.eval(Arrival, expr)
    end
end

text(x) = x === missing ? "" : strip(x isa AbstractVector ? join(string.(x)) : string(x))
asvector(x) = x isa AbstractArray ? vec(x) : [x]
bjt(x) = DateTime(astimezone(ZonedDateTime(text(x), "yyyy-mm-ddTHH:MM:SSzzzz"), BJT))
validobs(x) = isfinite(x) && 0 <= x < OBS_UPPER
validfc(x) = isfinite(x) && x >= 0
finite_mean(x) = (v = filter(isfinite, x); isempty(v) ? NaN : mean(v))
ratio_drop(b, r) = isfinite(b) && b > 0 && isfinite(r) ? 100 * (b - r) / b : NaN
timestamp(x) = x === nothing ? "" : Dates.format(x, dateformat"yyyy-mm-ddTHH:MM:SS")

function csvcell(x)
    s = x === nothing ? "" : string(x)
    return occursin(r"[\",\r\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function write_csv(path, rows)
    isempty(rows) && error("Refusing an empty table without a schema: $path")
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in rows
            propertynames(row) == columns || error("Inconsistent columns in $path")
            println(io, join(csvcell.([getproperty(row, k) for k in columns]), ","))
        end
    end
end

function read_station(path)
    a = matread(path)
    info = a["staInfo"]
    ids = text.(asvector(info["staID"]))
    length(unique(ids)) == length(ids) || error("Duplicate station IDs: $path")
    variables = uppercase.(text.(asvector(a["varName"])))
    variable = findfirst(==("PM10"), variables)
    variable === nothing && error("PM10 missing: $path")
    data = a["staDataAll"]
    ndims(data) == 4 || error("Unexpected MAT dimensions: $path")
    times = bjt.(asvector(a["time"]))
    leads = Int.(asvector(a["leadTimes"]))
    size(data) == (length(times), length(variables), length(ids), length(leads)) ||
        error("Inconsistent metadata: $path")
    return (ids=ids, cities=text.(asvector(info["city"])),
        types=text.(asvector(info["type"])), times=times, leads=leads,
        data=data, variable=variable)
end

function positions(ids, selected, label)
    lookup = Dict(id => i for (i, id) in enumerate(ids))
    all(haskey(lookup, id) for id in selected) || error("Missing required station: $label")
    return [lookup[id] for id in selected]
end

function discover(root; utc=false)
    isdir(root) || error("Missing input root: $root")
    rows = NamedTuple[]
    for dir in sort(readdir(root; join=true))
        isdir(dir) || continue
        stamp = basename(dir)
        occursin(r"^\d{10}$", stamp) || continue
        path = joinpath(dir, "staData.$stamp.mat")
        isfile(path) || continue
        init = DateTime(stamp, STAMP) + Hour(utc ? 8 : 0)
        push!(rows, (initialization=init, path=path, root=root))
    end
    return rows
end

function read_cycle(cycle, ids)
    a = read_station(cycle.path)
    a.times == [cycle.initialization] || error("Initialization/time-zone mismatch: $(cycle.path)")
    a.ids == ids || error("Station order changed: $(cycle.path)")
    return a
end

mutable struct Accumulator
    n::Int
    sse::Float64
    cycles::Int
end
Accumulator() = Accumulator(0, 0.0, 0)
rmse(a::Accumulator) = a.n == 0 ? NaN : sqrt(a.sse / a.n)
function add!(a, observation, forecast, mask)
    n = count(mask)
    n == 0 && return
    a.n += n
    a.sse += sum((forecast[mask] .- observation[mask]).^2)
    a.cycles += 1
end

function arrival_comparison(obs, fc)
    op, _ = Arrival.detect_starting_points(Float64.(obs), GRADIENT_THRESHOLD, CONCENTRATION_THRESHOLD)
    fp, _ = Arrival.detect_starting_points(Float64.(fc), GRADIENT_THRESHOLD, CONCENTRATION_THRESHOLD)
    # Correct order: forecast candidate indices precede observation candidates.
    oi, fi = Arrival.filter_starting_points(obs, fc, fp, op)
    ok = Arrival.judge_forcast(obs, fc, op, fp, TOLERANCE_H)
    expected = oi > 0 && fi > 0 && abs(oi - fi) <= TOLERANCE_H
    ok == Int(expected) || error("judge_forcast argument-order regression; inspect $DETECTOR")
    return (oi=oi, fi=fi, success=ok == 1, observed_event=!isempty(op),
        error_h=oi > 0 && fi > 0 ? abs(oi - fi) : NaN)
end

function peak_info(series)
    # A missing hour can conceal the event maximum; do not report an incomplete peak.
    all(isfinite, series) || return (index=0, value=NaN, ties=0)
    value = maximum(series)
    return (index=findfirst(==(value), series), value=value, ties=count(==(value), series))
end

function self_test()
    @testset "Manuscript arrival and RMSE definitions" begin
        obs = fill(100.0, 30); fc = fill(100.0, 30)
        op = [10, 11, 12]; fp = [16, 17, 18]
        obs[op] = [10, 20, 30]; fc[fp] = [10, 20, 30]
        obs[fp] = [30, 20, 10]; fc[op] = [10, 20, 30]
        @test Arrival.filter_starting_points(obs, fc, fp, op) == (10, 16)
        @test Arrival.filter_starting_points(obs, fc, op, fp) == (18, 10)
        @test Arrival.judge_forcast(obs, fc, op, fp, 6) == 1
        @test Arrival.judge_forcast(obs, fc, [10], [17], 6) == 0
        @test Arrival.judge_forcast(obs, fc, Int[], fp, 6) == 0
        @test Arrival.bool2leadtime([1 0 1; 0 0 0]) == reshape([3, 0], 2, 1)
        @test validobs(1999) && !validobs(2000) && !validobs(-1) && !validobs(NaN)
        @test peak_info([1.0, 3.0, 3.0]).ties == 2
        @test peak_info([1.0, NaN]).index == 0
        a = Accumulator()
        add!(a, [1.0, 2.0], [4.0, 6.0], [true, true])
        @test rmse(a) == sqrt(12.5)
        @test ratio_drop(100.0, 75.0) == 25.0
        @test isnan(ratio_drop(0.0, 0.0))
        @test bjt("2025-04-11T12:00:00+00:00") == DateTime(2025, 4, 11, 20)
    end
end

function options(args)
    opts = Dict(
        "--output-dir" => joinpath(@__DIR__, "output", Dates.format(now(), dateformat"yyyymmdd_HHMMSS")),
        "--observation-file" => OBS_DEFAULT,
        "--event-start" => "2025-04-11T00:00:00",
        "--event-end" => "2025-04-14T00:00:00")
    iseven(length(args)) || error("Options must be --name value pairs; or use --self-test")
    for i in 1:2:length(args)
        haskey(opts, args[i]) || error("Unknown option: $(args[i])")
        opts[args[i]] = args[i+1]
    end
    return opts
end

function main(args)
    "--self-test" in args && return self_test()
    opts = options(args)
    self_test()
    detector_hash = bytes2hex(sha256(read(DETECTOR)))
    output = abspath(opts["--output-dir"])
    ispath(output) && error("Output already exists; select a new --output-dir: $output")
    event_start, event_end = DateTime(opts["--event-start"]), DateTime(opts["--event-end"])
    event_end > event_start || error("event-end must be later than event-start")
    event_times = collect(event_start:Hour(1):event_end)
    event_times[end] == event_end || error("Event bounds must have an integer-hour separation")
    event_lookup = Dict(t => i for (i,t) in enumerate(event_times))
    all_cycles, nation_cycles = discover(ALL_ROOT), discover(NATION_ROOT)
    isempty(all_cycles) && error("No All-stations forecasts")
    isempty(nation_cycles) && error("No Nation-only forecasts")
    ruf_min = min(first(all_cycles).initialization, first(nation_cycles).initialization)
    ruf_max = max(last(all_cycles).initialization, last(nation_cycles).initialization)
    # Include daily baseline cycles within the calendar months covered by RUF;
    # retain ALL cycles in those months for the TOC-style independent-cycle result.
    months = Set((year(c.initialization), month(c.initialization)) for c in vcat(all_cycles, nation_cycles))
    baseline_cycles = NamedTuple[]; inventory = NamedTuple[]
    for root in BASE_ROOTS
        found = discover(root; utc=true)
        used = filter(c -> (year(c.initialization), month(c.initialization)) in months, found)
        append!(baseline_cycles, used)
        push!(inventory, (model="Baseline", root=root, discovered_cycles=length(found),
            used_cycles=length(used), first_bjt=isempty(found) ? "" : timestamp(first(found).initialization),
            last_bjt=isempty(found) ? "" : timestamp(last(found).initialization),
            selection="calendar months represented by RUF; other months excluded"))
    end
    sort!(baseline_cycles; by=c -> c.initialization)
    isempty(baseline_cycles) && error("No contemporaneous baseline cycles")
    length(unique(c.initialization for c in baseline_cycles)) == length(baseline_cycles) ||
        error("Duplicate baseline initialization across roots; resolve explicitly")
    for (name,root,cycles) in [("All-stations",ALL_ROOT,all_cycles),("Nation-only",NATION_ROOT,nation_cycles)]
        push!(inventory, (model=name,root=root,discovered_cycles=length(cycles),used_cycles=length(cycles),
            first_bjt=timestamp(first(cycles).initialization),last_bjt=timestamp(last(cycles).initialization),
            selection="all available 10-digit initialization folders with station MAT"))
    end
    println("Loading observation and station metadata")
    obs = read_station(opts["--observation-file"])
    at = read_station(first(all_cycles).path); nt = read_station(first(nation_cycles).path)
    bt = read_station(first(baseline_cycles).path)
    obslookup = Dict(t => i for (i,t) in enumerate(obs.times))
    length(obslookup) == length(obs.times) || error("Duplicate observation times")
    Set([0]) == Set(obs.leads) || error("Expected observations at lead 0")
    all(==("B"), bt.types) || error("Baseline station network is not exclusively national-control")
    obsset, nset, aset = Set(obs.ids), Set(nt.ids), Set(at.ids)
    station_ids = [id for id in bt.ids if id in obsset && id in nset && id in aset]
    isempty(station_ids) && error("No common stations for RMSE")
    model_ids = Dict("All-stations"=>at.ids,"Nation-only"=>nt.ids,"Baseline"=>bt.ids)
    model_pos = Dict(name=>positions(ids,station_ids,name) for (name,ids) in model_ids)
    obspos = positions(obs.ids,station_ids,"observation")
    # Both RUF experiments are evaluated at the same station network, including
    # county/township sites. 'Nation-only' describes assimilation, not evaluation.
    city_ids = [id for (id,city) in zip(at.ids,at.cities) if !isempty(city) && id in nset && id in obsset]
    city_ap = positions(at.ids,city_ids,"All city network")
    city_np = positions(nt.ids,city_ids,"Nation city network")
    city_op = positions(obs.ids,city_ids,"observation city network")
    city_labels = at.cities[city_ap]
    at.cities[city_ap] == nt.cities[city_np] || error("City labels differ between RUF experiments")
    cities = unique(city_labels)
    city_groups = [findall(==(city),city_labels) for city in cities]
    beijing = findfirst(==("北京"),cities)
    beijing === nothing && error("Beijing missing from city network")
    obs_city = fill(NaN,length(event_times),length(cities))
    for (ti,t) in enumerate(event_times)
        oi = get(obslookup,t,0); oi > 0 || error("Observation time missing: $t")
        vals = Float64.(obs.data[oi,obs.variable,city_op,1])
        vals[.!validobs.(vals)] .= NaN
        for (ci,group) in enumerate(city_groups)
            obs_city[ti,ci] = nanmean(vals[group])
        end
    end
    function observation(t)
        oi = get(obslookup,t,0)
        oi > 0 || error("Observation time missing: $t; use the correct observation file")
        return Float64.(obs.data[oi,obs.variable,obspos,1])
    end
    independent = Dict(name=>[Accumulator() for _ in 1:12] for name in keys(model_ids))
    matched = Dict(name=>[Accumulator() for _ in 1:12] for name in keys(model_ids))
    baseline_fields = Dict{Tuple{DateTime,Int},Vector{Float64}}()
    println("RMSE common national-control stations: ",length(station_ids),
        "; city evaluation stations: ",length(city_ids),"; cities: ",length(cities))
    input_manifest = NamedTuple[]
    function record_input(model,path,init)
        st=stat(path)
        push!(input_manifest,(model=model,path=path,initialization_bjt=timestamp(init),
            bytes=st.size,mtime_unix=st.mtime))
    end
    record_input("Observation",opts["--observation-file"],nothing)
    for (j,c) in enumerate(baseline_cycles)
        a=read_cycle(c,bt.ids); lp=Dict(l=>i for (i,l) in enumerate(a.leads))
        for lead in 1:12
            haskey(lp,lead) || error("Baseline missing lead $lead")
            f=Float64.(a.data[1,a.variable,model_pos["Baseline"],lp[lead]])
            o=observation(c.initialization+Hour(lead))
            add!(independent["Baseline"][lead],o,f,validobs.(o).&validfc.(f))
            baseline_fields[(c.initialization,lead)]=f
        end
        record_input("Baseline",c.path,c.initialization)
        (j==1 || j==length(baseline_cycles)) && println("Baseline $j/$(length(baseline_cycles))")
    end
    lookups=Dict("All-stations"=>Dict(c.initialization=>c for c in all_cycles),
        "Nation-only"=>Dict(c.initialization=>c for c in nation_cycles))
    fc_city=Dict(name=>fill(NaN,length(event_times),length(cities),72) for name in keys(lookups))
    city_pos=Dict("All-stations"=>city_ap,"Nation-only"=>city_np)
    cycle_union=sort(collect(union(keys(lookups["All-stations"]),keys(lookups["Nation-only"]))))
    for (j,init) in enumerate(cycle_union)
        fields=Dict{String,Matrix{Float64}}()
        for name in ("All-stations","Nation-only")
            haskey(lookups[name],init) || continue
            c=lookups[name][init]; a=read_cycle(c,model_ids[name])
            lp=Dict(l=>i for (i,l) in enumerate(a.leads))
            all(haskey(lp,l) for l in 1:72) || error("Incomplete 1:72 leads: $(c.path)")
            first12=Float64.(a.data[1,a.variable,model_pos[name],[lp[l] for l in 1:12]])
            fields[name]=first12
            for lead in 1:12
                o=observation(init+Hour(lead)); f=first12[:,lead]
                add!(independent[name][lead],o,f,validobs.(o).&validfc.(f))
            end
            for lead in 1:72
                ti=get(event_lookup,init+Hour(lead),0); ti==0 && continue
                f=Float64.(a.data[1,a.variable,city_pos[name],lp[lead]])
                f[.!validfc.(f)].=NaN
                for (ci,group) in enumerate(city_groups)
                    fc_city[name][ti,ci,lead]=nanmean(f[group])
                end
            end
            record_input(name,c.path,init)
        end
        if length(fields)==2 && haskey(baseline_fields,(init,1))
            for lead in 1:12
                o=observation(init+Hour(lead)); b=baseline_fields[(init,lead)]
                a=fields["All-stations"][:,lead]; n=fields["Nation-only"][:,lead]
                mask=validobs.(o).&validfc.(b).&validfc.(a).&validfc.(n)
                for (name,f) in [("Baseline",b),("All-stations",a),("Nation-only",n)]
                    add!(matched[name][lead],o,f,mask)
                end
            end
        end
        if j==1 || j%24==0 || j==length(cycle_union)
            println("RUF cycle $j/$(length(cycle_union)): $init");flush(stdout);GC.gc(false)
        end
    end
    metric_rows=NamedTuple[]; rmse_summary=NamedTuple[]
    for (scheme,accs) in [("independent_all_cycles",independent),("matched_initialization_and_samples",matched)]
        for name in ("All-stations","Nation-only")
            drops=Float64[]; br=Float64[]; rr=Float64[]
            for lead in 1:12
                b=accs["Baseline"][lead]; r=accs[name][lead]
                b.n>0 && r.n>0 || error("No samples: $scheme/$name/$lead")
                reduction=ratio_drop(rmse(b),rmse(r))
                push!(drops,reduction);push!(br,rmse(b));push!(rr,rmse(r))
                push!(metric_rows,(scheme=scheme,model=name,lead_hour=lead,
                    baseline_rmse=rmse(b),ruf_rmse=rmse(r),reduction_percent=reduction,
                    baseline_pairs=b.n,ruf_pairs=r.n,baseline_cycles=b.cycles,ruf_cycles=r.cycles,
                    common_station_count=length(station_ids)))
            end
            bp=sqrt(sum(a.sse for a in accs["Baseline"])/sum(a.n for a in accs["Baseline"]))
            rp=sqrt(sum(a.sse for a in accs[name])/sum(a.n for a in accs[name]))
            push!(rmse_summary,(scheme=scheme,model=name,first_lead=1,last_lead=12,
                mean_hourly_reduction_percent=mean(drops),
                reduction_of_mean_rmse_percent=ratio_drop(mean(br),mean(rr)),
                pooled_12h_reduction_percent=ratio_drop(bp,rp),
                mean_baseline_rmse=mean(br),mean_ruf_rmse=mean(rr),
                pooled_baseline_rmse=bp,pooled_ruf_rmse=rp))
        end
    end
    println("Calculating corrected arrival matching and Beijing peak timing")
    city_rows=NamedTuple[]; lead_rows=NamedTuple[]; fraction_rows=NamedTuple[]
    arrival_max=Dict(name=>zeros(Int,length(cities)) for name in keys(lookups))
    events=falses(length(cities))
    for (ci,city) in enumerate(cities)
        o=obs_city[:,ci]
        op,_=Arrival.detect_starting_points(o,GRADIENT_THRESHOLD,CONCENTRATION_THRESHOLD)
        events[ci]=!isempty(op)
        for name in ("All-stations","Nation-only")
            flags=zeros(Int,1,72)
            for lead in 1:72
                f=fc_city[name][:,ci,lead]
                result=arrival_comparison(o,f)
                flags[1,lead]=result.success
                push!(lead_rows,(model=name,city=city,lead_hour=lead,
                    observation_arrival_bjt=result.oi>0 ? timestamp(event_times[result.oi]) : "",
                    forecast_arrival_bjt=result.fi>0 ? timestamp(event_times[result.fi]) : "",
                    absolute_error_hour=result.error_h,success=result.success,
                    finite_forecast_hours=count(isfinite,f),expected_hours=length(f)))
            end
            maxlead=only(Arrival.bool2leadtime(flags))
            arrival_max[name][ci]=maxlead
            complete=count(l->all(isfinite,fc_city[name][:,ci,l]),1:72)
            push!(city_rows,(model=name,city=city,station_count=length(city_groups[ci]),
                observed_event=events[ci],max_warning_lead_hour=maxlead,
                successful_lead_count=sum(flags),complete_lead_series_count=complete))
        end
    end
    for name in ("All-stations","Nation-only"), threshold in (6,12,24)
        flags=arrival_max[name].>=threshold
        for (denom,mask) in [("all_evaluated_cities",trues(length(cities))),("observed_dust_cities",events)]
            total=count(mask);number=count(flags.&mask)
            push!(fraction_rows,(model=name,threshold_hour=threshold,denominator=denom,
                qualifying_cities=number,total_cities=total,percent=total>0 ? 100number/total : NaN))
        end
    end
    peak_rows=NamedTuple[]
    observed_peak=peak_info(obs_city[:,beijing])
    observed_peak.index>0 && observed_peak.value>CONCENTRATION_THRESHOLD || error("No complete Beijing dust peak")
    peak_max=Dict("All-stations"=>0,"Nation-only"=>0)
    for lead in 1:72
        peaks=Dict(name=>peak_info(fc_city[name][:,beijing,lead]) for name in keys(lookups))
        for name in ("All-stations","Nation-only")
            p=peaks[name]
            success=p.index>0 && p.value>CONCENTRATION_THRESHOLD && abs(p.index-observed_peak.index)<=TOLERANCE_H
            success && (peak_max[name]=lead)
        end
        a,n=peaks["All-stations"],peaks["Nation-only"]
        ae=a.index>0 ? abs(a.index-observed_peak.index) : NaN
        ne=n.index>0 ? abs(n.index-observed_peak.index) : NaN
        push!(peak_rows,(lead_hour=lead,observation_peak_bjt=timestamp(event_times[observed_peak.index]),
            observation_peak=observed_peak.value,observation_peak_ties=observed_peak.ties,
            all_peak_bjt=a.index>0 ? timestamp(event_times[a.index]) : "",
            nation_peak_bjt=n.index>0 ? timestamp(event_times[n.index]) : "",
            all_peak=a.value,nation_peak=n.value,all_peak_ties=a.ties,nation_peak_ties=n.ties,
            all_absolute_error_hour=ae,nation_absolute_error_hour=ne,
            nation_minus_all_peak_time_hour=a.index>0 && n.index>0 ? n.index-a.index : NaN,
            all_peak_success=a.value>CONCENTRATION_THRESHOLD && ae<=TOLERANCE_H,
            nation_peak_success=n.value>CONCENTRATION_THRESHOLD && ne<=TOLERANCE_H))
    end
    beijing_summary=[
        (metric="corrected_arrival_max_successful_lead",all_stations_hour=arrival_max["All-stations"][beijing],
            nation_only_hour=arrival_max["Nation-only"][beijing],
            all_minus_nation_hour=arrival_max["All-stations"][beijing]-arrival_max["Nation-only"][beijing]),
        (metric="raw_peak_max_successful_lead",all_stations_hour=peak_max["All-stations"],
            nation_only_hour=peak_max["Nation-only"],
            all_minus_nation_hour=peak_max["All-stations"]-peak_max["Nation-only"])]
    series_rows=NamedTuple[]
    for lead in 1:72, (ti,t) in enumerate(event_times)
        push!(series_rows,(valid_time_bjt=timestamp(t),lead_hour=lead,
            initialization_bjt=timestamp(t-Hour(lead)),observation=obs_city[ti,beijing],
            all_stations=fc_city["All-stations"][ti,beijing,lead],
            nation_only=fc_city["Nation-only"][ti,beijing,lead]))
    end
    # Verify sources did not change while the run was in progress.
    bytes2hex(sha256(read(DETECTOR)))==detector_hash || error("Arrival code changed during calculation")
    for row in input_manifest
        st=stat(row.path)
        st.size==row.bytes && st.mtime==row.mtime_unix || error("Input changed during run: $(row.path)")
    end
    mkpath(output)
    for (filename,rows) in [
        ("input_inventory.csv",inventory),("input_files.csv",input_manifest),
        ("rmse_by_lead.csv",metric_rows),("rmse_first12_summary.csv",rmse_summary),
        ("city_warning_leads.csv",city_rows),("city_warning_details.csv",lead_rows),
        ("city_warning_fractions.csv",fraction_rows),
        ("beijing_peak_by_lead.csv",peak_rows),("beijing_summary.csv",beijing_summary),
        ("beijing_timeseries.csv",series_rows)]
        write_csv(joinpath(output,filename),rows)
    end
    write_csv(joinpath(output,"rmse_stations.csv"),[(station_id=id,) for id in station_ids])
    write_csv(joinpath(output,"city_stations.csv"),
        [(station_id=id,city=city_labels[i],type=at.types[city_ap[i]]) for (i,id) in enumerate(city_ids)])
    metadata=Dict("event_start_bjt"=>timestamp(event_start),"event_end_bjt"=>timestamp(event_end),
        "tolerance_hours"=>TOLERANCE_H,"concentration_threshold"=>CONCENTRATION_THRESHOLD,
        "derivative_threshold"=>GRADIENT_THRESHOLD,"observation_upper_exclusive"=>OBS_UPPER,
        "detector_path"=>DETECTOR,"detector_sha256"=>detector_hash,
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"julia_version"=>string(VERSION),
        "generated_bjt"=>string(now(BJT)),"rmse_leads"=>collect(1:12),
        "rmse_common_station_count"=>length(station_ids),"city_count"=>length(cities),
        "observed_dust_city_count"=>count(events),"city_station_count"=>length(city_ids),
        "ruf_first_initialization_bjt"=>timestamp(ruf_min),"ruf_last_initialization_bjt"=>timestamp(ruf_max))
    open(joinpath(output,"run_metadata.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    cp(DETECTOR,joinpath(output,"change_point_detect_used.jl"))
    open(joinpath(output,"summary.txt"),"w") do io
        println(io,"RUF 论文统计（修正候选点参数顺序；时间均为北京时间）")
        println(io,"预警事件窗口：",timestamp(event_start)," 至 ",timestamp(event_end),"（含端点）")
        println(io,"RMSE：公共国控站 ",length(station_ids)," 个；逐时效汇总站点—起报配对。")
        println(io,"前12h=1:12h，不含0h；平均降低百分比=12个时效的100×(Baseline RMSE-RUF RMSE)/Baseline RMSE的算术平均。")
        println(io,"independent_all_cycles沿用TOC统计：模型各用全部起报，样本不配对。")
        println(io,"matched_initialization_and_samples：三组模型相同起报、站点、有效时间和共同有效样本。")
        println(io,"这两项均是相同预报时效比较，不是U1与U24更新频率拼接比较。")
        for r in rmse_summary
            @printf(io,"%s | %s | 平均降低 %.4f%%\n",r.scheme,r.model,r.mean_hourly_reduction_percent)
        end
        println(io,"\n预警：原图5点中心平滑、归一化正增量>0.1、平滑序列最大值>200；候选点筛选参数顺序已修正。")
        println(io,"每次比较取两组候选点前min(3,观测候选数,预报候选数)个，按各自原始浓度最小值选到达点。")
        println(io,"允许±6h，边界包含；1:72h中最大成功时效，不要求较短时效连续成功。")
        println(io,"两组RUF在相同全站评价网络上计算城市平均；Nation-only表示仅同化国控，不表示只在国控站评价。")
        println(io,"城市分母同时报告全部评价城市与观测检出沙尘城市；没有成功匹配记0h。")
        for r in fraction_rows
            @printf(io,"%s | ≥%dh | %s | %d/%d = %.4f%%\n",
                r.model,r.threshold_hour,r.denominator,r.qualifying_cities,r.total_cities,r.percent)
        end
        println(io,"\n北京观测原始城市均值峰值时间：",timestamp(event_times[observed_peak.index]))
        println(io,"峰值统计独立于上升阶段的到达时间：原始城市序列最大值，峰值须>200，时间误差≤6h；并列最大值取首次。")
        for r in beijing_summary
            println(io,r.metric,": 全站=",r.all_stations_hour,"h；国控=",r.nation_only_hour,
                "h；全站增加=",r.all_minus_nation_hour,"h。")
        end
        println(io,"以上时效均来自固定lead、逐小时不同起报拼接的序列，与原图相同；不等同于一份单一起报预报的最早发布时刻。")
        println(io,"beijing_peak_by_lead.csv另列每个lead下两组峰值时间及其差值；不能将峰值时间平移等同于预警提前量。")
        println(io,"BaseFolder_EPICC若无RUF同期月份，则仅列入数据清单，不混入不同事件。")
        println(io,"未使用旧图缓存；未改观测、模型输出、PM10定义、单位或空间插值。")
    end
    println(read(joinpath(output,"summary.txt"),String))
    println("OUTPUT_DIR=",output)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main(ARGS)
