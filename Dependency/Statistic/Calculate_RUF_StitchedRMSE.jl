#!/usr/bin/env julia
# julia --startup-file=no --project=/home/zsm/IETM/code/RUF_Module \
#   /home/zsm/IETM/code/PaperFigure_code/Statistic/Calculate_RUF_StitchedRMSE.jl
# Supersedes the independent-cycle RMSE comparison for the requested continuous
# 2025-04-09..2025-04-21 verification period. Arrival statistics are not rerun.
include(joinpath(@__DIR__, "Calculate_RUF_PaperStatistics.jl"))
using Distributed

const RMSE_START_BJT = DateTime(2025, 4, 9, 0)
const RMSE_END_BJT = DateTime(2025, 4, 21, 23)
const RMSE_LEADS = collect(1:12)
const JULIA_PROJECT = "/home/zsm/IETM/code/RUF_Module"

# Strictly earlier initialization: at an update boundary use the previous
# daily cycle's lead 24, not the new cycle's lead 0.
function select_daily_cycle(t, sorted_initializations)
    i = searchsortedfirst(sorted_initializations, t) - 1
    i > 0 || error("No baseline initialization strictly before $t")
    lead = Int(Dates.value(t - sorted_initializations[i]) ÷ 3_600_000)
    1 <= lead <= 24 || error("Daily baseline coverage gap at $t (lead=$lead)")
    sorted_initializations[i] + Hour(lead) == t || error("Non-hourly time")
    return i, lead
end

function stitched_tests()
    @testset "Continuous-time stitching" begin
        cycles = [DateTime(2025,4,8,20), DateTime(2025,4,9,20)]
        @test select_daily_cycle(DateTime(2025,4,9,0),cycles) == (1,4)
        @test select_daily_cycle(DateTime(2025,4,9,20),cycles) == (1,24)
        @test select_daily_cycle(DateTime(2025,4,9,21),cycles) == (2,1)
        @test_throws ErrorException select_daily_cycle(DateTime(2025,4,8,20),cycles)
        @test_throws ErrorException select_daily_cycle(DateTime(2025,4,11,0),cycles)
        @test RMSE_START_BJT-Hour(1) == DateTime(2025,4,8,23)
        @test length(RMSE_START_BJT:Hour(1):RMSE_END_BJT) == 312
        @test ratio_drop(sqrt(mean([9.0,16.0])),sqrt(mean([1.0,4.0]))) ≈ 55.278640450004204
        @test mean([ratio_drop(10.,5.),ratio_drop(10.,8.)]) == 35.
    end
end

function input_record(model, path, init)
    s = stat(path)
    return (model=model,path=path,initialization_bjt=timestamp(init),
        bytes=s.size,mtime_unix=s.mtime)
end

# One worker per experiment; each MAT file is read once, even though several
# leads from it contribute to the continuous verification window.
function load_stitched_ruf(name, cycles, ids, station_indices, times)
    data = fill(NaN,length(station_indices),length(times),length(RMSE_LEADS))
    lookup = Dict(t=>i for (i,t) in enumerate(times))
    filled = falses(length(times),length(RMSE_LEADS))
    records = NamedTuple[]
    for (j,c) in enumerate(cycles)
        push!(records,input_record(name,c.path,c.initialization))
        a = read_cycle(c,ids)
        lp = Dict(l=>i for (i,l) in enumerate(a.leads))
        for (li,lead) in enumerate(RMSE_LEADS)
            ti = get(lookup,c.initialization+Hour(lead),0)
            ti == 0 && continue
            haskey(lp,lead) || error("Missing lead $lead in $(c.path)")
            !filled[ti,li] || error("Duplicate RUF valid-time/lead")
            data[:,ti,li] .= a.data[1,a.variable,station_indices,lp[lead]]
            filled[ti,li] = true
        end
        if j==1 || j%48==0 || j==length(cycles)
            println("$name: $j/$(length(cycles)) files read");flush(stdout);GC.gc(false)
        end
    end
    all(filled) || error("Incomplete stitched RUF timeline: $name")
    return (data=data,records=records)
end

function stitched_main(args)
    args == ["--self-test"] && return stitched_tests()
    stitched_tests()
    output = if isempty(args)
        joinpath(@__DIR__,"output","stitched_"*Dates.format(now(BJT),dateformat"yyyymmdd_HHMMSS"))
    elseif length(args)==2 && args[1]=="--output-dir"
        abspath(args[2])
    else
        error("Usage: Calculate_RUF_StitchedRMSE.jl [--output-dir PATH | --self-test]")
    end
    ispath(output) && error("Refusing to overwrite output: $output")
    times = collect(RMSE_START_BJT:Hour(1):RMSE_END_BJT)
    source_hash = bytes2hex(sha256(read(@__FILE__)))
    helper_path = joinpath(@__DIR__,"Calculate_RUF_PaperStatistics.jl")
    helper_hash = bytes2hex(sha256(read(helper_path)))
    inventories = NamedTuple[]
    base_all = NamedTuple[]
    for root in BASE_ROOTS
        append!(base_all,discover(root;utc=true))
    end
    sort!(base_all;by=c->c.initialization)
    length(unique(c.initialization for c in base_all))==length(base_all) ||
        error("Duplicate baseline initializations across roots")
    base_inits = [c.initialization for c in base_all]
    selections = [select_daily_cycle(t,base_inits) for t in times]
    base_indices = sort(unique(first.(selections)))
    used_base = base_all[base_indices]
    all(diff([c.initialization for c in used_base]).==Hour(24)) ||
        error("Baseline cycles must be daily")
    cycle_maps = Dict(name=>Dict(c.initialization=>c for c in discover(root))
        for (name,root) in [("All-stations",ALL_ROOT),("Nation-only",NATION_ROOT)])
    required_inits = sort(unique(t-Hour(l) for t in times for l in RMSE_LEADS))
    selected_ruf = Dict{String,Vector{NamedTuple}}()
    for name in ("All-stations","Nation-only")
        all(haskey(cycle_maps[name],t) for t in required_inits) ||
            error("Missing required hourly RUF initializations: $name")
        selected_ruf[name] = [cycle_maps[name][t] for t in required_inits]
    end
    for root in BASE_ROOTS
        found=filter(c->c.root==root,base_all);used=filter(c->c.root==root,used_base)
        push!(inventories,(model="Baseline",root=root,available_cycles=length(found),
            used_cycles=length(used),selection="daily leads 1:24 covering the fixed valid-time window"))
    end
    for (name,root) in [("All-stations",ALL_ROOT),("Nation-only",NATION_ROOT)]
        push!(inventories,(model=name,root=root,available_cycles=length(cycle_maps[name]),
            used_cycles=length(selected_ruf[name]),selection="initialization = valid time - fixed lead"))
    end
    println("Loading observations and common-station metadata");flush(stdout)
    records=[input_record("Observation",OBS_DEFAULT,nothing)]
    obs = read_station(OBS_DEFAULT)
    bt=read_station(first(used_base).path)
    at=read_station(first(selected_ruf["All-stations"]).path)
    nt=read_station(first(selected_ruf["Nation-only"]).path)
    all(==("B"),bt.types) || error("Expected national-control baseline stations")
    available=intersect(Set(obs.ids),Set(at.ids),Set(nt.ids))
    stations=[id for id in bt.ids if id in available]
    isempty(stations) && error("No common evaluation stations")
    op=positions(obs.ids,stations,"observations")
    bp=positions(bt.ids,stations,"baseline")
    model_ids=Dict("All-stations"=>at.ids,"Nation-only"=>nt.ids)
    model_indices=Dict(name=>positions(ids,stations,name) for (name,ids) in model_ids)
    obs_times=Dict(t=>i for (i,t) in enumerate(obs.times))
    length(obs_times)==length(obs.times) || error("Duplicate observation times")
    all(haskey(obs_times,t) for t in times) || error("Incomplete observation timeline")
    obs.leads==[0] || error("Expected lead-0 observation data")
    o=permutedims(Float64.(obs.data[[obs_times[t] for t in times],obs.variable,op,1]),(2,1))
    baseline=fill(NaN,size(o))
    mapping=NamedTuple[]
    for gi in base_indices
        c=base_all[gi];push!(records,input_record("Baseline",c.path,c.initialization))
        a=read_cycle(c,bt.ids);lp=Dict(l=>i for (i,l) in enumerate(a.leads))
        for (ti,(chosen,lead)) in enumerate(selections)
            chosen==gi || continue
            haskey(lp,lead) || error("Missing baseline lead $lead")
            baseline[:,ti] .= a.data[1,a.variable,bp,lp[lead]]
        end
    end
    for (ti,t) in enumerate(times)
        chosen,lead=selections[ti];c=base_all[chosen]
        push!(mapping,(valid_time_bjt=timestamp(t),baseline_initialization_bjt=timestamp(c.initialization),
            baseline_lead_hour=lead,baseline_path=c.path))
    end
    println("Stitched baseline: $(length(times)) hours, $(length(used_base)) daily cycles, $(length(stations)) stations")
    obs=nothing;at=nothing;nt=nothing;bt=nothing;GC.gc()
    # Separate Julia processes avoid sharing MAT/HDF5 handles across threads.
    pids=addprocs(2;exeflags=`--startup-file=no --project=$JULIA_PROJECT --threads=1`)
    results=Dict{String,Any}()
    try
        for pid in pids
            remotecall_wait(Base.include,pid,Main,abspath(@__FILE__))
        end
        @sync for (pid,name) in zip(pids,("All-stations","Nation-only"))
            @async results[name]=remotecall_fetch(load_stitched_ruf,pid,name,
                selected_ruf[name],model_ids[name],model_indices[name],times)
        end
    finally
        rmprocs(pids)
    end
    common=validobs.(o).&validfc.(baseline)
    obs_valid=count(validobs,o)
    obs_base_valid=count(common)
    for name in ("All-stations","Nation-only")
        append!(records,results[name].records)
        for li in eachindex(RMSE_LEADS)
            common .&= validfc.(@view results[name].data[:,:,li])
        end
    end
    n=count(common);n>0 || error("No common finite samples")
    baseline_sse=sum((baseline[common].-o[common]).^2)
    baseline_rmse=sqrt(baseline_sse/n)
    println("Shared station-time samples: $n; baseline RMSE: $baseline_rmse");flush(stdout)
    rows=NamedTuple[];hourly=NamedTuple[];summaries=NamedTuple[];samples=NamedTuple[]
    for name in ("All-stations","Nation-only")
        reductions=Float64[];rmses=Float64[];total_sse=0.0
        for (li,lead) in enumerate(RMSE_LEADS)
            r=@view results[name].data[:,:,li]
            sse=sum((r[common].-o[common]).^2);value=sqrt(sse/n)
            reduction=ratio_drop(baseline_rmse,value)
            push!(reductions,reduction);push!(rmses,value);total_sse+=sse
            push!(rows,(model=name,lead_hour=lead,baseline_rmse=baseline_rmse,ruf_rmse=value,
                reduction_percent=reduction,paired_values=n,baseline_sse=baseline_sse,ruf_sse=sse,
                valid_time_count=length(times),station_count=length(stations)))
            for (ti,t) in enumerate(times)
                mask=common[:,ti]
                bs=sum((baseline[mask,ti].-o[mask,ti]).^2)
                rs=sum((r[mask,ti].-o[mask,ti]).^2)
                push!(hourly,(model=name,lead_hour=lead,valid_time_bjt=timestamp(t),
                    ruf_initialization_bjt=timestamp(t-Hour(lead)),
                    baseline_initialization_bjt=mapping[ti].baseline_initialization_bjt,
                    baseline_lead_hour=mapping[ti].baseline_lead_hour,
                    paired_values=count(mask),baseline_sse=bs,ruf_sse=rs))
                # One explicitly identified station provides a compact numerical trace.
                push!(samples,(model=name,lead_hour=lead,valid_time_bjt=timestamp(t),
                    station_id=stations[1],observation=o[1,ti],baseline=baseline[1,ti],
                    ruf=r[1,ti],included_in_common_mask=common[1,ti]))
            end
        end
        mean_drop=mean(reductions)
        @assert isapprox(mean_drop,ratio_drop(baseline_rmse,mean(rmses));atol=1e-12)
        push!(summaries,(model=name,first_lead=1,last_lead=12,
            mean_hourly_reduction_percent=mean_drop,baseline_rmse=baseline_rmse,
            mean_ruf_rmse=mean(rmses),
            pooled_12h_rmse=sqrt(total_sse/(12n)),
            pooled_12h_reduction_percent=ratio_drop(baseline_rmse,sqrt(total_sse/(12n))),
            paired_values_per_lead=n))
    end
    bytes2hex(sha256(read(@__FILE__)))==source_hash || error("Script changed during run")
    bytes2hex(sha256(read(helper_path)))==helper_hash || error("IO helper changed during run")
    for r in records
        st=stat(r.path)
        st.size==r.bytes && st.mtime==r.mtime_unix || error("Input changed during run: $(r.path)")
    end
    mkpath(output)
    for (file,data) in [("rmse_by_lead.csv",rows),("rmse_first12_summary.csv",summaries),
        ("baseline_stitch_map.csv",mapping),("hourly_error_sums.csv",hourly),
        ("trace_station_series.csv",samples),("input_inventory.csv",inventories),("input_files.csv",records)]
        write_csv(joinpath(output,file),data)
    end
    write_csv(joinpath(output,"evaluation_stations.csv"),[(station_id=id,) for id in stations])
    metadata=Dict("valid_start_bjt"=>timestamp(first(times)),"valid_end_bjt"=>timestamp(last(times)),
        "valid_hours"=>length(times),"baseline_leads"=>"1:24; initialization strictly before valid time",
        "ruf_leads"=>RMSE_LEADS,"station_count"=>length(stations),"possible_station_time_pairs"=>length(o),
        "valid_observation_pairs"=>obs_valid,"valid_obs_baseline_pairs"=>obs_base_valid,
        "shared_pairs_all_models_and_leads"=>n,"baseline_cycles_used"=>length(used_base),
        "ruf_cycles_used_per_model"=>length(required_inits),"observation_file"=>OBS_DEFAULT,
        "observation_validity"=>"finite, 0 <= PM10 < 2000","forecast_validity"=>"finite, PM10 >= 0",
        "evaluation"=>"pooled common station-time samples; same mask for both RUF models and all 12 leads",
        "script_sha256"=>source_hash,"helper_sha256"=>helper_hash,"julia_version"=>string(VERSION),
        "generated_bjt"=>string(now(BJT)))
    open(joinpath(output,"run_metadata.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    open(joinpath(output,"summary.txt"),"w") do io
        println(io,"连续有效时间序列 RMSE（替代上一版各自起报统计口径）")
        println(io,"时间：北京时间 ",timestamp(first(times))," 至 ",timestamp(last(times)),"，共 ",length(times)," h。")
        println(io,"基础预报：逐日取1—24h拼接，0h不参与；共 ",length(used_base)," 个批次。")
        println(io,"RUF：每个固定L=1:12h，以t-L起报的第L小时预报拼接；每组读取 ",length(required_inits)," 个批次。")
        println(io,"第1h示例：2025-04-08 23:00起报 → 2025-04-09 00:00有效时间。")
        println(io,"评价站点：共有的 ",length(stations)," 个国控站。两组模型和全部时效使用相同的 ",n," 个有效站点—时间样本。")
        println(io,"所有RMSE均相对于观测计算；没有将RUF与基础预报的差值当作RMSE。")
        @printf(io,"基础预报RMSE固定为 %.8f μg/m³。\n",baseline_rmse)
        println(io,"平均降幅 = mean_L[100×(RMSE_base-RMSE_RUF,L)/RMSE_base]；与12个时效合并求RMSE的降幅不同。")
        for r in summaries
            @printf(io,"%s：前12h平均RMSE %.8f μg/m³，平均降低 %.6f%%。\n",
                r.model,r.mean_ruf_rmse,r.mean_hourly_reduction_percent)
        end
        println(io,"BaseFolder_EPICC无覆盖该窗口的批次，未使用；有效基础数据来自paper_data/baseFolder。")
        println(io,"到达时间和峰值预警统计未重算。本次只更新连续时间序列RMSE；原始数据未修改。")
    end
    println(read(joinpath(output,"summary.txt"),String))
    println("OUTPUT_DIR=",output)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && stitched_main(ARGS)
