module QualityControlM
    using NaNStatistics#,MyPlotM
    include("locUtilities.jl")
    # include("convFilter.jl")
    include("checkCoarse.jl")
    include("checkTime.jl")
    include("checkSpace.jl")
    include("checkSpatiotemporal.jl")
    include("checkStagnant.jl")
    include("checkCalibration.jl")
    include("otherChecks.jl")
    include("QCraw.jl")

    export staObsQualityControl,backGroundEstimate


    """
    Set stdMutiplier to use a single probabilityLimit for all the checks. stdMutiplier for spatiotemporal check
    default  probabilityLimit:
    pdf_singleCoarse=>1e-7  pdf_st=>1e-7  pdf_stag=>1e-7  pdf_cali=>1e-7

    """
    function staObsQualityControl(staObs,staInfo,pollTypeS;probabilityLimit=1e-7,analysis=false,disp=false,ind2ref_all=axes(staObs,3),
        checkS=["obsRange","coarse","timeSpace","stagnant","calibration","PM25OutweightPM10","isolation"],
        )
        # time space quality control.option shoud include staInfo which is field with lon and lat DataFrame(
        # dimension of staObs should be (time,var,sta)
        # single probabilityLimit for all checks, adjust the stdRatio for checks and pollutants
        # typeName should be one of {'SO2'='SO2D',} {'NO2'='NO2D',}{"co"='COD',}{"o3"='O3H',}{'PM25'='PM25D',}{'PM10'='PM10D',}
        # 
        
        # updates:
        # 2016-09-21 original version
        
        # ****************** initialized *********************************
        if pollTypeS isa String;pollTypeS=[pollTypeS];end
        pollTypeS=[uppercase.(x) for x in pollTypeS]
        stalonlat=Matrix(staInfo[:,["lon","lat"]])
        # localizer=localizerM(stalonlat,stalonlat)
        staObsOut=copy(staObs)
        # ================== QC for observation ousite instrument range ===========================
        if "obsRange" in checkS
            if disp;    println("-------- obsRange -------- ")  ; end
            q_outRange=falses(size(staObs))
            for (n_poll,poll) in enumerate(pollTypeS)
                obsRange = getObsRange(poll)
                # upper and lower limit screen
                q_outRange[:,n_poll,:] = (staObs[:,n_poll,:] .< obsRange[1]) .| (staObs[:,n_poll,:] .> obsRange[2])
            end
            staObsOut[q_outRange] .= NaN
        end
        # ================== QC for observation with coarse error ==================================
        if "coarse" in checkS
            if disp;    println("--------  coarse -------- ");  end
            # do the coarse check for single coarse data ------------------------------
            q_singleCoarse=falses(size(staObsOut))
            # for n_poll in eachindex(pollTypeS)
            # Threads.@threads for n_poll in eachindex(pollTypeS)
            #     q_singleCoarse[:,n_poll,:]=checkSingleCoarse(view(staObsOut,:,n_poll,:);analysis=false)
            # end
            
            misfit_singleCoarse=similar(staObsOut)
            std_singleCoarse=similar(staObsOut)
            pdf_singleCoarse=similar(staObsOut)
            Threads.@threads for n_poll in eachindex(pollTypeS)
                q_singleCoarse[:,n_poll,:],misfit_singleCoarse[:,n_poll,:],std_singleCoarse[:,n_poll,:],pdf_singleCoarse[:,n_poll,:]=checkSingleCoarse(view(staObsOut,:,n_poll,:);analysis=true,probabilityLimit=1e-5)
            end
            staObsOut[q_singleCoarse] .= NaN
            # do the coarse check ------------------------------
            # q_daysCoarse,misfit_daysCoarse,std_daysCoarse,pdf_daysCoarse=mymapslices(x->
            # q_coarse,misfit_coarse,std_coarse,pdf_coarse=mymapslices(x->
            #     checkDaysCoarse(x;analysis=true)
            #     ,staObs
            # );

            # stdMutiplier=2
            # misfit_coarse=similar(staObsOut)
            # std_coarse=similar(staObsOut)
            # # set lower limit for std
            # Threads.@threads for n_poll in eachindex(pollTypeS)
            #     poll=pollTypeS[n_poll]
            #     misfit_coarse[:,n_poll,:]=medfilter(view(staObsOut,:,n_poll,:);edgeLen=24*5) .- view(staObsOut,:,n_poll,:);
            #     std_coarse[:,n_poll,:]=1.4826*medfilter(abs.(view(misfit_coarse,:,n_poll,:));edgeLen=24*5) .* stdMutiplier
            #     lowLimit_std=maximum(getObsRange(poll))/200
            #     replace!(x-> x < lowLimit_std ? lowLimit_std : x , view(std_coarse,:,n_poll,:))
            # end
            # # calculate pdf
            # pdf_coarse=pdf.(Normal(),-abs.(misfit_coarse./std_coarse));
            # pdf_coarse[isnan.(pdf_coarse)] .=1;
            # # get q_coarse
            # q_coarse = pdf_coarse .< probabilityLimit
            # # delete coarse obs
            # staObsOut[q_coarse ] .= NaN
        else
            # q_singleCoarse=falses(size(staObs))
            q_singleCoarse=falses(size(staObs))
        end
                
        # ================== QC for time and space continuit ==================================
        if "timeSpace" in checkS
            if disp;    println("-------- timeSpace -------- ");    end
            q_timeSpace,fit_time,fit_space,std_time,std_space,pdf_st=checkSpatiotemporal(staObsOut,stalonlat;analysis=true,ind2ref_all=ind2ref_all,disp=disp);
            misfit_space=fit_space-staObsOut
            misfit_time = fit_time - staObsOut
        else 
            q_timeSpace=falses(size(staObsOut))
        end

        # ================== QC for stagnant check ==========================
        if "stagnant" in checkS
            if disp;    println("-------- stagnant -------- "); end
            q_stag,misfit_stag,std_stag,pdf_stag=mymapslices((a,b,c)->
                checkStagnant(a,b,c;analysis=true,probabilityLimit=probabilityLimit),
                staObsOut,misfit_space,std_space
            );
        else
            q_stag=falses(size(staObsOut))
        end

        # ================== QC for calibration check ==========================    
        if "calibration" in checkS
            if disp;    println("-------- calibaration -------- "); end
            q_cali,misfit_cali,std_cali,pdf_cali=mymapslices(x->
                checkCalibration(x;analysis=true),
                staObsOut
            );
        else
            q_cali=falses(size(staObsOut))
        end
        # -------------- set staObsOut to NaN for those do not pass the stagant, spatiotemporal, and calibaration checks -------------
        staObsOut[q_stag .| q_timeSpace .| q_cali] .= NaN

        # ================== QC for PM2.5 > PM10 ==========================
        if "PM25OutweightPM10" in checkS
            if "PM25" in pollTypeS && "PM10" in pollTypeS
                if disp;    println("-------- PM25OutweightPM10-------- ");   end
                ind_pm25=indexin(["PM25"],pollTypeS)[1]
                ind_pm10=indexin(["PM10"],pollTypeS)[1]
                staObsOut_pm25,staObsOut_pm10=view(staObsOut,:,ind_pm25,:),view(staObsOut,:,ind_pm10,:)
                q_PM25OutweightPM10=staObsOut_pm25 .> staObsOut_pm10
                staObsOut_pm10[q_PM25OutweightPM10] .= NaN
            end
        end

        # ================== QC for isolation check ===================
        if "isolation" in checkS
            if disp;    println("-------- isolation -------- ");    end
            q_isolate=checkIsolation(staObsOut;edgeLen=12,validNum=4)
            # q_isolate=.!isnan(staObsOut) .& mapslice(x-> windowValidNumber(x;edgeLen=12),staObsOut,dims=1) .< 4
            # q_isolate=mapslices(x->checkIsolation(x;edgeLen=12,validNum=4),staObsOut,dims=1)
            staObsOut[q_isolate] .= NaN
        elseif analysis
            q_isolate=falses(size(staObs))
        end

        # ============================= output ================================
        # put analysis together
        if analysis 
            q_outlier = isnan.(staObsOut) .& .!isnan.(staObs)
            d_analysis =   Dict{String,Any}("q_outlier" => q_outlier)
            if "PM25OutweightPM10" in checkS
                if "PM25" in pollTypeS && "PM10" in pollTypeS
                    d_analysis["q_PM25OutweightPM10"] = q_PM25OutweightPM10
                end
            end
            if "obsRange" in checkS
                d_analysis["q_outRange"] = q_outRange
            end
            if "coarse" in checkS
                d_analysis["q_coarse"] = q_singleCoarse
                d_analysis["pdf_course"] = replace(pdf_singleCoarse , NaN => 1.0)

                # d_analysis["pdf_course"] = pdf_singleCoarse
            end
            
            if "timeSpace" in checkS
                d_analysis["q_timeSpace"] = q_timeSpace
                d_analysis["pdf_time"] = pdf.(Normal(),-abs.(misfit_time./std_time))
                d_analysis["pdf_space"] = pdf.(Normal(),-abs.(misfit_space./std_space))
                d_analysis["pdf_st"] = replace(pdf_st , NaN => 1.0)
                d_analysis["std_space"] = replace(std_space , NaN => 1.0)
                d_analysis["std_time"] = std_time
                d_analysis["estimate_time"] = staObs - misfit_time
                d_analysis["estimate_space"] = staObs - misfit_space
            end
            if "stagnant" in checkS
                d_analysis["q_stag"] = q_stag
                d_analysis["pdf_stag"] = replace(pdf_stag , NaN => 1.0)
            end
            if "calibration" in checkS
                d_analysis["q_cali"] = q_cali
                d_analysis["pdf_cali"] = replace(pdf_cali , NaN => 1.0)
            end
            if "isolation" in checkS
                d_analysis["q_isolate"] = q_isolate
            end
    
            return staObsOut, d_analysis
        else
            return staObsOut
        end
        
    end

    # function checkIsolation(inData::AbstractArray{T};edgeLen=12,validNum=4) where T
    #     # check observation with little valid observations around
    #     q_isolate = falses(size(inData))
    #     Threads.@threads   for n_var in axes(inData,2)
    #         validNum_series= Array{T}(undef,size(inData,1))
    #         for n_sta in axes(inData,3)
    #             validNum_series=movingValid!(validNum_series,view(inData,:,n_var,n_sta);edgeLen=edgeLen)
    #             q_isolate[:,n_var,n_sta] .= .!isnan.(view(inData,:,n_var,n_sta)) .& ( validNum_series.< validNum) 
    #         end
    #     end
    #     # q_isolate =  .!isnan.(inData) .& ( windowValidNumber(inData;edgeLen=edgeLen) .< validNum) 
    #     return q_isolate
    # end

    function backGroundEstimate(inData;edgeLen=24*15)
    # the output background has fields: quarter,median,threeQuarter,diffQuarter,diffDiffQuarter
    # dimensions of inData should be (time,var,sta)

        diffIn=diff(inData,dims=1)
        diffDiffIn=diff(diffIn,dims=1)
    
        background.median=mapslice(x-> medfilter(x;edgeLen=edgeLen),inData,dims=1)
        background.diffMedian=mapslice(x-> medfilter(x;edgeLen=edgeLen),diffIn,dims=1)
        background.diffDiffMedian=mapslice(x-> medfilter(x;edgeLen=edgeLen),diffDiffIn,dims=1)
        
        # steps and unique obs used for checks
        minPositive=mapslice(y-> minimum(filter(x->x>0,y)),inData,dims=1)
        background.minPositive=minPositive
        background.step=nanpctile(abs.(diff(inData,dims=1)),3,dims=1)
        return background
    end

    
    function QCserie(staObs,polType;probabilityLimit=1e-7,analysis=false,
        checkS=["obsRange","coarse","time","stagnant","calibration","isolation"],
        )
        staObsOut=copy(staObs)

        if "obsRange" in checkS
            q_outRange=flases(size(staObs))
            obsRange = getObsRange(poll)
            # upper and lower limit screen
            q_outRange= (staObs .< obsRange[1]) .| (staObs .> obsRange[2])
            staObsOut[q_outRange] .= NaN
        end
        
        if "coarse" in checkS
            q_coarse,misfit_coarse,std_coarse,pdf_coarse=checkDaysCoarse(inData;analysis=true);
            staObsOut[q_coarse ] .= NaN
        else
            q_coarse=falses(size(staObsOut))
        end

        if "time" in checkS
            q_time,misfit_time,std_time,pdf_time=checkTemporal(staObsOut;analysis=true);
        else 
            q_time=falses(size(staObsOut))
        end

        if "stagnant" in checkS
            q_lowDiff,q_lowDiffDiff=checkStagnant(inData)
            q_stag = q_lowDiff .| q_lowDiffDiff
        else
            q_stag=falses(size(staObsOut))
        end

        if "calibration" in checkS
            q_cali = calibrationCheck(inData;probabilityLimit=probabilityLimit)
        else
            q_cali=falses(size(staObsOut))
        end
        # -------------- set staObsOut to NaN for those do not pass the stagant, spatiotemporal, and calibaration checks -------------
        staObsOut[q_stag .| q_timeSpace .| q_cali] .= NaN

        if "isolation" in checkS
            q_isolate=checkIsolation(staObsOut;edgeLen=12,validNum=4)
            staObsOut[q_isolate] .= NaN
        elseif analysis
            q_isolate=falses(size(staObs))
        end
    
        
        if analysis 
            d_analysis = Dict("q_outlier" => isnan.(staObsOut) .& .!isnan.(staObs),   "q_coarse" => q_coarse,   "q_time" => q_time,  
                "q_lowDiff" => q_lowDiff,      "q_lowDiffDiff" => q_lowDiffDiff,    "q_stag" => q_stag,       "q_cali" => q_cali,
                "q_isolate" => q_isolate,       "misfit_time" => misfit_time,     "std_time" => std_time,
            )

            return staObsOut, d_analysis
        else
            return staObsOut
        end
    end
end # module QualityControlM
