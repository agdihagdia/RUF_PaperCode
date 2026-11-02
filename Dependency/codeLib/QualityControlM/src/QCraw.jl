using Distributions
export QCraw,QCraw_v2

function setStdLowerLimit(stdData,stdLowLimit_all)    
    for (n_poll,lowerLimit_std) in enumerate(stdLowLimit_all) # set lower limit for std
        std_curr=view(stdData,:,n_poll,:)
        std_curr[std_curr .< lowerLimit_std] .= lowerLimit_std
    end
    return stdData
end
function lowerStd_forNegativeMisfit(stdData,misfit;adjustRatio=0.6)
    # set lower limit for outliers with lower value
    q=misfit .< 0
    stdData[q] .= stdData[q] .* adjustRatio
    return stdData
end
"""
Set stdMutiplier to use a single probabilityLimit for all the checks. stdMutiplier for spatiotemporal check
default  probabilityLimit:
pdf_singleCoarse=>1e-7  pdf_st=>1e-7  pdf_stag=>1e-7  pdf_cali=>1e-7

"""
function QCraw(staObs,staInfo,pollTypeS;probabilityLimit=1e-10,basePollName=nothing,analysis=false,disp=false,ind2ref_all=axes(staObs,3)  )
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

    stdLowLimit_all=[isequal(poll,"CO") ? 0.1 : 3 for poll in pollTypeS]

    if isnothing(basePollName)
        basePollName=copy(pollTypeS)
        for (n_poll,poll) in enumerate(pollTypeS)
            if poll in ["PM10","SO2","NO2","CO"]
                basePollName[n_poll]="PM25"
            elseif poll in ["PM25"]
                basePollName[n_poll]="PM10"
            elseif poll in ["O3"]
                basePollName[n_poll]="NO2"
            else
                error("Poll type not supported: $poll")
            end
        end
    end

    # ================== Elimilate for Good ===========================
    # check out of instruments range-----------
    q_outRange=falses(size(staObs))
    for (n_poll,poll) in enumerate(pollTypeS)
        obsRange = getObsRange(poll)
        # upper and lower limit screen
        q_outRange[:,n_poll,:] = (staObs[:,n_poll,:] .< obsRange[1]) .| (staObs[:,n_poll,:] .> obsRange[2])
    end
    staObsOut[q_outRange] .= NaN
    
    # course ===========================
    # calculate fit_singleCoarse, and misfit_singleCoarse,diff_neighbour
    fit_singleCoarse=similar(staObsOut)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        fit_singleCoarse[:,n_poll,:]=medfilter(view(staObsOut,:,n_poll,:);edgeLen=1)
    end
    misfit_singleCoarse =  staObsOut - fit_singleCoarse
    # get stdBase_curr
    diffData=abs.(diff(staObs;dims=1));diffData=cat(dims=1,diffData[1:1,:,:],diffData)
    stdBase_curr=similar(diffData)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        stdBase_curr[:,n_poll,:]=1.4826*medfilter(view(diffData,:,n_poll,:);edgeLen=12*3)
    end    
    stdBase_curr=setStdLowerLimit(stdBase_curr,stdLowLimit_all)   

    # get diff_neighbour
    diff_neighbour=fill(NaN,size(staObs))
    diff_neighbour[2:end-1,:,:]=abs.(staObs[3:end,:,:]-staObs[1:end-2,:,:]);
    
    # elimilate by singleCoarse ------------
    stdMultipler=5;
    pdf_singleCoarse=pdf.(Normal(),-abs.(misfit_singleCoarse./(stdMultipler .* stdBase_curr)));
    q_singleCoarse=.!isnan.(staObs) .& ( pdf_singleCoarse .< probabilityLimit ) .& (abs.(misfit_singleCoarse) .> diff_neighbour)
    staObsOut[q_singleCoarse] .= NaN

    # q_bigErr= q_outRange .| q_singleCoarse
    # ================== Elimilate Together ===========================
    pollBase=staObsOut[:,indexin(basePollName,pollTypeS),:]
    pollRatio=staObsOut ./ pollBase
    # check spatial ===============
    q_space,fit_space,std_space,pdf_space=checkSpatial(staObsOut,stalonlat;analysis=true,ind2ref_all=ind2ref_all,probabilityLimit=probabilityLimit,
        edgeLen_misfit=24*15,edgeLen_std_space=24*3, numRef = 10) # varing cutoff that ensure 10 sites to reference
    misfit_space= staObsOut - fit_space
    std_space=lowerStd_forNegativeMisfit(std_space,misfit_space;adjustRatio=0.7)
    pdf_space=pdf.(Normal(),-abs.(misfit_space./std_space));
    q_space=  pdf_space .< probabilityLimit 
    # fit_space=misfit_space .+ staObs
    # fill the NaN with fit_space ------------
    staObsOut_filled=copy(staObsOut)
    q_fill=isnan.(staObsOut_filled)
    staObsOut_filled[q_fill] = fit_space[q_fill]
    # check single course, but wait for elimilate later -----------    
    Threads.@threads for n_poll in eachindex(pollTypeS)
        fit_singleCoarse[:,n_poll,:]=medfilter(view(staObsOut_filled,:,n_poll,:);edgeLen=1)
    end
    misfit_singleCoarse = staObsOut_filled - fit_singleCoarse
    diff_neighbour[2:end-1,:,:]=abs.(staObsOut_filled[3:end,:,:]-staObsOut_filled[1:end-2,:,:]);    
    stdMultipler=1.2;
   
    diffData=abs.(diff(staObsOut_filled;dims=1));diffData=cat(dims=1,diffData[1:1,:,:],diffData)
    std_course=similar(diffData)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        std_course[:,n_poll,:]=stdMultipler .* 1.4826*medfilter(view(diffData,:,n_poll,:);edgeLen=12*3)
    end        
    std_course=setStdLowerLimit(std_course,stdLowLimit_all)      
    std_course=lowerStd_forNegativeMisfit(std_course,misfit_singleCoarse;adjustRatio=0.8)
    pdf_singleCoarse=pdf.(Normal(),-abs.(misfit_singleCoarse./std_course));
    q_singleCoarse=.!isnan.(staObs) .& ( pdf_singleCoarse .< probabilityLimit ) .& (abs.(misfit_singleCoarse) .> diff_neighbour)

    # check temporal ===============
    fit_time=nanFilter(staObsOut_filled; responsetype = Lowpass(1/(60*60*8); fs=1/60/60)) 
    misfit_time= staObsOut - fit_time
    std_time=1.3 .* movingMean(abs.(misfit_time);edgeLen=12*3)
    # movingStd(abs.(fit_time - staObsOut);edgeLen=12*3) 
    std_time=setStdLowerLimit(std_time,stdLowLimit_all)  
    std_time=lowerStd_forNegativeMisfit(std_time,misfit_time;adjustRatio=0.7)
    pdf_time=pdf.(Normal(),-abs.(misfit_time./std_time));
    q_time= pdf_time .< probabilityLimit
        
    # check pollRatio ===============    
    fit_pollRatio=similar(pollRatio)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        fit_pollRatio[:,n_poll,:]=medfilter(view(pollRatio,:,n_poll,:);edgeLen=12*3)
    end
    misfit_pollRatio= pollRatio - fit_pollRatio
    std_ratioMisfit=similar(pollRatio)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        std_ratioMisfit[:,n_poll,:]=1.4826*medfilter(abs.(view(misfit_pollRatio,:,n_poll,:));edgeLen=12*3)
    end      
    std_ratioMisfit[1:3,:,:] .= std_ratioMisfit[4:4,:,:] # deal with boundary
    std_ratioMisfit[end-3:end,:,:] .= std_ratioMisfit[end-4:end-4,:,:]# deal with boundary
    std_ratioMisfit=lowerStd_forNegativeMisfit(std_ratioMisfit,misfit_pollRatio;adjustRatio=0.6)
    pdf_pollRatio=pdf.(Normal(),-abs.(misfit_pollRatio./(1.0 .* std_ratioMisfit)));
    q_pollRatio= pdf_pollRatio .< probabilityLimit 

    # std_ratioMisfit=movingStd(abs.(fit_pollRatio - pollRatio);edgeLen=12*7) 
    # q_pollRatio,fit_pollRatio,std_ratioMisfit,pdf_pollRatio=mymapslices((x,fit,std)->
    #     checkTemporal(x;analysis=true,edgeLen_std=12*3,probabilityLimit=probabilityLimit,fit_curr=fit,std_curr=std,lowLimit_std=0.0)
    #     ,pollRatio,fit_pollRatio,std_ratioMisfit);

    # QC for stagnant check ==========================
    q_stag,misfit_stag,std_stag,pdf_stag=mymapslices((a,b,c)->
            checkStagnant(a,b,c;analysis=true,probabilityLimit=probabilityLimit),
            staObsOut,misfit_space,std_space
        );
    
    # -------------- set staObsOut to NaN for those do not pass the stagant, spatiotemporal, and calibaration checks -------------
    # staObsOut[q_stag .| (q_time .& q_space) .| (q_singleCoarse .& q_space) .| (pdf_space .< 1e-10 *  probabilityLimit)] .= NaN
    q_combine= q_stag .| (q_time .& q_space) .| (q_singleCoarse .& q_space) .| 
            (pdf_space .< 1e-10 *  probabilityLimit) .|  (q_pollRatio .& q_space)
    # special treatment for selected pollutants
    for (n_poll,poll) in enumerate(pollTypeS)
        if poll in ["SO2","NO2"]
            q_combine[:,n_poll,:].=(pdf_space[:,n_poll,:] .< 1e-10 *  probabilityLimit) .& (pdf_time[:,n_poll,:] .< 1e-10 *  probabilityLimit) .& (pdf_singleCoarse[:,n_poll,:] .< 1e-10 *  probabilityLimit) 
        end
    end
    staObsOut[q_combine] .= NaN
    
    # ================== Elimilate after ===========================
    # QC for isolation check ===================
    edgeLen=3
    num_noValid=getOutlierNum(staObsOut;edgeLen=edgeLen)
    # stdDecreaseRatio= max.(0.0,num_noValid .- edgeLen .+ 2) ./ (edgeLen .+ 3)
    stdDecreaseRatio= 0.7 .* max.(0.0,num_noValid ) ./ (2*edgeLen .+ 1)
    std_isolation=std_space .* (1 .- stdDecreaseRatio)
    pdf_isolate=pdf.(Normal(),-abs.(misfit_space./std_isolation));    
    q_isolate= pdf_isolate .< probabilityLimit
    staObsOut[q_isolate] .= NaN

    # ============================= output ================================
    if analysis 
        q_outlier= isnan.(staObsOut) .& .!isnan.(staObs)
        d_analysis =   Dict{String,Any}("q_outlier" => q_outlier,
            "q_outRange" => q_outRange,
            "q_singleCoarse" => q_singleCoarse,
            "q_time" => q_time,
            "q_space" => q_space,
            "q_stag" => q_stag,
            "q_isolate" => q_isolate,
            "misfit_singleCoarse" => misfit_singleCoarse,
            "fit_time" => fit_time,
            "fit_space" => fit_space,
            "pollRatio" => pollRatio,
            "fit_pollRatio" => fit_pollRatio,
            "misfit_stag" => misfit_stag,
            "std_time" => std_time,
            "std_space" => std_space,
            "std_stag" => std_stag,
            "std_ratioMisfit" => std_ratioMisfit,
            "pdf_singleCoarse" => pdf_singleCoarse,
            "pdf_time" => pdf_time,
            "pdf_space" => pdf_space,
            "pdf_stag" => pdf_stag,
            "pdf_pollRatio" => pdf_pollRatio,
            "pdf_isolate" => pdf_isolate,
            "fit_singleCoarse" => fit_singleCoarse
        )

        return staObsOut,d_analysis
    else
        return staObsOut
    end
end


function QCraw_v2(staObs,staInfo,pollTypeS;probabilityLimit=1e-10,basePollName=nothing,analysis=false,disp=false,ind2ref_all=axes(staObs,3)  )
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

    stdLowLimit_all=[isequal(poll,"CO") ? 0.1 : 3 for poll in pollTypeS]

    if isnothing(basePollName)
        basePollName=copy(pollTypeS)
        for (n_poll,poll) in enumerate(pollTypeS)
            if poll in ["PM10","SO2","NO2","CO"]
                basePollName[n_poll]="PM25"
            elseif poll in ["PM25"]
                basePollName[n_poll]="PM10"
            elseif poll in ["O3"]
                basePollName[n_poll]="NO2"
            else
                error("Poll type not supported: $poll")
            end
        end
    end

    # ================== Elimilate for Good ===========================
    # check out of instruments range-----------
    q_outRange=falses(size(staObs))
    for (n_poll,poll) in enumerate(pollTypeS)
        obsRange = getObsRange(poll)
        # upper and lower limit screen
        q_outRange[:,n_poll,:] = (staObs[:,n_poll,:] .< obsRange[1]) .| (staObs[:,n_poll,:] .> obsRange[2])
    end
    staObsOut[q_outRange] .= NaN
    
    # course ===========================
    # calculate fit_singleCoarse, and misfit_singleCoarse,diff_neighbour
    fit_singleCoarse=similar(staObsOut)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        fit_singleCoarse[:,n_poll,:]=medfilter(view(staObsOut,:,n_poll,:);edgeLen=1)
    end
    misfit_singleCoarse =  staObsOut - fit_singleCoarse
    # get stdBase_curr
    diffData=abs.(diff(staObs;dims=1));diffData=cat(dims=1,diffData[1:1,:,:],diffData)
    stdBase_curr=similar(diffData)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        stdBase_curr[:,n_poll,:]=1.4826*medfilter(view(diffData,:,n_poll,:);edgeLen=12*3)
    end    
    stdBase_curr=setStdLowerLimit(stdBase_curr,stdLowLimit_all)   

    # get diff_neighbour
    diff_neighbour=fill(NaN,size(staObs))
    diff_neighbour[2:end-1,:,:]=abs.(staObs[3:end,:,:]-staObs[1:end-2,:,:]);
    
    # elimilate by singleCoarse ------------
    stdMultipler=5;
    pdf_singleCoarse=pdf.(Normal(),-abs.(misfit_singleCoarse./(stdMultipler .* stdBase_curr)));
    q_singleCoarse=.!isnan.(staObs) .& ( pdf_singleCoarse .< probabilityLimit ) .& (abs.(misfit_singleCoarse) .> diff_neighbour)

    for (n_poll,poll) in enumerate(pollTypeS)
        if isequal(poll,"SO2") # do not perform single coarse check for so2
            q_singleCoarse[:,n_poll,:].=false
        end
    end
    staObsOut[q_singleCoarse] .= NaN

    # q_bigErr= q_outRange .| q_singleCoarse
    # ================== Elimilate Together ===========================
    # check spatial ===============
    q_space,fit_space,std_space,pdf_space=checkSpatial(staObsOut,stalonlat;analysis=true,ind2ref_all=ind2ref_all,probabilityLimit=probabilityLimit,
        edgeLen_misfit=24*15,edgeLen_std_space=24*3, numRef = 10) # varing cutoff that ensure 10 sites to reference
        # edgeLen_misfit=24*15,edgeLen_std_space=24*3, numRef = 10) # varing cutoff that ensure 10 sites to reference
    misfit_space= staObsOut - fit_space
    std_space=lowerStd_forNegativeMisfit(std_space,misfit_space;adjustRatio=0.5)
    pdf_space=pdf.(Normal(),-abs.(misfit_space./std_space));
    q_space=  pdf_space .< probabilityLimit 
    # fit_space=misfit_space .+ staObs
    # fill the NaN with fit_space ------------
    staObsOut_filled=copy(staObsOut)
    q_fill=isnan.(staObsOut_filled)
    staObsOut_filled[q_fill] = fit_space[q_fill]
    # check single course, but wait for elimilate later -----------    
    Threads.@threads for n_poll in eachindex(pollTypeS)
        fit_singleCoarse[:,n_poll,:]=medfilter(view(staObsOut_filled,:,n_poll,:);edgeLen=1)
    end
    misfit_singleCoarse = staObsOut_filled - fit_singleCoarse
    diff_neighbour[2:end-1,:,:]=abs.(staObsOut_filled[3:end,:,:]-staObsOut_filled[1:end-2,:,:]);    
    stdMultipler=1.2;
   
    diffData=abs.(diff(staObsOut_filled;dims=1));diffData=cat(dims=1,diffData[1:1,:,:],diffData)
    std_course=similar(diffData)
    Threads.@threads for n_poll in eachindex(pollTypeS)
        std_course[:,n_poll,:]=stdMultipler .* 1.4826*medfilter(view(diffData,:,n_poll,:);edgeLen=12*3)
    end        
    std_course=setStdLowerLimit(std_course,stdLowLimit_all)      
    std_course=lowerStd_forNegativeMisfit(std_course,misfit_singleCoarse;adjustRatio=0.8)
    pdf_singleCoarse=pdf.(Normal(),-abs.(misfit_singleCoarse./std_course));
    q_singleCoarse=.!isnan.(staObs) .& ( pdf_singleCoarse .< probabilityLimit ) .& (abs.(misfit_singleCoarse) .> diff_neighbour)

    # check temporal ===============
    # q_pm= in.(pollTypeS,[["PM25","PM10"]])
    q_fast = in.(pollTypeS,[["SO2","NO2"]])
    fit_time=similar(staObsOut)
    if any(q_fast)
        fit_time[:,q_fast,:]=nanFilter(staObsOut_filled[:,q_fast,:]; responsetype = Lowpass(1/5)) 
    end
    if any(.!q_fast)
        fit_time[:,.!q_fast,:]=nanFilter(staObsOut_filled[:,.!q_fast,:]; responsetype = Lowpass(1/8)) 
    end
    misfit_time= staObsOut - fit_time
    std_time=similar(misfit_time)
    if any(q_fast)
        std_time[:,q_fast,:]=1.3 .* movingMean(abs.(misfit_time[:,q_fast,:]);edgeLen=12)
    end
    if any(.!q_fast)
        std_time[:,.!q_fast,:]=1.3 .* movingMean(abs.(misfit_time[:,.!q_fast,:]);edgeLen=12*3)
    end
    # movingStd(abs.(fit_time - staObsOut);edgeLen=12*3) 
    std_time=setStdLowerLimit(std_time,stdLowLimit_all)  
    std_time=lowerStd_forNegativeMisfit(std_time,misfit_time;adjustRatio=0.5)
    pdf_time=pdf.(Normal(),-abs.(misfit_time./std_time));
    q_time= pdf_time .< probabilityLimit
        
    # check pollRatio ======================================= 
    fit_base=movingMean(abs.(staObsOut);edgeLen=24*6)
    ratioFit = staObsOut ./ fit_base
    ratioFit_ref = ratioFit[:,indexin(basePollName,pollTypeS),:]
    misfit_ratioFit = ratioFit - ratioFit_ref
    std_ratioMisfit=0.8 .* movingMean(abs.(misfit_ratioFit);edgeLen=24*6)
    # std_ratioMisfit=similar(misfit_ratioFit)
    # Threads.@threads for n_poll in eachindex(pollTypeS)
    #     std_ratioMisfit[:,n_poll,:]=1.4826*medfilter(abs.(view(misfit_ratioFit,:,n_poll,:));edgeLen=12*6)
    # end      
    # std_ratioMisfit[1:3,:,:] .= std_ratioMisfit[4:4,:,:] # deal with boundary
    # std_ratioMisfit[end-3:end,:,:] .= std_ratioMisfit[end-4:end-4,:,:]# deal with boundary
    # std_ratioMisfit=lowerStd_forNegativeMisfit(std_ratioMisfit,misfit_pollRatio;adjustRatio=0.6)
    pdf_pollRatio=pdf.(Normal(),-abs.(misfit_ratioFit./std_ratioMisfit));
    q_pollRatio= pdf_pollRatio .< probabilityLimit 

    # std_ratioMisfit=movingStd(abs.(fit_pollRatio - pollRatio);edgeLen=12*7) 
    # q_pollRatio,fit_pollRatio,std_ratioMisfit,pdf_pollRatio=mymapslices((x,fit,std)->
    #     checkTemporal(x;analysis=true,edgeLen_std=12*3,probabilityLimit=probabilityLimit,fit_curr=fit,std_curr=std,lowLimit_std=0.0)
    #     ,pollRatio,fit_pollRatio,std_ratioMisfit);

    # QC for stagnant check ==========================
    q_stag,misfit_stag,std_stag,pdf_stag=mymapslices((a,b,c)->
            checkStagnant(a,b,c;analysis=true,probabilityLimit=probabilityLimit),
            staObsOut,misfit_space,std_space
        );
    
    # -------------- set staObsOut to NaN for those do not pass the stagant, spatiotemporal, and calibaration checks -------------
    # staObsOut[q_stag .| (q_time .& q_space) .| (q_singleCoarse .& q_space) .| (pdf_space .< 1e-10 *  probabilityLimit)] .= NaN
    q_combine= q_stag .| # stagnant
                (q_time .& q_space) .| (q_singleCoarse .& q_space) .|   (q_pollRatio .& q_space  ) .|  ( q_pollRatio .& ( pdf_singleCoarse .< sqrt(probabilityLimit) ) .&  ( pdf_time .< sqrt(probabilityLimit) ) ) 
    # q_combine= q_stag .| (q_time .& q_space) .| (q_singleCoarse .& q_space) .|    (pdf_space .< 1e-10 *  probabilityLimit) .|  (q_pollRatio .& q_space)
    # special treatment for selected pollutants
    for (n_poll,poll) in enumerate(pollTypeS)
        if poll in ["SO2","NO2"]
            q_combine[:,n_poll,:].= (pdf_space[:,n_poll,:] .< 1e-10 *  probabilityLimit) .& 
                 ((pdf_time[:,n_poll,:] .< probabilityLimit)    .& (pdf_singleCoarse[:,n_poll,:] .< 1e-10 *  probabilityLimit) )
        end
        
        if poll in ["PM25","PM10"]
            q_rawPollRatio=pdf_pollRatio[:,n_poll,:] .< sqrt(probabilityLimit)
            q_combine[:,n_poll,:] = q_stag[:,n_poll,:] .| (q_rawPollRatio .& q_time[:,n_poll,:] .& q_space[:,n_poll,:]) .| 
                (( pdf_time[:,n_poll,:] .<  probabilityLimit)  .& ( pdf_space[:,n_poll,:] .< 1e-8 * probabilityLimit) ) .| # big spatial temporal error
                (q_rawPollRatio .& q_singleCoarse[:,n_poll,:] .& q_space[:,n_poll,:]) .|    (q_pollRatio[:,n_poll,:] .& q_space[:,n_poll,:]  ) .|  ( q_pollRatio[:,n_poll,:] .& ( pdf_singleCoarse[:,n_poll,:] .< sqrt(probabilityLimit) ) .&  ( pdf_time[:,n_poll,:] .< sqrt(probabilityLimit) ) ) 
            # q_combine[:,n_poll,:] .& (pdf_pollRatio[:,n_poll,:] .< 1e-2)
            # check pm25 outweight pm10
            ind_pm=indexin( ["PM25","PM10"],["PM25","PM10","SO2","NO2","CO","O3"])
            if all(.!isnothing(ind_pm)) # both pm25 and pm10 are avaliable
                q_outweight=staObsOut[:,ind_pm[1],:] .> staObsOut[:,ind_pm[2],:]
                q_combine[:,n_poll,:] = q_combine[:,n_poll,:] .| (q_outweight .& q_time[:,n_poll,:] .& q_space[:,n_poll,:]) 
            end
        end
        # if isequal(poll,"SO2")
        #     q_combine[:,n_poll,:].=false
        # end
    end
    staObsOut[q_combine] .= NaN
    
    # ================== Elimilate after ===========================
    # QC for isolation check ===================
    edgeLen=5
    num_noValid=getOutlierNum(staObsOut;edgeLen=edgeLen)
    # stdDecreaseRatio= max.(0.0,num_noValid .- edgeLen .+ 2) ./ (edgeLen .+ 3)
    stdDecreaseRatio= 0.5 .* max.(0.0,num_noValid ) ./ (2*edgeLen .+ 1)
    std_isolation=std_space .* (1 .- stdDecreaseRatio)
    pdf_isolate=pdf.(Normal(),-abs.(misfit_space./std_isolation));    
    q_isolate= ( num_noValid .> 0) .& (pdf_isolate .< probabilityLimit) .& (pdf_space .< 1e3 * probabilityLimit) .& ((pdf_pollRatio .<  1e3 * probabilityLimit ) .| (pdf_singleCoarse .< probabilityLimit) .| (pdf_time .< 1e2 * probabilityLimit))
    # extra dealling for so2
    q_so2= in.(pollTypeS,[["SO2"]])
    q_isolate[:,q_so2,:] .= false
    # # extra dealling for pm
    # q_pm= in.(pollTypeS,[["PM25","PM10"]])
    # q_isolate[:,q_pm,:] .= q_isolate[:,q_pm,:] .&  q_pollRatio[:,q_pm,:]
    staObsOut[q_isolate] .= NaN

    # ============================= output ================================
    if analysis 
        q_outlier= isnan.(staObsOut) .& .!isnan.(staObs)
        d_analysis =   Dict{String,Any}("q_outlier" => q_outlier,
            "q_outRange" => q_outRange,
            "q_singleCoarse" => q_singleCoarse,
            "q_time" => q_time,
            "q_space" => q_space,
            "q_stag" => q_stag,
            "q_isolate" => q_isolate,
            "misfit_singleCoarse" => misfit_singleCoarse,
            "fit_time" => fit_time,
            "fit_space" => fit_space,
            "ratioFit" => ratioFit,
            "misfit_ratioFit" => misfit_ratioFit,
            "misfit_stag" => misfit_stag,
            "std_time" => std_time,
            "std_space" => std_space,
            "std_stag" => std_stag,
            "std_ratioMisfit" => std_ratioMisfit,
            "pdf_singleCoarse" => pdf_singleCoarse,
            "pdf_time" => pdf_time,
            "pdf_space" => pdf_space,
            "pdf_stag" => pdf_stag,
            "pdf_pollRatio" => pdf_pollRatio,
            "pdf_isolate" => pdf_isolate,
            "fit_singleCoarse" => fit_singleCoarse
        )

        return staObsOut,d_analysis
    else
        return staObsOut
    end
end