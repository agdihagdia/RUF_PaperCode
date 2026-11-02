using Distributions,NaNStatistics#,MyPlotM
# include("locUtilities.jl")
# include("checkCoarse.jl")
export checkCalibration,calibrationCheck, getCaliHour,topStepEstimate,enhanceDiurnal
# using the new checkCalibration instead of calibrationCheck
"""
check calibaration using checkSingleCoarse

q_cali = calibrationCheck(inData;probabilityLimit=1e-7,edgeLen = 3*24,step=24,stdMutiplier=1.5)
"""
function checkCalibration(inData;probabilityLimit=1e-7,stdMutiplier=1,analysis=false,edgeLen = 3*24,step=24)
    # inData should be a vector
    validIn=filter(!isnan,inData)
    if isempty(validIn) # return if all inputs are nan
        q_curr,misfit_enhanced , std_enhanced,pdf_enhanced= falses(size(inData)) , fill(NaN,size(inData)), fill(NaN,size(inData)), fill(NaN,size(inData))
    else
        diffIn=diff(inData)
        inRange = maximum(validIn) - minimum(validIn) 
        step_obs = minimum(filter(>(0),abs.(diffIn));init=Inf)
        if step_obs < inRange/64  # don't use stagnant check for station with large step /little unique observations
            q_ori,misfit_ori,std_ori,pdf_ori=checkSingleCoarse(inData;analysis=true,probabilityLimit=probabilityLimit,stdMutiplier=stdMutiplier);
            # q_ori=checkSingleCoarse(inData;probabilityLimit=probabilityLimit,stdMutiplier=stdMutiplier);
            inData_enhanced = enhanceDiurnal(inData;edgeLen=edgeLen,step=step)
            q_enhanced,misfit_enhanced,std_enhanced,pdf_enhanced=checkSingleCoarse(inData_enhanced;analysis=true);
            q_enhanced=checkSingleCoarse(inData_enhanced;probabilityLimit=probabilityLimit,stdMutiplier=stdMutiplier);
            q_curr= q_ori .& q_enhanced
        else # ignore the checkCalibration for observations with low resolutions
            q_curr,misfit_enhanced , std_enhanced,pdf_enhanced = falses(size(inData)) , fill(NaN,size(inData)), fill(NaN,size(inData)), fill(NaN,size(inData))
        end
    end
    if analysis 
        # return q_curr,misfit_ori,std_ori,pdf_ori,misfit_enhanced,std_enhanced,pdf_enhanced
        return q_curr , misfit_enhanced , std_enhanced,pdf_enhanced
    else
        return q_curr
    end
end
# quickPlot(inData,pdf_ori,pdf_enhanced)
function enhanceDiurnal(inData;edgeLen=3*24,step=24)
    outData=fill(NaN,size(inData))
    for n in eachindex(outData)
        # calulate the ind of start and end to mean
        if n > edgeLen
            s = n-edgeLen
        else
            s = n - step * fld(n-1,step)
        end
        if n+edgeLen <= length(inData)
            e = n+edgeLen
        else
            e = n + step * fld(length(inData) - n,step)
        end
        # println([n,edgeLen,s,step,e])
        # outData[n] = mean(filter(!isnan,inData[s:step:e]))
        outData[n] = nanmean(inData[s:step:e])
    end
    return outData
end
# # # test code for enhanceDiurnal
# inData=5 .*sin.(LinRange(0,2*pi,24*20))+rand(24*20)
# inData[8*24:24:15*24] .+= 5
# outData=enhanceDiurnal(inData;edgeLen=3*24,step=24)
# using MyPlotM
# plt=get_plt()
# plt.close("all")
# # plt.plot(inData,label="original")
# # plt.plot(outData,label="enhanced")
# plt.plot( hcat(inData,outData))
# # plt.legend()
# plt.gcf()

function calibrationCheck(inData;days2mean=14,probabilityLimit=1e-4)
    # 2017-03-02 add a method to check calibration by probability
    
    if rem(length(inData),24) !=0;error(["number of time in the observation should be times of 24 if you want to use calibration check"]);end
    q_cali=falses(size(inData));
    for row_startIni=1:Int(24*days2mean/2):size(inData,1)            
        row_end=row_startIni+24*days2mean -1;
        row_start=row_startIni;
        if row_end>size(inData,1)                
            row_start=row_start-(row_end-size(inData,1));
            if row_start<1;row_start=1;end
            row_end=size(inData,1;)
        end
        durnalMean=fill(NaN,24,size(inData,2));
        windowData=inData[row_start:row_end];
        if sum(isnan.(windowData))>length(windowData)/2;continue;end
        durnalMean=[nanpctile(x,50) for x in eachrow(reshape(windowData,24,:))]
        n_hour=getCaliHour(durnalMean;probabilityLimit=probabilityLimit);
        # println(typeof(n_hour))
        if !isnan(n_hour)
            q_cali[row_start-1+n_hour:24:row_end] .= true;
        end
    end
    return q_cali
end

function getCaliHour(inData;probabilityLimit=1e-4)
    # n_hour will be in 1 to 24, or nan
    # in should be a vector of 24 element
    # 2017-03-06 add support for nan in

    cycIn=[inData[end];inData;inData[1]];
    medFit_cycle=movingFrame(nanmedian,cycIn;edgeLen=1);
    medFit=medFit_cycle[2:end-1];
    misfit=abs.(inData - medFit);
    # println(nanmaximum(misfit),misfit)

    # candidate tor calibration hour
    n_hour=indexin(nanmaximum(misfit),misfit)[1]; 
    sorted_misfit=filter(!isnan,sort(misfit));

    # return normal if too many nan in inData
    if length(sorted_misfit)<20;n_hour=NaN;return n_hour;end# to many nan
    
    # println(sorted_misfit)
    estimate_maxMisfit=sorted_misfit[end-1];
    steps=nanminimum(abs.(diff(unique(inData))));
    estimate_maxMisfitStd=maximum([ mean(sorted_misfit[end-6:end-1]),#part of misfit
                                (nanmaximum(medFit)-nanminimum(medFit))/6,#sixth cycle
                                steps]);#minimum step
    if cdf( Normal(0,estimate_maxMisfitStd),-abs.(sorted_misfit[end]-estimate_maxMisfit) ) >  probabilityLimit
        n_hour=NaN;#mark normal
    end
    println(cdf( Normal(0,estimate_maxMisfitStd),-abs.(sorted_misfit[end]-estimate_maxMisfit) ))
    println(n_hour)

    return n_hour
end

function topStepEstimate(inData,n_hour)
    cycIn=[inData[end-1:end];in;inData[1:2]];
    cycStep=diff(cycIn);
    out=abs(cycStep(n_hour)-cycStep(n_hour+3))/2;
    return out
end

# # # test code for stagnantCheck
# inData=5 .*sin.(LinRange(0,2*pi,24*20))+rand(24*20)
# inData[3:24:end] .+= 1
# # n_hour=getCaliHour(durnalMean);
# q_cali=calibrationCheck(inData;probabilityLimit=1e-4);
# hcat(inData,q_cali)'
# plt=get_plt()
# plt.close("all")
# plt.plot(inData)
# xx=[1:length(inData);]
# plt.plot(xx[q_cali],inData[q_cali],"*r")

# plt.gcf()