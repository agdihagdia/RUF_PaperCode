
using Distributions,NaNStatistics,LocalizerM, Polyester
# include("locUtilities.jl")
export checkSpatial, indexOfAgreement, movingAgreement, movingAgreement!, estimateSpaceStd, estimateSpace, estimateSpace_old,estimateSpace_mean

"""
lonlat=Matrix(staInfo[:,["lon","lat"]])
localizer=localizerM(lonlat,lonlat)
q_curr,misfit,std_curr,pdf_curr=checkSpatial(inData,localizer;analysis=true);

cutoff_spatial = getVaringSpaceCutOff(localizer;validStaNum=10) /2
q_curr,misfit,std_curr,pdf_curr=checkSpatial(inData,localizer;cutoff_spatial=cutoff_spatial);
"""
function checkSpatial(inData::AbstractArray,localizer::LocalizerStr;edgeLen_estimate=24*15,edgeLen_std=24*3,lowLimit_std=1,probabilityLimit=1e-7,analysis=false,stdMutiplier=1,ind2ref_all=axes(inData,2),
    numRef = 10 # varing cutoff that ensure 10 sites to reference
    )
    # dimension of inData should be (Times,sites)

    # estimate by spatial continuity
    fit_curr = estimateSpace(inData,localizer;numRef=numRef,edgeLen=edgeLen_estimate,ind2ref_all=ind2ref_all)
    misfit = fit_curr - inData 
    std_curr=estimateSpaceStd(inData,misfit,localizer;edgeLen=edgeLen_std,numRef=numRef,ind2ref_all=ind2ref_all) .* stdMutiplier
    
    # # deal with std_curr
    # std_curr[std_curr .< lowLimit_std] .=lowLimit_std # set lower limit for std

    # calculate pdf_curr
    pdf_curr=pdf.(Normal(),-abs.(misfit./std_curr));
    replace!(pdf_curr , NaN => 1.0)  # set NaN to 1
    q_curr = pdf_curr .< probabilityLimit
    if analysis
        return q_curr,fit_curr,std_curr,pdf_curr
    else
        return q_curr
    end
end



function checkSpatial(inData::AbstractArray{T},stalonlat::AbstractArray;edgeLen_misfit=24*15,edgeLen_std_space=24*3,probabilityLimit=1e-7,analysis=false,disp=false,ind2ref_all=axes(inData,3),
    numRef = 10 # varing cutoff that ensure 10 sites to reference
    ) where T
    # dimension of inData should be (Time,varName,site)

    # prepare q2ref_all
    q2ref_all=falses(size(inData,3));q2ref_all[ind2ref_all] .= true;

    # check space
    if disp;    println("spatial...");    end
    q_space=falses(size(inData))
    fit_curr=fill(convert(T,NaN),size(inData))
    std_curr=fill(convert(T,NaN),size(inData))
    pdf_curr=fill(convert(T,NaN),size(inData))
    # Threads.@threads for n_var in axes(inData,2)
    for n_var in axes(inData,2)
        data_aVar=view(inData,:,n_var,:)
        q_validSta=vec(sum(.!isnan.(data_aVar),dims=1)) .> 0 # valid stations
        q_validRef=q_validSta .& q2ref_all
        localizer=localizerM(stalonlat[q_validSta,:],stalonlat[q_validRef,:])
        data_aVar_valid=view(inData,:,n_var,q_validSta)
        q_space[:,n_var,q_validSta],fit_curr[:,n_var,q_validSta],std_curr[:,n_var,q_validSta],pdf_curr[:,n_var,q_validSta] =
        checkSpatial(data_aVar_valid,localizer;edgeLen_estimate=edgeLen_misfit,edgeLen_std=edgeLen_std_space,analysis=true,ind2ref_all=findall(q2ref_all[q_validSta]),
            numRef = numRef ) 
    end
    misfit=fit_curr-inData

    # calculate pdf_curr
    pdf_curr=pdf.(Normal(),-abs.(misfit./std_curr));
    replace!(pdf_curr , NaN => 1.0)  # set NaN to 1
    q_curr = pdf_curr .< probabilityLimit

    # output
    if analysis 
        return q_curr,fit_curr,std_curr,pdf_curr
    else
        return q_curr
    end
end

function indexOfAgreement(x,y)
# x is simulation; y is observation
    yMean=nanmean(y);
    index=1-nansum(abs.(y-x))./nansum(abs.(x-yMean) + abs.(y-yMean));
    return index
end
# function movingAgreement(x,y;edgeLen=24*15)
#     x_loc,y_loc=copy.([x,y])
#     q_invalid= isnan.(x_loc) .| isnan.(y_loc)
#     x_loc[q_invalid] .= NaN
#     y_loc[q_invalid] .= NaN
#     yMean=movingMean(y;edgeLen=edgeLen)
#     index=1 .- movingMean(abs.(y-x);edgeLen=edgeLen) ./ movingMean(abs.(x-yMean) + abs.(y-yMean));
#     return index
# end

"""
x=rand(100);y=rand(100)
@code_warntype movingAgreement(x,y)


"""
function movingAgreement(x,y;edgeLen=24*15)
    yMean=movingMean(y;edgeLen=edgeLen)
    index=1 .- movingMean(abs.(y-x);edgeLen=edgeLen) ./ movingMean(abs.(x-yMean) + abs.(y-yMean);edgeLen=edgeLen);
    return index
end

function movingAgreement!(index,x,y;edgeLen=24*15,yMean=movingMean(y;edgeLen=edgeLen))
    T=eltype(index)
    
    # for length(inData) < edgeLen -1
    if length(index) < edgeLen + 1
        index = movingAgreement(x,y;edgeLen=edgeLen)
        return index
    end

    # curr_meanSum=nanmean(inData[1:edgeLen])
    x_y,x_yMean,y_yMean,sampleN=zero(T),zero(T),zero(T),zero(T),zero(T)
    # initialize ----------------------
    for e=1:edgeLen
        if !isnan(x[e]) && !isnan(y[e]) 
            sampleN +=1
            x_y += abs(x[e] - y[e])
            x_yMean += abs(x[e] - yMean[e])
            y_yMean += abs(y[e] - yMean[e])
            # out_curr += inData[i] .^2
        end
    end
    # moving -----------------------
    for n in eachindex(index)
        e=n+edgeLen
        sPrevious=n-edgeLen -1

        # add the data at the end point
        if e <= length(index)
            if !isnan(x[e]) && !isnan(y[e]) 
                sampleN += 1
                x_y += abs(x[e] - y[e])
                x_yMean += abs(x[e] - yMean[e])
                y_yMean += abs(y[e] - yMean[e])
            end
        end

        # remove the data before the start point
        if sPrevious > 0
            if !isnan(x[sPrevious]) && !isnan(y[sPrevious]) 
                sampleN -= 1
                if sampleN > 0
                    x_y -= abs(x[sPrevious] - y[sPrevious])
                    x_yMean -= abs(x[sPrevious] - yMean[sPrevious])
                    y_yMean -= abs(y[sPrevious] - yMean[sPrevious])
                else
                    x_y,x_yMean,y_yMean=zero(T),zero(T),zero(T)
                end
            end
        end

        # calculate the output
        if sampleN >1
            index[n] = 1 - x_y / (x_yMean + y_yMean)

        else
            index[n] = NaN
        end

        # println([n,sampleN,out_curr,outData[n]])
    end
    return index
end


function estimateSpaceStd(obs::AbstractArray{T},misfit,localizer;edgeLen=24*3,numRef=10,ind2ref_all=axes(inData,2)) where T
    # number of nearby station may not be enough for reliable std, use multiple time as well
    # estimate the observation by sites nearby
    # dimensions of inData should be in (time,sta)

    # input check ----------------------
    if size(obs,2) - 1 < numRef
        @warn  "number of reference stations should be larger than number of stations + 1, changing the numRef to $(size(obs,2) - 1)"
        numRef = size(obs,2) - 1
    end

    estimate_spaceStd=Array{T}(undef,size(misfit))
    movingStd_spaceMisfit=Array{T}(undef,size(misfit))
    misfit .= abs.(misfit)
    for n_sta in axes(misfit,2)
        movingStd_curr=view(movingStd_spaceMisfit,:,n_sta)
        @views movingStd_curr = movingStd!(movingStd_curr,misfit[:,n_sta];edgeLen=edgeLen) 
    end
    # movingStd_spaceMisfit=mymapslices_par(x-> movingStd(abs.(x);edgeLen=edgeLen),misfit)
    # for n_sta in axes(misfit,2)
    ind_inValidRef=Vector{Int}(undef,numRef)
    ind_nearby=Vector{Int}(undef,numRef)
    refRatioMatrix=Array{T}(undef,size(obs,1),numRef)
    yMean=Array{T}(undef,size(obs,1))
    @views for n_sta in axes(misfit,2)
    # @views Threads.@threads  for n_sta in axes(misfit,2)
        ind_inValidRef .= getNearest(localizer,n_sta,2:numRef+1); # exclude self
        ind_nearby .= ind2ref_all[ind_inValidRef]
        yMean=movingMean!(yMean,view(obs,:,n_sta);edgeLen=edgeLen)
        for (n_ind,ind) in enumerate(ind_nearby)
            refRatio_curr=view(refRatioMatrix,:,n_ind )
            @views refRatio_curr = movingAgreement!(refRatio_curr,obs[:,ind],obs[:,n_sta];edgeLen=edgeLen,yMean=yMean) 
        end
        refRatioMatrix = refRatioMatrix ./ nansum(refRatioMatrix,dims=2) # normalize
        estimate_spaceStd[:,n_sta] = nansum(movingStd_spaceMisfit[:,ind_nearby] .* refRatioMatrix,dims=2)
    end
    # to deal with fast change of concentration during trasport( heavy polution,sand storm)
    # estimate std by variation of spatial fit during a small time spand
    fit_spatial= obs + misfit
    function stepHight(x;edgeLen=3)
        return abs(nansum(x[end+1-edgeLen:end]) .- nansum(x[1:edgeLen]))/edgeLen
    end
    # movingFrame_basic(fun,inData;edgeLen=2)
    std_range=mapslices(x->movingFrame_basic(snip->stepHight(snip;edgeLen=3),x;edgeLen=3),fit_spatial;dims=1) ./ 10
    estimate_spaceStd .= max.(estimate_spaceStd,std_range)

    return estimate_spaceStd
end

# function estimateSpaceStd(obs::AbstractArray{T},misfit,localizer;edgeLen=24*3,numRef=10,ind2ref_all=axes(inData,2)) where T
# # number of nearby station may not be enough for reliable std, use multiple time as well
#     # estimate the observation by sites nearby
#     # dimensions of inData should be in (time,sta)

#     # input check ----------------------
#     if size(obs,2) - 1 < numRef
#         @warn  "number of reference stations should be larger than number of stations + 1, changing the numRef to $(size(obs,2) - 1)"
#         numRef = size(obs,2) - 1
#     end

#     estimate_spaceStd=Array{T}(undef,size(misfit))
#     movingStd_spaceMisfit=Array{T}(undef,size(misfit))
#     misfit .= abs.(misfit)
#     for n_sta in axes(misfit,2)
#         movingStd_curr=view(movingStd_spaceMisfit,:,n_sta)
#         @views movingStd_curr = movingStd!(movingStd_curr,misfit[:,n_sta];edgeLen=edgeLen) 
#     end
#     # movingStd_spaceMisfit=mymapslices_par(x-> movingStd(abs.(x);edgeLen=edgeLen),misfit)
#     # for n_sta in axes(misfit,2)
#     ind_inValidRef=Vector{Int}(undef,numRef)
#     ind_nearby=Vector{Int}(undef,numRef)
#     refRatioMatrix=Array{T}(undef,size(obs,1),numRef)
#     yMean=Array{T}(undef,size(obs,1))
#     @views for n_sta in axes(misfit,2)
#     # @views Threads.@threads  for n_sta in axes(misfit,2)
#         ind_inValidRef .= getNearest(localizer,n_sta,2:numRef+1); # exclude self
#         ind_nearby .= ind2ref_all[ind_inValidRef]
#         yMean=movingMean!(yMean,view(obs,:,n_sta);edgeLen=edgeLen)
#         for (n_ind,ind) in enumerate(ind_nearby)
#             refRatio_curr=view(refRatioMatrix,:,n_ind )
#             @views refRatio_curr = movingAgreement!(refRatio_curr,obs[:,ind],obs[:,n_sta];edgeLen=edgeLen,yMean=yMean) 
#         end
#         refRatioMatrix = refRatioMatrix ./ nansum(refRatioMatrix,dims=2) # normalize
#         estimate_spaceStd[:,n_sta] = nansum(movingStd_spaceMisfit[:,ind_nearby] .* refRatioMatrix,dims=2)
#     end
#     return estimate_spaceStd
# end
        
"""
lonlat=Matrix(staInfo[:,["lon","lat"]])
localizer=localizerM(lonlat,lonlat)
estimate_space= estimateSpace(inData,localizer;cutoff_spatial=20e3)
"""
function estimateSpace(inData::AbstractArray{T},localizer;edgeLen=24*15,numRef=10,ind2cal_all=axes(inData,2),ind2ref_all=axes(inData,2)) where T
    # estimate the observation by sites nearby

    # input check ----------------------
    if size(inData,2) - 1 < numRef
        @warn  "number of reference stations should be larger than number of stations + 1, changing the numRef to $(size(inData,2) - 1)"
        numRef = size(inData,2) - 1
    end

    # dimensions of inData should be in (time,sta)
    # inData_moveMean=mapslices(x-> movingMean(x;edgeLen=edgeLen),inData,dims=1)    
    inData_moveMean=Array{T}(undef,size(inData))
    for n_sta in axes(inData_moveMean,2)
        movingMean_curr=view(inData_moveMean,:,n_sta)
        movingMean_curr = movingMean!(movingMean_curr,view(inData,:,n_sta);edgeLen=edgeLen) 
    end
    # inData_moveMean=mymapslices_par(x-> movingMean(x;edgeLen=edgeLen),inData)

    estimate_space=Array{T}(undef,size(inData))
    # 每个站点的邻近站索引、参考权重和时间均值都属于本次迭代的
    # 可变工作区。不能在Threads.@threads外共享，否则多个线程会
    # 同时覆写这些数组，导致空间估计结果不可复现。
    Threads.@threads for n_sta in ind2cal_all
        # 以下对象均为当前迭代私有；不同线程只共享只读输入，
        # 并且只写estimate_space中互不重叠的站点列。
        ind_inValidRef = getNearest(localizer,n_sta,2:numRef+1) # exclude self
        ind_nearby = ind2ref_all[ind_inValidRef]
        refRatioMatrix = Array{T}(undef,size(inData,1),numRef)
        yMean = view(inData_moveMean,:,n_sta)

        for (n_ind,ind) in enumerate(ind_nearby)
            refRatio_curr=view(refRatioMatrix,:,n_ind )
            @views movingAgreement!(refRatio_curr,inData[:,ind],inData[:,n_sta];edgeLen=edgeLen,yMean=yMean)
        end

        refRatioMatrix .= refRatioMatrix ./ nansum(refRatioMatrix,dims=2) # normalize
        @views estimate_space[:,n_sta] = inData_moveMean[:,n_sta] .+ nansum((inData[:,ind_nearby] .- inData_moveMean[:,ind_nearby])  .* refRatioMatrix,dims=2)
    end
    return estimate_space
end

"""
lonlat=Matrix(staInfo[:,["lon","lat"]])
localizer=localizerM(lonlat,lonlat)
estimate_space= estimateSpace(inData,localizer;cutoff_spatial=20e3)
"""
function estimateSpace_mean(inData::AbstractArray{T},localizer;edgeLen=24*15,numRef=10,ind2cal_all=axes(inData,2),ind2ref_all=axes(inData,2)) where T
    # estimate the observation by sites nearby

    # input check ----------------------
    if size(inData,2) - 1 < numRef
        @warn  "number of reference stations should be larger than number of stations + 1, changing the numRef to $(size(inData,2) - 1)"
        numRef = size(inData,2) - 1
    end

    # dimensions of inData should be in (time,sta)
    # inData_moveMean=mapslices(x-> movingMean(x;edgeLen=edgeLen),inData,dims=1)    
    inData_moveMean=Array{T}(undef,size(inData))
    for n_sta in axes(inData_moveMean,2)
        movingMean_curr=view(inData_moveMean,:,n_sta)
        movingMean_curr = movingMean!(movingMean_curr,view(inData,:,n_sta);edgeLen=edgeLen) 
    end
    # inData_moveMean=mymapslices_par(x-> movingMean(x;edgeLen=edgeLen),inData)

    estimate_space=Array{T}(undef,size(inData))
    ind_inValidRef=Vector{Int}(undef,numRef)
    ind_nearby=Vector{Int}(undef,numRef)
    refRatioMatrix=Array{T}(undef,size(inData,1),numRef)
    yMean=Array{T}(undef,size(inData,1))
    # for n_sta in ind2cal_all
    # Threads.@threads for n_sta in ind2cal_all
    for n_sta in ind2cal_all
        ind_inValidRef .= getNearest(localizer,n_sta,2:numRef+1); # exclude self
        ind_nearby .= ind2ref_all[ind_inValidRef]
        # yMean=movingMean!(yMean,view(inData,:,n_sta);edgeLen=edgeLen)
        yMean .= inData_moveMean[:,n_sta]
        for (n_ind,ind) in enumerate(ind_nearby)
            refRatio_curr=view(refRatioMatrix,:,n_ind )
            @views refRatio_curr = movingAgreement!(refRatio_curr,inData[:,ind],inData[:,n_sta];edgeLen=edgeLen,yMean=yMean) 
        end
        refRatioMatrix .= refRatioMatrix ./ nansum(refRatioMatrix,dims=2) # normalize
        @views estimate_space[:,n_sta] =  nansum( inData[:,ind_nearby]  .* refRatioMatrix,dims=2)
    end
    return estimate_space
end

function estimateSpace_old(inData,localizer;edgeLen=24*15,cutoff_spatial=20e3,ind2cal_all=axes(inData,2),excludeSelf=false)
    # estimate the observation by sites nearby
    # dimensions of inData should be in (time,sta)
    estimate_space=fill(NaN,size(inData))
    for n_sta in ind2cal_all
        if cutoff_spatial isa Vector
            cutoff_curr=cutoff_spatial[n_sta]
        else
            cutoff_curr=cutoff_spatial
        end
        q_loc,loc_coef=getLocCoef(localizer,n_sta,cutoff_curr;nargout=2);
        if excludeSelf; replace!(loc_coef,1 => 0);end # exclue self
        ind_nearby=findall(q_loc)
        refRatioMatrix=fill(NaN,size(inData,1),length(ind_nearby))
        for (n_ind,ind) in enumerate(ind_nearby)
            refRatioMatrix[:,n_ind] = movingAgreement(inData[:,ind],inData[:,n_sta];edgeLen=edgeLen) 
        end        
        refRatioMatrix = refRatioMatrix .* vec(loc_coef)' 
        refRatioMatrix = refRatioMatrix ./ nansum(refRatioMatrix,dims=2) # normalize
        estimate_space[:,n_sta] =  nansum(inData[:,q_loc] .* refRatioMatrix,dims=2)
    end
    return estimate_space
end
