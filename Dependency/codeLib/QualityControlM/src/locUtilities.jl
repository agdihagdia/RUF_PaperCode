
using NaNStatistics,Polyester,LinearAlgebra,StaticArrays,DSP#,FastRunningMedian#,PyCall#,MyPlotM
export mymapslices, mymapslices! , mymapslices_par,  medfilter, medfilter! , setLowerLimit!,misfitStd, getObsRange, extendTrueAround, nanFilter
export movingFrame, windowValidNumber, movingMean, movingMean!, movingStd,movingStd!, movingVar,movingFrame_basic
using RunningMedianM

function mymapslices(fun,inData::Array{T,3}) where T
    # return fun(data[:,1,1])
    # q_curr=BitArray(undef,size(data))
    # misfit_curr=Array{Float64,3}(undef,size(data))
    # std_curr=Array{Float64,3}(undef,size(data))
    # pdf_curr=Array{Float64,3}(undef,size(data))
    q_curr=falses(size(inData))
    misfit_curr=similar(inData)
    std_curr=similar(inData)
    pdf_curr=similar(inData)
    Threads.@threads for n_var in axes(inData,2)
        for n_sta in axes(inData,3)
    # Threads.@threads for n_sta in axes(inData,3)
    #     for n_var in axes(inData,2)
            # println(n_var,"  ",n_sta)
            q_curr[:,n_var,n_sta],misfit_curr[:,n_var,n_sta],std_curr[:,n_var,n_sta],pdf_curr[:,n_var,n_sta] =
                fun(@view inData[:,n_var,n_sta])
        end
    end
    return q_curr,misfit_curr,std_curr,pdf_curr
end



function mymapslices!(fun!,inData::Array{T,3}) where T
    q_curr=falses(size(inData))
    misfit_curr=similar(inData)
    std_curr=similar(inData)
    pdf_curr=similar(inData)
    Threads.@threads for n_sta in axes(inData,3)
    # for n_sta in axes(inData,3)
        for n_var in axes(inData,2)
            # println(n_var,"  ",n_sta)
            Tana_curr=( misfit_curr=view(misfit_curr,:,n_var,n_sta) , std_curr=view(std_curr,:,n_var,n_sta) , pdf_curr=view(pdf_curr,:,n_var,n_sta) )
            q_curr[:,n_var,n_sta],Tana_curr = fun!( view(q_curr,:,n_var,n_sta) , view(inData,:,n_var,n_sta) , Tana_curr )
        end
    end
    Tanalysis = ( misfit_curr=misfit_curr , std_curr=std_curr , pdf_curr=pdf_curr )
    return q_curr , Tanalysis
end


function mymapslices_par(fun,inData)
    outData=similar(inData)
    @batch per=thread  for n_sta in axes(inData,3)
    # for n_sta in axes(inData,3)
        for n_var in axes(inData,2)
            # println(n_var,"  ",n_sta)
            outData[:,n_var,n_sta] =  fun(@view inData[:,n_var,n_sta])
        end
    end
    return outData
end


function mymapslices(fun,inData,B,C)
    q_curr=falses(size(inData))
    misfit_curr=similar(inData)
    std_curr=similar(inData)
    pdf_curr=similar(inData)
    Threads.@threads for n_sta in axes(inData,3)
        for n_var in axes(inData,2)
            q_curr[:,n_var,n_sta],misfit_curr[:,n_var,n_sta],std_curr[:,n_var,n_sta],pdf_curr[:,n_var,n_sta] =
                fun(view(inData,:,n_var,n_sta),view(B,:,n_var,n_sta),view(C,:,n_var,n_sta))
            # fun(inData[:,n_var,n_sta],B[:,n_var,n_sta],C[:,n_var,n_sta])
        end
    end
    return q_curr,misfit_curr,std_curr,pdf_curr
end



"""
Fast version to calculate the pdf of 2-dimensional multivariate normal distribution

Test codes:
μ, Σ = [0, 0], [1 .4 ;.4 1]
xAll=[randn(2) for i=1:1000]
@time a=pdf_2dMvNormal(μ, Σ,xAll)
using Distributions
@time b=pdf.([MvNormal(μ, Σ )],eachcol(xAll))


"""
function pdf_2dMvNormal(μ, Σ,xAll::AbstractArray{T}) where T
    # dimension of the first dimension in xAll should be 2
    # dimension of μ should be (2,)
    # dimension of Σ should be (2,2)
    # return the pdf of x, which is a vector of length n_sample
    sizeIn=size(xAll)
    xAll=reshape(xAll,2,:)
    tmp=Array{T}(undef,size(xAll)[2:end])
    invΣ=@MMatrix zeros(2,2)
    invΣ .= Σ^-1
    x= @MVector zeros(2)
    @inbounds for n in eachindex(tmp)
        x .= view(xAll,:,n) .- μ
        tmp[n]=-0.5 * (x[1]^2 * invΣ[1,1] + x[2]^2 * invΣ[2,2] + x[1]*x[2]*(invΣ[1,2]+invΣ[2,1]) )
    end
    if length(sizeIn) >1
        tmp=reshape(tmp,sizeIn[2:end])
    else
        tmp=tmp[1]
    end
    return (2*pi)^-1 .* det(Σ)^-0.5 .* exp.(tmp)
end
"""
support for both vector and array a, filter alone the first dimension

"""
function nanFilter(a;
    responsetype = Lowpass(1/8),
    # responsetype = Lowpass(1/(60*60*8); fs=1/60/60),
    designmethod = Butterworth(24)
    ) 
    sizeIn=size(a)
    a=reshape(a,size(a,1),:)
    # filter that ignore nan, support for mutidimentional input (a). Filter along the first dimension
    loc_filter=digitalfilter(responsetype, designmethod)    

    a=copy(a)
    # Flat function with NaNs for comparison.
    on = ones(size(a));
    # Find all the NaNs in the input.
    n = isnan.(a);
    # Replace NaNs with zero, both in 'a' and 'on'.
    a[n] .= 0
    on[n] .= 0;

    c=similar(a)
    local flat
    for i in axes(c,2)
        # If all the values are NaN, skip this column.
        if sum(n[:,i]) == length(n[:,i])
            c[:,i] .= NaN
            continue            
        end
        # Calculate what a 'flat' function looks like after convolution.
        flat = filtfilt(loc_filter,on[:,i])
        # Do the actual convolution
        c[:,i] = filtfilt(loc_filter,a[:,i])./flat;
    end
    c=reshape(c,sizeIn)
    return c
end


function medfilter(inData;edgeLen=24*2)
    # running_median that ignore nan
    q_nan=isnan.(inData)
    # x=copy(inData)
    inData[q_nan] .= -Inf
    a=running_median(inData, edgeLen * 2 + 1)
    inData[q_nan] .= Inf
    b=running_median(inData, edgeLen * 2 + 1)
    inData[q_nan] .= NaN
    # window without valid data will return NaN, because Inf - Inf = NaN
    return (a + b)./2
end

function medfilter!(outData,inData;edgeLen=24*2)
    # running_median that ignore nan
    q_nan=isnan.(inData)
    # x=copy(inData)
    inData[q_nan] .= -Inf
    a=running_median(inData, edgeLen * 2 + 1)
    inData[q_nan] .= Inf
    b=running_median(inData, edgeLen * 2 + 1)
    inData[q_nan] .= NaN
    
    a=running_median(replace(nData,NaN=>-Inf), edgeLen * 2 + 1)
    b=running_median(replace(nData,NaN=>Inf), edgeLen * 2 + 1)
    # window without valid data will return NaN, because Inf - Inf = NaN
    @. outData = (a +b)/2
    outData .= ( running_median(replace(nData,NaN=>-Inf), edgeLen * 2 + 1) .+ running_median(replace(nData,NaN=>Inf), edgeLen * 2 + 1) )./2
    return outData
    # return 
end

# function medfilter(inData;edgeLen=24*2)
#     sizeIn=size(inData)
#     inMatrix=reshape(inData,sizeIn[1],:)

#     # rolling accros the first dimension, two dimensional matrix is supported
#     pd=pyimport("pandas")
#     df=pd.DataFrame(inMatrix)
#     out=df.rolling(edgeLen*2+1,center=true,min_periods=1).median().to_numpy()
   
#     return reshape(out,sizeIn)
# end

# test code 
# a=5 .*sin.(LinRange(0,20*pi,1000))+rand(1000)
# a[14:27] .= -Inf
# a[:4] .= -Inf
# # b= medfilter(a;edgeLen=5)
# b= running_median(a, 5 * 2 + 1)
# using PyPlot
# plt=get_plt()
# plt.clf()
# plt.plot(a,labels="input")
# plt.plot(b,labels="output")

# plt.xlim([0,50])
# plt.gcf()

function setLowerLimit!(inData;lowerLimit=fill(0.3,size(inData[1,:,:])))
    # set value lower than limit or NaN to lowerLimit
    for n_var in axes(lowerLimit,2)
        for n_sta in axes(lowerLimit,3)
            q_lower = .! inData[:,n_var,n_sta] .>= lowerLimit[1,n_var,n_sta]
            inData[q_lower,n_var,n_sta ] .= lowerLimit[1,n_var,n_sta]
        end
    end
end

# function movingFrame(fun,inData;edgeLen=2)
#     outData=fill(NaN,size(inData))
#     for n in eachindex(inData)
#         s,e=max(1,n-edgeLen),min(length(inData),n+edgeLen)
#         @inbounds outData[n]=fun(view(inData,s:e))
#     end
#     return outData
# end
    
# # test code for stagnantCheck
# using NaNStatistics,MyPlotM
# plt=get_plt()
# plt.close("all")
# inData=5 .*sin.(LinRange(0,2*pi,24*20))+rand(24*20)
# inData[6:12] .= NaN
# inData[15:16] .= NaN
# @time q_lowDiffDiff=movingFrame(nanmedian,inData);
# plt.plot(hcat(inData,q_lowDiffDiff))
# plt.xlim([0,30])
# plt.gcf()



function movingMean(misfit;edgeLen=50)
    fun_newElement(x)=x
    sumVar,sampleNum=movingFrame(misfit,fun_newElement;edgeLen=edgeLen)
    # sampleNum[sampleNum .== 0] .= NaN # assigning NaN to elements with not valid sample
    replace!(sampleNum,0 => NaN) # assigning NaN to elements with not valid sample
    return sumVar./sampleNum
end


function movingValid!(outData,inData;edgeLen=50) 
    fun_newElement(x)= NaN
    fun_window(x,n)=n
    outData=movingFrame!(outData,inData,fun_newElement,fun_window;edgeLen=edgeLen)
    return outData
end

function movingMean!(outData,inData;edgeLen=50) 
    fun_newElement(x)=x
    fun_window(x,n)=x/n
    outData=movingFrame!(outData,inData,fun_newElement,fun_window;edgeLen=edgeLen)
    return outData
end
"""
a=rand(1000)
b=similar(a)
@time movingMean!(b,a)
@time c=movingMean(a);
c == b
"""

function movingStd!(outData,inData;edgeLen=50) 
    fun_newElement(x)=x^2
    # fun_window(x,n)=sqrt(x/n)# x could be negative, debug latter
    fun_window(x,n)=sqrt(max(0.0,x)/n)
    outData=movingFrame!(outData,inData,fun_newElement,fun_window;edgeLen=edgeLen)
    return outData
end
# function movingFrame!(outData::Vector{T},inData::Vector{T},fun_newElement;edgeLen=50,sampleMean=true) where T

#     # for length(inData) < edgeLen -1
#     if length(inData) < edgeLen + 1
#         if sampleMean
#             sampleN=length(filter(!isnan,inData))
#         else
#             sampleN=one(T)
#         end
#         outData .= sum(fun_newElement.(filter(!isnan,inData))) ./sampleN
#         return outData
#     end

#     # curr_meanSum=nanmean(inData[1:edgeLen])
#     out_curr=zero(T)
#     sampleN=zero(T)
#     # initialize ----------------------
#     for e=1:edgeLen
#         if !isnan(inData[e])
#             sampleN +=1
#             out_curr += fun_newElement(inData[e])
#             # out_curr += inData[i] .^2
#         end
#     end
#     # moving -----------------------
#     for n in eachindex(inData)
#         e=n+edgeLen
#         sPrevious=n-edgeLen -1

#         # add the data at the end point
#         if e <= length(inData)
#             if !isnan(inData[e])
#                 # curr_meanSum += inData[e]
#                 out_curr += fun_newElement(inData[e])
#                 # out_curr += inData[e] .^2
#                 sampleN += 1
#             end
#         end

#         # remove the data before the start point
#         if sPrevious > 0
#             if !isnan(inData[sPrevious])
#                 # curr_meanSum -= inData[sPrevious]
#                 out_curr -= fun_newElement(inData[sPrevious])
#                 # out_curr -= inData[sPrevious].^2
#                 sampleN -= 1
#             end
#         end

#         # calculate the output
#         if sampleN >1
#             if sampleMean
#                 outData[n] = out_curr /sampleN
#             else
#                 outData[n] = out_curr
#             end
#         else
#             outData[n] = NaN
#         end

#         # println([n,sampleN,out_curr,outData[n]])
#     end
#     return outData
# end

function movingFrame(inData,fun_newElement;edgeLen=50)
    T=eltype(inData)
    if length(inData) < edgeLen + 1
        sampleN=length(filter(!isnan,inData))
        out_sampleN = fill(convert(eltype(inData),sampleN) ,size(inData))
        outData = fill(sum(fun_newElement.(filter(!isnan,inData))) ,size(inData))
        return outData,out_sampleN
    end

    # curr_meanSum=nanmean(inData[1:edgeLen])
    out_curr=zero(T)
    sampleN=zero(T)
    outData=Array{T}(undef,size(inData))
    out_sampleN=Array{T}(undef,size(inData)) # intitialize as float64 to avoid future type convertion

    # initialize ----------------------
    for e=1:edgeLen
        if !isnan(inData[e])
            sampleN +=1
            out_curr += fun_newElement(inData[e])
        end
    end
    # moving -----------------------
    for n in eachindex(inData)
        e=n+edgeLen
        sPrevious=n-edgeLen -1

        # add the data at the end point
        if e <= length(inData)
            if !isnan(inData[e])
                out_curr += fun_newElement(inData[e])
                sampleN += 1
            end
        end

        # remove the data before the start point
        if sPrevious > 0
            if !isnan(inData[sPrevious])
                out_curr -= fun_newElement(inData[sPrevious])
                sampleN -= 1
            end
        end

        # calculate the output
        out_sampleN[n] = sampleN
        if sampleN >1
            outData[n] = out_curr
        else
            outData[n] = NaN
        end

        # println(n,"  ",sampleN,"   ",out_curr)
    end
    return outData,out_sampleN
end
"""
fun_newElement to add the new element to value_cur
fun_window to calculate the output from value_cur and sampleN
sample usage for movingMean
    fun_newElement(x)=x
    fun_window(x,n)=x/n
    outData=movingFrame!(outData,inData,fun_newElement,fun_window;edgeLen=edgeLen)

"""
function movingFrame!(outData,inData,fun_newElement,fun_window;edgeLen=50)
    T=eltype(inData)
    if length(inData) < edgeLen + 1
        sampleN=length(filter(!isnan,inData))
        outData .= fun_window(sum(fun_newElement.(filter(!isnan,inData))),sampleN)
        return outData
    end

    # curr_meanSum=nanmean(inData[1:edgeLen])
    value_cur=zero(Float64)
    sampleN=zero(T)

    # initialize ----------------------
    for e=1:edgeLen
        if !isnan(inData[e])
            sampleN +=1
            value_cur += fun_newElement(inData[e])
        end
    end
    # moving -----------------------
    for n in eachindex(inData)
        e=n+edgeLen
        sPrevious=n-edgeLen -1

        # add the data at the end point
        if e <= length(inData)
            if !isnan(inData[e])
                value_cur += fun_newElement(inData[e])
                sampleN += 1
            end
        end

        # remove the data before the start point
        if sPrevious > 0
            if !isnan(inData[sPrevious])
                value_cur -= fun_newElement(inData[sPrevious])
                sampleN -= 1
            end
        end

        # calculate the output
        if sampleN >=1
            outData[n] = fun_window(value_cur,sampleN)
        else
            # value_cur=zero(T)# reset the value_cur
            outData[n] = NaN
        end
        # println(n,"  ",sampleN,"   ",out_curr)
    end
    return outData
end

function windowValidNumber(inData;edgeLen=12)
    sumVar,sampleNumber=movingFrame(inData,x->NaN;edgeLen=edgeLen)
    return sampleNumber
end

function movingStd(misfit;edgeLen=50)
    fun_newElement(x)=x^2
    sumVar,sampleNum=movingFrame(misfit,fun_newElement;edgeLen=edgeLen)
    sumVar[sumVar .< 0 ] .= 0 # set lower limit, negative value might occur due to resolution error
    # sampleNum[sampleNum .==0] .= NaN # assigning NaN to elements with not valid sample
    replace!(sampleNum,0 => NaN) # assigning NaN to elements with not valid sample
    return sqrt.(sumVar./sampleNum)
end
function movingFrame_basic(fun,inData;edgeLen=2)
    outData=similar(inData)
    for n in eachindex(inData)
        s,e=max(1,n-edgeLen),min(length(inData),n+edgeLen)
        outData[n]=fun(inData[s:e])
    end
    return outData
end
# misfit=vcat(randn(200),randn(200) .*5,randn(200).*2,randn(200),randn(200) .*5,randn(200).*2)
# using MyPlotM,NaNStatistics
# plt=get_plt()
# plt.close("all")
# @time stdall=movingStd(misfit;edgeLen=50);
# # @time stdall=movingVar(misfit;edgeLen=50);
# plt.plot(hcat(abs.(misfit),stdall))
# plt.gcf()


function movingVar(misfit::Vector;edgeLen=50)
    # curr_meanSum=nanmean(misfit[1:edgeLen])
    if edgeLen > length(misfit);edgeLen=length(misfit);end
    curr_varSum=nansum(misfit[1:edgeLen].^2)
    sampleN=sum(.! isnan.(misfit[1:edgeLen]))
    out_var=fill(NaN,size(misfit))
    for n in eachindex(misfit)
        e=n+edgeLen
        sPrevious=n-edgeLen -1

        # add the data at the end point
        if e <= length(misfit)
            if !isnan(misfit[e])
                # curr_meanSum += misfit[e]
                curr_varSum += misfit[e].^2
                sampleN += 1
            end
        end

        # remove the data before the start point
        if sPrevious > 0
            if !isnan(misfit[sPrevious])
                # curr_meanSum -= misfit[sPrevious]
                curr_varSum -= misfit[sPrevious].^2
                sampleN -= 1
            end
        end

        # calculate the variance
        if sampleN >1
            out_var[n] = curr_varSum/(sampleN -1)
        else
            out_var[n] = 0
        end
        # println(sampleN,"   ",curr_varSum)
    end
    return out_var
end

# misfit=vcat(randn(200),randn(200) .*5,randn(200).*2,randn(200),randn(200) .*5,randn(200).*2)
# using NaNStatistics,MyPlotM
# plt=get_plt()
# plt.close("all")
# @time varall=sqrt.(movingVar(misfit;edgeLen=50));
# plt.plot(hcat(abs.(misfit),varall))
# plt.gcf()


function misfitStd(misfit)
    # std without subtract by mean, deal with NaN
    num=sum(.!isnan.(misfit));
    if num<2;return NaN;end
    std=sqrt(nansum(misfit.^2)./(num-1));
    return std
end


function getObsRange(varName)
    # convert the variable name to lowercase
    varName = uppercase.(varName)
    # use a dictionary to store the ranges for each variable
    if varName in ["SO2","SO2D"] 
        ranges = [0.0, 1071.0]
    elseif varName in ["NO2","NO2D"]
        ranges = [0.0, 1026.0]
    elseif varName in ["CO","COD"]
        ranges = [0.0, 62.5]
    elseif varName in ["O3","O3H"]
        ranges = [0.0, 1428.0]
    elseif varName in ["PM25","PM25D"]
        ranges = [0.0, 1000.0]
    elseif varName in ["PM10","PM10D"]
        # ranges = [0.0, 1000.0]
        ranges = [0.0, 10000.0] # some instruments have a upper limit of 10000
    end

    return ranges
end
function extendTrueAround(inData,extendSteps)
    out=inData;
    for n_extend=1:extendSteps
        out[1:end-n_extend]=inData[1+n_extend:end]|out[1:end-n_extend];
        out[1+n_extend:end]=inData[1:end-n_extend]|out[1+n_extend:end];
    end
    return out
end
"""
calculate the pdf for 2d Multivariate_normal_distribution according to  https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    
sample usage:
μ, Σ = [0, 0], [1 .4 ;.4 1]
xAll=randn(2,1000)
a=pdf_2dMvNormal(μ, Σ,xAll)

"""