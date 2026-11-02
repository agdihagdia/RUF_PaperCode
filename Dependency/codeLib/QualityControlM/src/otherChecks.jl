export checkIsolation
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

function getOutlierNum(inData::AbstractArray{T};edgeLen=12) where T
    # check observation with little valid observations around
    outlierData= ones(size(inData));
    outlierData[.!isnan.(inData)] .= NaN
    outlierNum = zeros(size(outlierData))
    Threads.@threads   for n_var in axes(outlierData,2)
        for n_sta in axes(outlierData,3)
            outlierNum_series=view(outlierNum,:,n_var,n_sta)
            outlierNum_series=movingValid!(outlierNum_series,view(outlierData,:,n_var,n_sta);edgeLen=edgeLen)
        end
    end
    # outlierNum= edgeLen * 2 + 1 .- validNum
    # q_isolate =  .!isnan.(inData) .& ( windowValidNumber(inData;edgeLen=edgeLen) .< validNum) 
    return outlierNum
end
