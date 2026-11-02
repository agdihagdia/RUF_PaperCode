# include("locUtilities.jl")
export checkTemporal


"""
q_time,misfit_time,std_time,pdf_time=checkTemporal(inData;analysis=true);
"""
# function checkTemporal(inData;edgeLen_std=12*3,lowLimit_std=0.1,probabilityLimit=1e-7,analysis=false,stdMutiplier=1,
#     misfit=nanFilter(inData; responsetype = Lowpass(1/(60*60*8); fs=1/60/60)) - inData ,
#     std_curr=movingStd(abs.(misfit);edgeLen=edgeLen_std) .* stdMutiplier
#     )

#     # deal with std_curr
#     std_curr[std_curr .< lowLimit_std] .=lowLimit_std # set lower limit for std. zero std might occur due to stagnant error

#     # calculate pdf_curr
#     pdf_curr=pdf.(Normal(),-abs.(misfit./std_curr));
#     replace!(pdf_curr , NaN => 1.0)  # set NaN to 1

#     q_curr = pdf_curr .< probabilityLimit
#     if analysis 
#         return q_curr,misfit,std_curr,pdf_curr
#     else
#         return q_curr
#     end
# end
function checkTemporal(inData;edgeLen_std=12*3,lowLimit_std=0.1,probabilityLimit=1e-7,analysis=false,stdMutiplier=1,
    fit_curr=nanFilter(inData; responsetype = Lowpass(1/(60*60*8); fs=1/60/60))  ,
    std_curr=movingStd(abs.(fit_curr - inData);edgeLen=edgeLen_std) .* stdMutiplier
    )
    misfit_curr=fit_curr - inData

    # deal with std_curr
    std_curr[std_curr .< lowLimit_std] .=lowLimit_std # set lower limit for std. zero std might occur due to stagnant error

    # calculate pdf_curr
    pdf_curr=pdf.(Normal(),-abs.(misfit_curr./std_curr));
    replace!(pdf_curr , NaN => 1.0)  # set NaN to 1

    q_curr = pdf_curr .< probabilityLimit
    if analysis 
        return q_curr,fit_curr,std_curr,pdf_curr
    else
        return q_curr
    end
end