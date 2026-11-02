using Distributions
# include("locUtilities.jl")

export checkSingleCoarse, checkDaysCoarse,checkDaysCoarse!, getMajority
"""
q_curr,misfit,std_curr,pdf_curr=checkSingleCoarse(inData;analysis=true);
"""
function checkSingleCoarse(inData::AbstractVector;edgeLen_std=12*3,lowLimit_std=1,probabilityLimit=1e-7,analysis=false,stdMutiplier=2.5,
    misfit=medfilter(inData;edgeLen=1) - inData ,
    std_curr=1.4826*medfilter(abs.(diff(inData));edgeLen=edgeLen_std) .* stdMutiplier
    )
    
    # get misfit and its std ---------------
    append!(std_curr,std_curr[end]) # fill the last std with the previous one
    replace!(x-> x<lowLimit_std ? lowLimit_std : x , std_curr)  # set lower limit for std
    # get diff_neighbour
    diff_neighbour=fill(NaN,size(inData))
    diff_neighbour[2:end-1]=abs.(inData[3:end]-inData[1:end-2]);
    pdf_curr=pdf.(Normal(),-abs.(misfit./std_curr));
    replace!(pdf_curr , NaN => 1.0)  # set NaN to 1

    # is coarse error when low probability and greater than the diff of obs nearby
    q_curr=( abs.(misfit) .> diff_neighbour) .&  (pdf_curr .< probabilityLimit);
    if analysis 
        return q_curr,misfit,std_curr,pdf_curr
    else
        return q_curr
    end
end

# function checkSingleCoarse(inData::AbstractArray{T};edgeLen_std=12*3,lowLimit_std=1,probabilityLimit=1e-7,analysis=false,stdMutiplier=1.5,
#     misfit=medfilter(inData;edgeLen=1) - inData ,
#     std_curr=1.4826*medfilter(abs.(diff(inData;dims=1));edgeLen=edgeLen_std) .* stdMutiplier
#     ) where T
function checkSingleCoarse(inData::AbstractArray{T};edgeLen_std=12*3,lowLimit_std=3,probabilityLimit=1e-7,analysis=false,stdMutiplier=1.8 ,
    misfit=medfilter(inData;edgeLen=1) - inData ,
    # std_curr=movingStd(abs.(diff(inData;dims=1));edgeLen=edgeLen_std) .* stdMutiplier
        std_curr=1.4826*medfilter(abs.(diff(inData;dims=1));edgeLen=edgeLen_std) .* stdMutiplier
    ) where T
    
    # get misfit and its std ---------------
    std_curr=vcat(std_curr,std_curr[[end],:]) # fill the last std with the previous one
    replace!(x-> x<lowLimit_std ? lowLimit_std : x , std_curr)  # set lower limit for std
    # get diff_neighbour
    diff_neighbour=fill(NaN,size(inData))
    diff_neighbour[2:end-1,:]=abs.(inData[3:end,:]-inData[1:end-2,:]);
    pdf_curr=pdf.(Normal(),-abs.(misfit./std_curr));
    replace!(pdf_curr , NaN => 1.0)  # set NaN to 1

    # is coarse error when low probability and greater than the diff of obs nearby
    q_curr=( abs.(misfit) .> diff_neighbour) .&  (pdf_curr .< probabilityLimit);
    if analysis 
        return q_curr,misfit,std_curr,pdf_curr
    else
        return q_curr
    end
end

# # test for checkSingleCoarse ------------------------------------------------------------
# using Revise,DataBaseM,MyPlotM
# plt=get_plt();plt.close("all")
# Times=DateTime([2016  3 15]...):Hour(1):DateTime([2016  5 20  ]...);
# polType="no2";staID2analysis="440100063";
# obsData=vec(getObsPoll(Times,polType,staID2analysis))
# q_coarse,misfit,std_curr,pdf_curr=checkSingleCoarse(obsData;analysis=true)
# plt.plot(Times,obsData,"-*")
# plt.plot(Times[q_coarse],obsData[q_coarse],"r*")
# # plt.xlim([DateTime(2016,4,15),DateTime(2016,5,1)])
# plt.gcf()

# quickPlot(obsData,misfit,std_curr,log.(pdf_curr))


"""
q_coarse,misfit_coarse,std_coarse,pdf_coarse=checkDaysCoarse(inData;analysis=true);
"""
function checkDaysCoarse(inData;edgeLen=24*5,lowLimit_std=10,probabilityLimit=1e-7,analysis=false,stdMutiplier=2,
    misfit=medfilter(inData;edgeLen=edgeLen)-inData,
    std_curr=1.4826*medfilter(abs.(misfit);edgeLen=edgeLen) .* stdMutiplier
    )
    replace!(x-> x < lowLimit_std ? lowLimit_std : x , std_curr) # set lower limit for std
    # std_curr[std_curr .< lowLimit_std] .= lowLimit_std
    # get diff_neighbour
    pdf_curr=pdf.(Normal(),-abs.(misfit./std_curr));
    pdf_curr[isnan.(pdf_curr)] .=1;

    # is coarse error when low probability and greater than the diff of obs nearby
    q_coarse = pdf_curr .< probabilityLimit
    if analysis
        return q_coarse,misfit,std_curr,pdf_curr
    else
        return q_coarse
    end
end

"""
No working, do not use this
"""
function checkDaysCoarse!(q_coarse,inData;edgeLen=24*5,lowLimit_std=10,probabilityLimit=1e-7,analysis=false,stdMutiplier=2,
    Tanalysis = ( misfit_curr=similar(inData) , std_curr=similar(inData) , pdf_curr=similar(inData) ),
    )
    misfit_curr,std_curr,pdf_curr = Tanalysis
    misfit_curr = medfilter!(misfit_curr,inData;edgeLen=edgeLen) - inData;
    std_curr = medfilter!(std_curr , abs.(misfit_curr) ; edgeLen=edgeLen) 
    @. std_curr = 1.4826 * std_curr .* stdMutiplier
    replace!(x-> x < lowLimit_std ? lowLimit_std : x , std_curr) # set lower limit for std
    std_curr[std_curr .< lowLimit_std] .= lowLimit_std
    # get diff_neighbour
    pdf_curr .= pdf.(Normal(),-abs.(misfit_curr./std_curr));
    pdf_curr[isnan.(pdf_curr)] .=1;
    # is coarse error when low probability and greater than the diff of obs nearby
    @. q_coarse = pdf_curr .< probabilityLimit
    
    if analysis 
        return q_coarse,Tanalysis
    else
        return q_coarse
    end
end
"""
function f(x)
    @. x = 1
    return x
end
x=rand(10)
@time f(x)

"""

function getMajority(inData,rangLimit;edgeLen=24*15,
    estimate_median=medfilter(x;edgeLen=edgeLen))
    # return true when data witin estimate +- rangeLimit

    q_major=@. abs.(  estimate_median - inData ) <= rangLimit;
    return q_major
end

