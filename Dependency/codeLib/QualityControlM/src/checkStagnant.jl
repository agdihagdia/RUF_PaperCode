using NaNStatistics#,MyPlotM
export diffCheck, diffdiffCheck, stagnantCheck, sequenceAbsMean,checkStagnant
function diffCheck(inData;sameNlimit=5,diffLimit=nanpctile(abs.(diff(inData)),20)/8)
    # inData  should be a vector
    n_start=NaN;nanspace=0;valueBefore=0
    q_lowDiff=falses(size(inData))
    for n_step=1:(length(inData)-1)
        # check start of low diff ------------------------
        if  isnan(n_start)
            diffCurr=abs(inData[n_step+1]-inData[n_step])
            if !isnan(diffCurr) &&   diffCurr  <=  diffLimit# same with next
                n_start=n_step; # mark abnormal start
                nanspace=0;
                valueBefore=inData[n_step]
            end
        else
            diffCurr=abs(inData[n_step+1]-valueBefore)            
            # check ongoning of low diff -----------------------
            if  isnan(diffCurr)
                # mark nan as low diff;
                nanspace=nanspace+1;
            else # not nan for diffCurr
                if diffCurr <= diffLimit  && (n_step < length(inData)-1)
                    # lowe diff continue
                    valueBefore=inData[n_step]
                else
                # check end of low diff -----------------------
                    if n_step >= length(inData)-1
                        n_end_candidate=n_step +1
                    else
                        n_end_candidate=n_step
                    end
                    n_end=NaN # mark end, skip nan
                    for n=n_end_candidate:-1:n_start
                        if !isnan(inData[n])
                            n_end=n
                            break
                        end
                    end
                    if n_end-n_start+1>=sameNlimit+nanspace
                        q_lowDiff[n_start:n_end] .= true;
                    end
                    n_start=NaN;# mark normal
                end
            end
        end
    end
    return q_lowDiff
end

# function diffCheck(inData;diffIn=diff(inData),sameNlimit=5,diffLimit=nanpctile(abs.(diffIn),20)/8)
#     # inData  should be a vector

#     q_lowDiff=falses(size(inData));
#     n_start=NaN;nanspace=0;
#     for n_step=1:(length(diffIn)-1)
#         if  abs.(diffIn[n_step])  <=  diffLimit# same with next
#             if isnan(n_start) # normal
#                 n_start=n_step; # mark abnormal start
#                 nanspace=0;
#             end
#         elseif  isnan(diffIn[n_step])
#             # mark nan as low diff
#             if !isnan(n_start) && isnan(diffIn[n_step]);nanspace=nanspace+1;continue;end
#         else# not the same with next
#             if !isnan(n_start) # end of a same period
#                 n_end=n_step; # mark abnormal end. add 1 because of diff 
#                 if n_end-n_start+1>=sameNlimit+nanspace
#                     q_lowDiff[n_start:n_end] .= true;
#                 end
#                 n_start=NaN;# mark normal
#             end
#         end
#     end
#     return q_lowDiff
# end

# # test code for diffCheck
# inData=5 .*sin.(LinRange(0,2*pi,100))+rand(100)
# inData[4:10] .= 0
# diffCheck(inData)'


function diffdiffCheck(inData;sameNlimit=9,diffIn=diff(inData),
    diffDiffLimit=nanpctile(abs.(diff(diffIn)),25)/2)
    temp_lowDiffDiff=diffCheck(diffIn;sameNlimit=sameNlimit,diffLimit=diffDiffLimit)
    q_lowDiffDiff=falses(size(inData));
    for n in eachindex(q_lowDiffDiff)
        if n==1
            q_lowDiffDiff[1] = temp_lowDiffDiff[1]
        elseif n == length(q_lowDiffDiff)
            q_lowDiffDiff[end] = temp_lowDiffDiff[end]
        else
            if temp_lowDiffDiff[n] || temp_lowDiffDiff[n-1]
                q_lowDiffDiff[n]=true
            else
                q_lowDiffDiff[n]=false
            end
        end
    end
    return q_lowDiffDiff
end


# function diffdiffCheck(inData;sameNlimit=9,diffIn=diff(inData),
#     diffDiffLimit=nanpctile(abs.(diff(diffIn)),25)/2)


# q_lowDiffDiff=falses(size(inData));
# n_start=NaN;nanspace=0;
# for n_step=1:(length(diffIn)-1)        
# if  abs.(diffIn[n_step]-diffIn[n_step+1])  <  diffDiffLimit# is low diffdiff
#     if isnan(n_start) # normal
#         n_start=n_step; # mark abnormal start
#         nanspace=0;
#     end
#     if nanspace>0 # check if we can concatenate two piece of stagnant
#         if abs.(inData[n_step]-inData[n_step-nanspace+1])  > nanspace*diffDiffLimit # can not concatenate
#             n_end=n_step-1; # mark abnormal end. add 2 because of diff and compare to next diffInf
#             if n_end-n_start+1>=sameNlimit+nanspace
#                 q_lowDiffDiff[n_start:n_end] .= true;
#             end
#             n_start=NaN;# mark normal
#         end
#         nanspace=0;
#     end
# else# not the same with next
#     if !isnan(n_start) && isnan(abs.(diffIn[n_step]-diffIn[n_step+1]));nanspace=nanspace+1;continue;end
#     if !isnan(n_start) # end of a same period
#         if nanspace>0
#             n_end=n_step-1;
#         else
#             n_end=n_step+1; # mark abnormal end. add 2 because of diff and compare to next diffInf
#         end
#         if sum(.!isnan.(inData[n_start:n_end]))>=sameNlimit
#             q_lowDiffDiff[n_start:n_end] .= true;
#         end
#         n_start=NaN;# mark normal
#         nanspace=0;
#     end
# end
# end
# return q_lowDiffDiff
# end

# # test code for diffdiffCheck
# inData=5 .*sin.(LinRange(0,2*pi,100))+rand(100)
# inData[4:13] .= [LinRange(0,1,10);]
# q_lowDiffDiff=diffdiffCheck(inData)
# hcat(inData,q_lowDiffDiff)'

function checkStagnant(inData,misfit_space,std_space;lowLimit_std=0.1,probabilityLimit=1e-2,analysis=false)
    # inData should be a vector
    validIn=filter(!isnan,inData)
    if isempty(validIn) # return if all inputs are nan
        q_curr,misfit,std_curr,pdf_curr = falses(size(inData)) , fill(NaN,size(inData)), fill(NaN,size(inData)), fill(NaN,size(inData))
    else
        diffIn=diff(inData)
        inRange = maximum(validIn) - minimum(validIn) 
        step_obs = minimum(filter(>(0),abs.(diffIn));init=Inf)

        if step_obs < inRange/64 
            q_lowDiff=diffCheck(inData; diffLimit= max(inRange/5000,step_obs) )
            q_lowDiffDiff=diffdiffCheck(inData;diffIn=diffIn, diffDiffLimit= max(inRange/5000,step_obs));

            # get the misft of stagant periods ------------
            misfit,seqNum=sequenceAbsMean(q_lowDiffDiff .| q_lowDiff,misfit_space)
            # get the std of stagnant periods -------
            std_curr=std_space./sqrt.(seqNum)
            std_curr[std_curr .< lowLimit_std] .=lowLimit_std # set lower limit for std. zero std might occur due to stagnant error
            # check the probability -------
            normalizedStagDistance = misfit./std_curr
            pdf_curr = pdf.(Normal(),normalizedStagDistance)
            q_curr=pdf_curr .< probabilityLimit
        else# ignore the stagnantCheck for observations with low resolutions
            q_curr,misfit,std_curr,pdf_curr = falses(size(inData)) , fill(NaN,size(inData)), fill(NaN,size(inData)), fill(NaN,size(inData))
        end
    end

    # output ------------
    if analysis 
        return q_curr,misfit,std_curr,pdf_curr
    else
        return q_curr
    end
end

"""
q_lowDiff,q_lowDiffDiff=checkStagnant(inData)

"""
function checkStagnant(inData)
    # inData should be a vector
    diffIn=diff(inData)

    # ignore the stagnantCheck for observations with low resolutions
    q_lowDiff=diffCheck(inData[:,n_column];diffIn=diffIn)
    q_lowDiffDiff=diffdiffCheck(inData[:,n_column];diffIn=diffIn);

    return q_lowDiff,q_lowDiffDiff
end

# # test code for stagnantCheck
# inData=5 .*sin.(LinRange(0,2*pi,100))+rand(100)
# inData[4:13] .= 0
# inData[end-12:end-3] .= [LinRange(0,1,10);]
# q_lowDiffDiff=stagnantCheck(inData)
# hcat(inData,q_lowDiffDiff)'


function sequenceAbsMean(q_seq,data2mean)
# out, q_seq, data2mean has the same size, and should be vectors
# 2017-1-12 add output of seqNum
# 2017-3-3 calculate mean absolute distance istead of mean distance
    out=fill(NaN,size(q_seq));
    seqNum=fill(NaN,size(q_seq));
    n_start=nothing;#n_end=NaN;n_step=1;
    for n_step=1:size(out,1)-1
        if q_seq[n_step]# into sequence
            if isnothing(n_start) 
                n_start=n_step; # mark start
            end
            if !q_seq[n_step+1] || n_step==size(out,1)-1 # not in sequence next time or end
                if n_step==size(out,1)-1 && q_seq[n_step+1]
                    n_end=n_step+1; # mark end to the 
                else
                    n_end=n_step; # mark end
                end
                ## do my stuff
                seqNum[n_start:n_end] .= sum(.!isnan.(data2mean[n_start:n_end]));
                out[n_start:n_end] .= nanmean(abs.(data2mean[n_start:n_end]));
                ## reset 
                n_start=nothing;# mark normal
            end
        end
    end
    return out,seqNum
end
