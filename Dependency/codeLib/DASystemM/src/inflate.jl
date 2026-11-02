
using LinearAlgebra
export inflate
"""

passed: WB_locCorr
need to check: MLE_state", "MLE_state_locCorr", "WB", "WB_loc",  "WB_locB
"""
function inflate(State, observation, obsSim, localizer,space_cutoff_all; inflateType = "WB_locCorr",obsErrFun= obs-> 0.1 .* obs,timeCutoff=nothing,
    nGrid2inflate = 1:size(State, 1), lowerBound_varS2varE = 0.4^2 ,upperBound_varS2varE = 20^2)
    # dimensions of input State(grid,member), observation(station), obsSim(station,member)
    # localizer is a class implement getLocCoef method, which may be created by localizer4D
    # obsSim is HA, simulation after projection

    # prepare variables (State,observation,obsSimBias,obsSubSim)
    # *********************************************************************
    memberNumber = size(obsSim, 2)
    updatedState = fill!(similar(State), NaN)  # allocate for output
    obsSimMean = mean(obsSim, dims=2)
    obsSimBias = obsSim .- obsSimMean
    obsSubSim = observation .- obsSimMean
    if !(space_cutoff_all isa Vector);           space_cutoff_all=fill(space_cutoff_all,size(State,1));       end


    # preparation for different methods *********
    if any(inflateType .== ["MLE_state", "MLE_state_locCorr", "WB", "WB_loc", "WB_locCorr", "WB_locB"])
        inflateFactor = ones(size(State, 1))
    elseif inflateType == "MLE_state_obs"
        inflateFactor = ones(size(State, 1), 2)
    end
    if any(inflateType .== ["MLE_state_locCorr", "WB_locCorr"])
        useCoefMatrix = 1
    else
        useCoefMatrix = 0
    end

    # main iteration, calculate the inflation factor for each state
    # *********************************************************************
    # Threads.@threads     for nGrid_curr in nGrid2inflate
    for nGrid_curr in nGrid2inflate
        space_cutoff_curr=space_cutoff_all[nGrid_curr]
        # disp(nGrid_curr)% comment out later
        # get localization coeficient++++++++++
        if !isnothing(localizer)
            if useCoefMatrix > 0
                q_loc, loc_coef, loc_coefMatrix = getLocCoef(localizer,nGrid_curr, space_cutoff_curr, time_cutoff=timeCutoff,nargout=3)
            else
                q_loc, loc_coef = getLocCoef(localizer,nGrid_curr, space_cutoff_curr, time_cutoff=timeCutoff,nargout=2)
            end
        else
            q_loc = trues(size(observation))
            loc_coef = ones(size(observation))
            loc_coefMatrix = diagm(loc_coef)
        end
        if length(loc_coef)>0
            # # deal with nan observation -------------
            q_loc=vec(q_loc);## reduce station + time to a single dimension for algorithm consistency with 3D EnKF
            loc_coef=vec(loc_coef);
            q_validObs=.!isnan.(observation[q_loc]);
            if sum(q_validObs)==length(q_validObs) #no nan
            elseif sum(q_validObs)==0# all nan
                continue
            else
                loc_coef=loc_coef[q_validObs];## get available only
                if useCoefMatrix > 0
                    loc_coefMatrix=loc_coefMatrix[q_validObs,q_validObs]
                end
                q_loc[q_loc]=q_validObs;## clear corresponding q_loc
            end
            # calculate inflation factore+++++++++++
            loc_observation = observation[q_loc]
            loc_obsSimBias = obsSimBias[q_loc, :]
            loc_obsSubSim = obsSubSim[q_loc]
            loc_obsErr = obsErrFun(loc_observation)
            # loc_obsErr= diag(obsErrFun(loc_observation ).^2);

            # calculate the bound_factor
            # varE2varS = mean(loc_obsErr.^2) / mean(var(loc_obsSimBias, dims=2))
            varE2varS = mean(loc_obsErr.^2) / var(State[nGrid_curr,:])
            # varE2varS=mean(loc_obsErr.^2)/mean(var(obsSimBias,0,2));
            lowerBound_factor = lowerBound_varS2varE * varE2varS
            upperBound_factor = upperBound_varS2varE * varE2varS
            if inflateType == "WB"
                inflateFactor[nGrid_curr] = WBinflation(loc_obsErr, loc_obsSimBias, loc_obsSubSim)
            elseif inflateType == "WB_loc"
                inflateFactor[nGrid_curr] = WBinflation(loc_obsErr, loc_obsSimBias, loc_obsSubSim .* loc_coef)
            elseif inflateType == "WB_locCorr"
                inflateFactor[nGrid_curr] = WBinflation_locCorr(loc_obsErr, loc_obsSimBias, loc_obsSubSim, loc_coef, loc_coefMatrix)
            elseif inflateType == "WB_locB"
                inflateFactor[nGrid_curr] = WBinflation_locB(loc_obsErr, loc_obsSimBias, loc_obsSubSim, loc_coef, option)
            elseif inflateType == "MLE_state"
                inflateFactor[nGrid_curr] = maximumLikelihoodInflationEstimate(loc_obsErr, loc_obsSimBias, loc_obsSubSim, option)
            elseif inflateType == "MLE_state_locCorr"
                inflateFactor[nGrid_curr] = maximumLikelihoodInflationEstimate_locCorr(loc_obsErr, loc_obsSimBias, loc_obsSubSim, loc_coef, loc_coefMatrix, option)
            elseif inflateType == "MLE_state_obs"
                inflateFactor[nGrid_curr, 1], inflateFactor[nGrid_curr, 2] = maximumLikelihoodInflationEstimate(loc_obsErr, loc_obsSimBias, loc_obsSubSim, option)
            else
                error("option.inflateType should be one of: MLE_state,MLE_state_obs, but", option.EnKFtype)
            end
            # check the bound_factor
            if inflateType in ["WB","WB_loc","WB_locB","WB_locCorr"]
                if inflateFactor[nGrid_curr] <lowerBound_factor
                    inflateFactor[nGrid_curr]=lowerBound_factor  
                elseif inflateFactor[nGrid_curr]>upperBound_factor
                    inflateFactor[nGrid_curr]=upperBound_factor
                end      
            end

        end
    end
    # figure;mesh(reshape(sqrt(inflateFactor),432,339)');
    # inflate the state
    # *********************************************************************
    return inflateFactor
    # updatedState = inflateMembers(State, sqrt.(inflateFactor[:, 1]))
    # return inflateFactor, updatedState
end
function maximumLikelihoodInflationEstimate(obsErr, obsSimBias, obsSubSim, option)
    # R is the observation error covariance
    # algorithm based on "Maximum likelihood estimation of inflation factors on error covariance matrices for ensemble Kalman filter assimilation" page 3
    if haskey(option, "fmOption")
        fmOption = option.fmOption
    else
        fmOption = Optim.Options(g_tol=1e-2, x_tol=1e-2, f_tol=1e-2, display=:none)
    end
    R = diagm(obsErr.^2)
    HPHt = (obsSimBias * obsSimBias') / (size(obsSimBias, 2) - 1)
    if nargout() == 1 # inflate the state only
        fun = lamda -> log(det(lamda * HPHt + R)) + obsSubSim' / (lamda * HPHt + R) * obsSubSim
        stateInflate = optimize(fun, 1, option.lowerBound_factor, option.upperBound_factor, method = :constrained_NelderMead, options = fmOption).minimizer
    else # inflate state and obsError
        fun = lamda -> log(det(lamda[1] * HPHt + lamda[2] * R)) + obsSubSim' / (lamda[1] * HPHt + lamda[2] * R) * obsSubSim
        lamda = optimize(fun, [1, 1], [option.lowerBound_factor, option.lowerBound_factor], [option.upperBound_factor, option.upperBound_factor], method = :constrained_NelderMead, options = fmOption).minimizer
        stateInflate = lamda[1]
        obsInflate = lamda[2]
    end
    return stateInflate, obsInflate
end

function maximumLikelihoodInflationEstimate_locCorr(obsErr, obsSimBias, obsSubSim, loc_coef, loc_coefMatrix, option)
    # R is the observation error covariance
    # algorithm based on "Maximum likelihood estimation of inflation factors on error covariance matrices for ensemble Kalman filter assimilation" page 3
    if haskey(option, "fmOption")
        fmOption = option.fmOption
    else
        fmOption = Optim.Options(g_tol=1e-2, x_tol=1e-2, f_tol=1e-2, display=:none)
    end
    R = diagm(obsErr.^2)
    obsSubSim = sqrt.(loc_coef) .* obsSubSim
    HPHt = loc_coefMatrix .* (obsSimBias * obsSimBias') / (size(obsSimBias, 2) - 1)
    if nargout() == 1 # inflate the state only
        fun = lamda -> log(det(lamda * HPHt + R)) + obsSubSim' / (lamda * HPHt + R) * obsSubSim
        stateInflate = optimize(fun, 1, option.lowerBound_factor, option.upperBound_factor, method = :constrained_NelderMead, options = fmOption).minimizer
    else # inflate state and obsError
        fun = lamda -> log(det(lamda[1] * HPHt + lamda[2] * R)) + obsSubSim' / (lamda[1] * HPHt + lamda[2] * R) * obsSubSim
        lamda = optimize(fun, [1, 1], [option.lowerBound_factor, option.lowerBound_factor], [option.upperBound_factor, option.upperBound_factor], method = :constrained_NelderMead, options = fmOption).minimizer
        stateInflate = lamda[1]
        obsInflate = lamda[2]
    end
    return stateInflate, obsInflate
end
function WBinflation(obsErr, obsSimBias, obsSubSim)
    # obsSubsim must no contain nan
    # the algorithm here follows "Maximum likelihood estimation of inflation factors on error covariance matrices for ensemble Kalman filter assimilation" page 3
    # locolization coefcient is no considered
    diagRInv = 1 ./ obsErr
    HPHt = (obsSimBias * obsSimBias') / (size(obsSimBias, 2) - 1)
    inflateFactor = ((diagRInv .* obsSubSim)' * (diagRInv .* obsSubSim) - length(obsSubSim)) / tr((diagRInv * diagRInv') .* HPHt)
    return inflateFactor
end

function WBinflation_locCorr(obsErr, obsSimBias, obsSubSim, loc_coef, loc_coefMatrix)
    # WBinflation with B type localization
    # obsSubsim must no contain nan
    # the algorithm here follows "Maximum likelihood estimation of inflation factors on error covariance matrices for ensemble Kalman filter assimilation" page 3
    # locolization coefcient is no considered
    obsSubSim = sqrt.(loc_coef) .* obsSubSim
    diagRInv = 1 ./ obsErr
    HPHt = loc_coefMatrix .* (obsSimBias * obsSimBias') / (size(obsSimBias, 2) - 1)
    inflateFactor=((diagRInv.*obsSubSim)'  * (diagRInv.*obsSubSim) - sum(sqrt.(loc_coef)) )   / tr((diagRInv*diagRInv').*HPHt);
    return inflateFactor
end

function WBinflation_locB(obsErr, obsSimBias, obsSubSim, loc_coef, option)
    # WBinflation with B type localization
    # obsSubsim must no contain nan
    # the algorithm here follows "Maximum likelihood estimation of inflation factors on error covariance matrices for ensemble Kalman filter assimilation" page 3
    # locolization coefcient is no considered
    # obsSubSim=sqrt(loc_coef).*obsSubSim;
    diagRInv = 1 ./ obsErr
    HPHt = (obsSimBias * obsSimBias') / (size(obsSimBias, 2) - 1)
    inflateFactor = ((diagRInv .* obsSubSim)' * (loc_coef .* diagRInv .* obsSubSim) - sum(loc_coef)) / (diag(HPHt) ./ obsErr.^2) * loc_coef
    if inflateFactor < option.lowerBound_factor
        inflateFactor = option.lowerBound_factor
    end
    if inflateFactor > option.upperBound_factor
        inflateFactor = option.upperBound_factor
    end
    return inflateFactor
end