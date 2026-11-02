
# using Statistics,LocalizerM,Infiltrator,LinearAlgebra,MKL
# using MKL will speed this up, but will cause pycall error for julia 1.7
# see the following link for detail: https://github.com/JuliaLinearAlgebra/MKL.jl/issues/97

# updates
# 2022-11-08 add support for obsErr vector. Sample usage: EnKF(State,observation,staSim,localizer,space_cutoff_all;obsErr=myObsErr)
using Statistics,LocalizerM,LinearAlgebra,JLD2
using MKL
export EnKF
function EnKF(State,observation,staSim,localizer,space_cutoff_all;
    EnKFtype="DEnKF_locCorr",obsErrFun=obs-> 0.1*obs,    obsErr=nothing,timeCutoff=nothing,inflateFactor=nothing,nGrids2da=nothing,par=false)
    # # dimensions of input State(grid,member), observation(station,time), staSim(station,member,time)
    # # localizer is a class implement getLocCoef method, which may be created by localizerM
    # # staSim is HA, simulation after projection
    # Δ is the difference between members the their mean
    if isnothing(inflateFactor);spreadInflateFactor=nothing;else;spreadInflateFactor=sqrt.(inflateFactor);end
    if isnothing(obsErr);obsErr=obsErrFun(observation);end
    if !(space_cutoff_all isa Vector);           space_cutoff_all=fill(space_cutoff_all,size(State,1));       end
    if in(EnKFtype,["DEnKF_locCorr","EnSRKF_locCorr"]);   useCoefMatrix=true;   else;    useCoefMatrix=false;     end
    # # prepare variables (State,observation,staSimΔ,obsSubSim)
    # # *********************************************************************
    # # reduce station + time to a single dimension for algorithm consistency with 3D EnKF
    # # inflate state
    if !isnothing(spreadInflateFactor);State=inflateMembers(State,spreadInflateFactor);end
    observation=observation[:];
    memberNumber=size(State,2);
    if ndims(staSim)>2
        staSim=reshape(permutedims(staSim,(1,3,2)),Int(length(staSim)/memberNumber),memberNumber);#doest not used directly in main iteration
    end
    staSimΔ=staSim .- mean(staSim,dims=2);
    
    # # preparation for different methods
    if EnKFtype in ["DEnKF","DEnKF_inflate","DEnKF_locCorr","ETKF","LETKF","EnSRKF","EnSRKF_locCorr"]
        obsSubSim=observation-mean(staSim,dims=2);## one dimension, a vector
    end
    # # main iteration update the state for each grid/row/augmentstate
    # # *********************************************************************
    function DA_aGrid(State_aGrid,observation,staSimΔ,obsSubSim,obsErr,  localizer,nGrid_curr,space_cutoff_curr,timeCutoff,useCoefMatrix,spreadInflateFactor) 
        if useCoefMatrix>0
            q_loc,loc_coef,loc_coefMatrix=getLocCoef(localizer,nGrid_curr,space_cutoff_curr,time_cutoff=timeCutoff,nargout=3);
        else
            q_loc,loc_coef=getLocCoef(localizer,nGrid_curr,space_cutoff_curr,time_cutoff=timeCutoff,nargout=2);
        end
        # println(nGrid_curr,"  ",length(nGrids2da)," staNear:",length(loc_coef))# comment out later
        if length(loc_coef)>0
            # # deal with nan observation -------------
            q_loc=vec(q_loc);## reduce station + time to a single dimension for algorithm consistency with 3D EnKF
            loc_coef=vec(loc_coef);
            q_validObs=.!isnan.(observation[q_loc]);
            if sum(q_validObs)==length(q_validObs) #no nan
            elseif sum(q_validObs)==0# all nan
                return fill(NaN,size(State_aGrid))
            else # some nan observation
                loc_coef=loc_coef[q_validObs];## get available only
                if useCoefMatrix>0  ## clear corresponding loc_coefMatrix
                    loc_coefMatrix=loc_coefMatrix[q_validObs,q_validObs]
                end
                q_loc[q_loc]=q_validObs;## clear corresponding q_loc
            end
            # # inflation for staSim--------------------
            if isnothing(spreadInflateFactor)
                locstaSimΔ=staSimΔ[q_loc,:];
            else
                locstaSimΔ=staSimΔ[q_loc,:].*spreadInflateFactor[nGrid_curr]
            end
                
            # # ********************* EnKF Assimilation ************
            if  EnKFtype == "DEnKF_locCorr"
                updatedState_aGrid=DEnKF_locCorr(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc],loc_coef,loc_coefMatrix,obsErr[q_loc]);
            elseif EnKFtype == "original"
                updatedState_aGrid=originalEnKF(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc,:],loc_coef);
            elseif EnKFtype == "DEnKF"
                updatedState_aGrid=DEnKF(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc],loc_coef,obsErr[q_loc]);
            elseif EnKFtype == "EnKF_lowRank"
                updatedState_aGrid=EnKF_lowRank(State_aGrid,obsBias[q_loc,:],locstaSimΔ,obsSubSim[q_loc,:],loc_coef);
            elseif EnKFtype == "ETKF"
                updatedState_aGrid=ETKF(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc],obsErr[q_loc]);
            elseif EnKFtype == "LETKF"
                updatedState_aGrid=LETKF(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc],loc_coef,obsErr[q_loc]);
            elseif EnKFtype == "EnSRKF"
                updatedState_aGrid=EnSRKF(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc],loc_coef,obsErr[q_loc])
            elseif EnKFtype == "EnSRKF_locCorr"
                updatedState_aGrid=EnSRKF_locCorr(State_aGrid,observation[q_loc],locstaSimΔ,obsSubSim[q_loc],loc_coef,loc_coefMatrix,obsErr[q_loc])
            else
                @error "EnKFtype of $EnKFtype is not supported yet"
            end
            # for debug
            # if mean(updatedState_aGrid)<0;       println("nGrid:",nGrid_curr," nValid:" ,length(loc_coef), " updated:",mean(updatedState_aGrid));        end
        else
            updatedState_aGrid=similar(State_aGrid)
            updatedState_aGrid .= NaN
        end
        return updatedState_aGrid
    end
    
    updatedState=fill(NaN,size(State));# allocate for output
    if isnothing(nGrids2da);        nGrids2da=1:size(updatedState,1);    end
    if par
        Threads.@threads for nGrid_curr in axes(updatedState,1)
            if nGrid_curr in nGrids2da
                space_cutoff_curr=space_cutoff_all[nGrid_curr]
                State_aGrid=State[nGrid_curr,:]';
                updatedState[nGrid_curr,:]=DA_aGrid(State_aGrid,observation,staSimΔ,obsSubSim,obsErr,localizer,
                    nGrid_curr,space_cutoff_curr,timeCutoff,useCoefMatrix,spreadInflateFactor);
            end
        end
    else
        for n in eachindex(nGrids2da)
            nGrid_curr=nGrids2da[n]
            space_cutoff_curr=space_cutoff_all[nGrid_curr]
            State_aGrid=State[nGrid_curr,:]';
            updatedState[nGrid_curr,:]=DA_aGrid(State_aGrid,observation,staSimΔ,obsSubSim,obsErr,localizer,
                nGrid_curr,space_cutoff_curr,timeCutoff,useCoefMatrix,spreadInflateFactor);
        end
    end
    return updatedState
end
# # different method of EnKF
# # *************************************************************************

# # base on determistic EnKF, but localizing in correlation
function DEnKF_locCorr(State,loc_observation,loc_staSimΔ,loc_obsSubSim,loc_coef,loc_coefMatrix,loc_obsErr)
    # # do not support multiple states, size(State,1) must equal to 1
    stateMean=mean(State);## specify dimension will speed up a little
    StateΔ=State.-stateMean;
	denominator=(loc_coefMatrix.*(loc_staSimΔ*loc_staSimΔ') .+ (size(State,2)-1)*diagm(loc_obsErr.^2));
    gainMatrix=loc_coef'.*(StateΔ*loc_staSimΔ')/denominator;
	updatedStateMean=stateMean.+gainMatrix*loc_obsSubSim;
	# updatedStateMean=stateMean.+loc_coef'.*(StateΔ*loc_staSimΔ')/(loc_coefMatrix.*(loc_staSimΔ*loc_staSimΔ') .+ (size(State,2)-1)*diagm(loc_obsErr.^2))*loc_obsSubSim;
    updatedState=updatedStateMean.+(StateΔ.-gainMatrix*loc_staSimΔ/2);
    return updatedState
end

# # originalEnKF
function originalEnKF(State,loc_observation,loc_staSimΔ,loc_obsSubSim,loc_coef)
    StateΔ=State.-mean(State);
    denominator=(loc_staSimΔ*loc_staSimΔ'.+(size(State,2)-1)*diagm((0.1* loc_observation ).^2));
    updatedState=State.+StateΔ*loc_staSimΔ'/denominator*diagm(loc_coef)*loc_obsSubSim;
    return updatedState
end

# # determistic EnKF
function DEnKF(State,loc_observation,loc_staSimΔ,loc_obsSubSim,loc_coef,loc_obsErr)
    # # do not support multiple states, size(State,1) must equal to 1
    stateMean=mean(State);## specify dimension will speed up a little
    StateΔ=State.-stateMean;
    denominator=(loc_staSimΔ*loc_staSimΔ'.+(size(State,2)-1)*diagm(loc_obsErr.^2));
    gainMatrix=StateΔ*loc_staSimΔ'/denominator*diagm(loc_coef);
    updatedStateMean=stateMean.+gainMatrix*loc_obsSubSim;
    updatedState=updatedStateMean.+(StateΔ.-gainMatrix*loc_staSimΔ/2);
    # @save "debug.jld" State loc_observation loc_staSimΔ loc_obsSubSim loc_coef loc_obsErr gainMatrix updatedStateMean
    return updatedState
end

function inflateMembers(inData,inflateFactor)
    # # inflate in with inflateFactor, the last dimension of in should be ensembles
    # # the inflateFactor here should be sqrt(correlation inflate factor)
    sizeIn=size(inData);
    memberNumber=sizeIn[end];
    rspIn=reshape(inData,:,memberNumber);
    rspIn_mean=mean(rspIn,dims=2);
    rspOut=rspIn_mean .+ (rspIn .- rspIn_mean) .* inflateFactor
    outData=reshape(rspOut,sizeIn);
    return outData
end

# # square root low-rank EnKF
function EnSRKF(State,loc_observation,loc_staSimΔ,loc_obsSubSim,loc_coef,loc_obsErr)
    # # algorithm in page 4 (section 3.2) of 'Sampling strategies and square root analysis schemes for the EnKF'
    # ΔA is A' in the paper
    N=size(State,2)
    Ψᶠ_mean = dropdims(mean(State,dims=2),dims=2)
    ΔA= State .- Ψᶠ_mean
    S=loc_staSimΔ
    R=Diagonal(loc_obsErr.^2)

    # update the mean ------------------------------------------------------------
    C=S * S' + (N-1) * R
    Λ,Z=eigen(C)    
    # ignore negative eigenvalues since the diagonal values in the symmetric matrix C are positive
    q_valid = Λ.>0
    Λ,Z = Λ[q_valid] , Z[:,q_valid]

    Λ=Diagonal(Λ)
    Ψᵃ_mean = Ψᶠ_mean +  ΔA * S' * Z * Λ^(-1) * Z' * (loc_obsSubSim .* loc_coef)
    # update the spread ----------------------------------------------------------
    X₂ = Λ^(-1/2) * Z' * S
    U₂,Σ₂,V₂=svd(X₂;full=true)
    
    # set upper limit of Σ₂ to 1, other wise error will occur when calculating the square root
    Σ₂[Σ₂ .> 1] .= 1

    if length(Σ₂)<N;append!(Σ₂,zeros(N-length(Σ₂)));end
    Σ₂ = Diagonal(Σ₂)
    I=Diagonal(ones(N))
    # Θ,~,~ = svd(rand(N,N)) # do not use this for efficiency, consider obtain this from a input variable
    # Θ = I;
    ΔAᵃ = ΔA *  V₂ * (I - Σ₂' * Σ₂  )^(1/2)

    State_assimilated = Ψᵃ_mean .+ ΔAᵃ 
    return State_assimilated
end
# # square root low-rank EnKF using coveriance locolization
function EnSRKF_locCorr(State,loc_observation,loc_staSimΔ,loc_obsSubSim,loc_coef,loc_coefMatrix,loc_obsErr)
    # # algorithm in page 4 (section 3.2) of 'Sampling strategies and square root analysis schemes for the EnKF'
    # ΔA is A' in the paper
    N=size(State,2)
    Ψᶠ_mean = dropdims(mean(State,dims=2),dims=2)
    ΔA= State .- Ψᶠ_mean
    S=loc_staSimΔ
    R=Diagonal(loc_obsErr.^2)
    # C=(S * S') + (N-1) * R
    C=Symmetric(loc_coefMatrix .* (S * S') + (N-1) * R)
    Λ,Z=eigen(C)

    # ignore negative eigenvalues since the diagonal values in the symmetric matrix C are positive
    q_valid = Λ.>0
    Λ,Z = Λ[q_valid] , Z[:,q_valid]

    Λ=Diagonal(Λ)
    # update the mean 
    Ψᵃ_mean = Ψᶠ_mean +  loc_coef' .* (ΔA * S' ) * Z * Λ^(-1) * Z' * loc_obsSubSim
    
    # update the spread
    X₂ = Λ^(-1/2) * Z' * S
    U₂,Σ₂,V₂=svd(X₂;full=true)

    # set upper limit of Σ₂ to 1, other wise error will occur when calculating the square root
    Σ₂[Σ₂ .> 1] .= 1

    if length(Σ₂)<N;append!(Σ₂,zeros(N-length(Σ₂)));end
    Σ₂ = Diagonal(Σ₂)
    I=Diagonal(ones(N))
    # Θ = I;
    # negtive eigenvalues for 
    ΔAᵃ = ΔA *  V₂ * sqrt(I - Σ₂' * Σ₂)

    State_assimilated = Ψᵃ_mean .+ ΔAᵃ 
    return State_assimilated
end
# function EnSRKF_lowRank(State,loc_observation,loc_staSimΔ,loc_obsSubSim,loc_coef,loc_obsErr)
# # # algorithm in page 16 of 'Sampling strategies and square root analysis schemes for the EnKF'
# #     # do not support multiple states, size(State,1) must equal to 1
# #     #loc_staSimΔ for S;  loc_obsSubSim for (d-mean(HPhi,2))
# #     # loc_obsΔ for E;
# #     # use bsxfun to perform diagonal matrix production
#     stateMean=mean(State,dims=2);
#     StateΔ=State-stateMean;
#     memberNumber=size(State,2);
#     if size(loc_observation,2)==memberNumber
#         loc_obsΔ=loc_observation;
#     elseif size(loc_observation,2)==1
#         loc_obsΔ=getObsBias(loc_observation,memberNumber);
#         loc_obsΔ=0.1*loc_observation .* compressSample(randn(numel(loc_observation),memberNumber*20),memberNumber)
#     else
#         error("size(loc_observation,2) should be size(State,2) or 1")
#     end
    
#     U0,Sigma0,~=svd(loc_staSimΔ);
#     if size(loc_staSimΔ,1)==1
#         diag_Sigma0=Sigma0[1];
#     else
#         diag_Sigma0=diag(Sigma0);
#     end
#     q_valid_diag_Sigma0=diag_Sigma0 > sum(diag_Sigma0)*0.00001;## avoid inverse of small value
#     if length(q_valid_diag_Sigma0)>=memberNumber
#         q_valid_diag_Sigma0[end]=0;## only use N-1 element, said by the paper, don't know why
#             # diag_Sigma0 has N element while diag_Sigma0plus less than N-1
#     end
#     diag_Sigma0plus=1 ./ diag_Sigma0[q_valid_diag_Sigma0];
#     # X0=Sigma0plus*U0'*E
#     X0=zeros(memberNumber);
#     X0[q_valid_diag_Sigma0,:]=(diag_Sigma0plus .* U0[:,q_valid_diag_Sigma0]')   * loc_obsΔ;# N*N
#         # X0[q_valid_diag_Sigma0,:]=( bsxfun(@times,diag_Sigma0plus,U0(:,q_valid_diag_Sigma0)') )... 
#         #     *loc_obsΔ;# N*N
#     # svd for X0 
#     [U1,Sigma1,~]=svd(X0);# all N*N
#     # X1=U0*Sigma0plus*U1  # X1 is m*N
#     X1=U0[:,q_valid_diag_Sigma0] *  (diag_Sigma0plus .* U1[q_valid_diag_Sigma0,:]);
#     # update state mean
#     diag_dominator=1+diag(Sigma1).^2;#diag(I+Sigma1^2) N*N
# 	updatedStateMean=stateMean+StateΔ * loc_staSimΔ' * X1 * (1 ./(diag_dominator) .* X1') * (loc_coef[:].*loc_obsSubSim);
#     X2=( diag_dominator.^-1/2 .* X1') * loc_staSimΔ;
#     # svd for X2
#     [~,Sigma2,V2]=svd(X2);# all N*N
#     updatedState_aGrid=updatedStateMean+StateΔ * ( V2 .* sqrt(1-diag(Sigma2)'.^2) );
#     return updatedState_aGrid
# end
