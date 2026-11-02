using Distances,LinearAlgebra
using utilities
# using StaDataTool

export OptimalInterpolation,getCovRatio
function OptimalInterpolation(simMeshData,observation,obs_mapper;   R_obs=nothing, CovRatio_HBHt=nothing,CharaLength=50e3,
    f_simErr=nothing,par=false)# (sim_x,sim_y) ->  0.3^2*(sim_x*sim_y))
    # f_simErr= (sim_x,sim_y) -> (0.1*mean(simMeshData)).^2 .+ 0.3^2*(sim_x*sim_y))
    
# write this form of input latter: optimalInterpolation(simMeshData,observation,obs_mapper,option)
# updates:
# 2020-03-11 change the order of input variables, add default calculation for R_obs,CovRatio_HBHt,CovRatio_BHt

    if isnothing(CovRatio_HBHt)
        # calculate this part outside this function and reuse CovRatio_HBHt
        # CovRatio_BHt, because they don't change for assimilating stations
        # in optimal interpolation
        lonV,latV=getNearestLonLat(obs_mapper);    #shift lon lat to grid center to mitigate unrealistic intorpolation for non obs areas
        obsLonLat=hcat(lonV,latV)
        # get CovRatio_HBHt
        d_HBHt=pairwise(Haversine(),obsLonLat,dims=1);
        CovRatio_HBHt=getCovRatio(d_HBHt,CharaLength);
    end

    q_validObs=.!isnan.(observation);
    # get R_obs
    if isnothing(R_obs)
        # not input
        R_obs=diagm(ones(sum(q_validObs)).*(0.1*mean(observation[q_validObs]))^2);  
    else
        # get the valid input
        R_obs=R_obs[q_validObs,q_validObs];
    end

    # set default f_simErr
    backgroundVar=(0.1*mean(simMeshData)).^2 
    function f_simErr_default(sim_x,sim_y)
        backgroundVar .+ 0.3^2*(sim_x*sim_y)
    end
    if isnothing(f_simErr)
        f_simErr=f_simErr_default
    end

    # get staSim
    staSim=map2sta(obs_mapper,simMeshData);
    staSim=staSim[q_validObs]
    # get Cov_HBHt
    priVar=f_simErr(staSim,staSim')
    Cov_HBHt=CovRatio_HBHt[q_validObs,q_validObs].*priVar;
    Ky=((Cov_HBHt+R_obs)^-1)*(observation[q_validObs] -staSim);
    lonV,latV=getNearestLonLat(obs_mapper);    #shift lon lat to grid center to mitigate unrealistic intorpolation for non obs areas
    obsLonLat=hcat(lonV[q_validObs],latV[q_validObs])
    out=fill(NaN,size(simMeshData))
    if par
    Threads.@threads for n_grid in eachindex(out)
            out[n_grid]=OI_aGrid(obs_mapper.lonM[n_grid],obs_mapper.latM[n_grid],simMeshData[n_grid],obsLonLat,CharaLength,staSim,Ky,f_simErr)
        end
    else
        for n_grid in eachindex(out)
            out[n_grid]=OI_aGrid(obs_mapper.lonM[n_grid],obs_mapper.latM[n_grid],simMeshData[n_grid],obsLonLat,CharaLength,staSim,Ky,f_simErr)
        end
    end
    return out
end
function OI_aGrid(lon_grid,lat_grid,state_grid,obsLonLat,CharaLength,staSim,Ky,f_simErr)
    normalizedDist=colwise(Haversine(),obsLonLat',[lon_grid;lat_grid])'./CharaLength;
    Cov_BHt=@. (1+normalizedDist)*exp(-normalizedDist) * f_simErr(state_grid,staSim')
    return state_grid + (Cov_BHt*Ky)[1];
end
function getCovRatio(d_Matrix,CharaLength)
    #   get the coverance ratio
        normalizedDist=d_Matrix./CharaLength;
        CovRatio=@. (1+normalizedDist)*exp(-normalizedDist);
        return CovRatio
end
        
