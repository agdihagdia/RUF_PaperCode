module LocalizerM
using Base,Interpolations,Distances,utilities
export localizerM,getLocCoef,getVaringSpaceCutOff,getNearest,LocalizerStr
mutable struct LocalizerStr
    stateLonLat#::Matrix{Float64}
    obsLonLat#::Matrix{Float64}
    stateTime
    obsTime # must be a row for matix multiplication in getLocCoef
    locFun::ScaledInterpolation
    spaceCoef_cutoff
    coefMatrix
end

# ======= modification =============


function getNearest(obj::LocalizerStr,n_grid,nOut=1;distOut=false)
    # spaceDistance=Haversine().([obj.stateLonLat[n_grid,:]],obj.obsLonLat)
    spaceDistance=geoDist(obj.stateLonLat[n_grid,:],obj.obsLonLat)
    ind=partialsortperm(spaceDistance, nOut)
    if distOut; 
        return ind,spaceDistance[ind];
    else
        return ind
    end
end

function Base.copy(loc::LocalizerStr)
    return LocalizerStr(loc.stateLonLat,loc.obsLonLat,loc.stateTime,loc.obsTime,loc.locFun,loc.spaceCoef_cutoff,loc.coefMatrix)
end
function localizerM(stateLonLat,obsLonLat;stateTime=nothing,obsTime=nothing,ini_coefMatrix=false,spaceCoef_cutoff=nothing,locType="Gaspari_Cohn")
    # obsLonLat=collect(eachrow(obsLonLat))
    # prepare locFun ---------------------------------------------------
    if locType=="Gaspari_Cohn"
        GS_y=[1,0.994610462525819,0.978978819934086,0.953992031959236,0.920630623448450,0.879949799105578,0.833060558235049,0.781110809485783,0.725266485595112,0.666692658132689,0.606534652244403,0.545899161396294,0.485835362118470,0.427316028749016,0.371218648177911,0.318306534590944,0.269209944213627,0.224407190055107,0.184205810337121,0.148729919113402,0.117939471697685];
        # GS_y=Float16.([1,0.994610462525819,0.978978819934086,0.953992031959236,0.920630623448450,0.879949799105578,0.833060558235049,0.781110809485783,0.725266485595112,0.666692658132689,0.606534652244403,0.545899161396294,0.485835362118470,0.427316028749016,0.371218648177911,0.318306534590944,0.269209944213627,0.224407190055107,0.184205810337121,0.148729919113402,0.117939471697685]);
        GS_x=LinRange(0,2,length(GS_y));
        itp = interpolate(GS_y, BSpline(Cubic(Line(OnGrid()))))
        itp = extrapolate(itp,0)   # g
        locFun = scale(itp,GS_x)
    elseif locType=="speLine"
        GS_y=Float32.([1, 1, 2/3, 1/3 , 0]);
        GS_x=LinRange(0,2,length(GS_y));
        itp = interpolate(GS_y, BSpline(Linear()))
        itp = extrapolate(itp,0)   # g
        locFun = scale(itp,GS_x)
    elseif locType=="Slide"
        GS_y=Float32.([1, 1, 0]);
        GS_x=LinRange(0,2,length(GS_y));
        itp = interpolate(GS_y, BSpline(Linear()))
        itp = extrapolate(itp,0)   # g
        locFun = scale(itp,GS_x)
    elseif locType=="Step"
        GS_y=Float32.([1,1]);
        GS_x=LinRange(0,2,length(GS_y));
        itp = interpolate(GS_y, BSpline(Constant()))
        itp = extrapolate(itp,0)   # g
        locFun = scale(itp,GS_x)

    # else
    #     locFun=@(normalizedDistance) calc_loccoeffs(1,option.locType , normalizedDistance);
    end
    # prepare coefMatrix to prevent repetitive computation----------------------
    if ini_coefMatrix
        # dist=pairwise(Haversine(),obsLonLat,dims=1);
        dist=geoDist(obsLonLat);
        coefMatrix=@. locFun(dist/spaceCoef_cutoff);
    else
        coefMatrix=nothing
    end
    return LocalizerStr(stateLonLat,obsLonLat,stateTime,obsTime,locFun,spaceCoef_cutoff,coefMatrix)
end

"""

cutoff_spatial = getVaringSpaceCutOff(localizer;validStaNum=1)
"""
function getVaringSpaceCutOff(obj::LocalizerStr;validStaNum=1,par=false)
    function getCutOff(stateLonLat,obsLonLat,validStaNum)
        spaceDistance=geoDist(stateLonLat,obsLonLat)
        if 1==validStaNum
            cutoff=minimum(spaceDistance)
        else
            cutoff=spaceDistance[partialsortperm(spaceDistance,validStaNum)]
        end
        return cutoff
    end
    cutoff_out=fill(NaN,size(obj.stateLonLat,1))
    if par
        Threads.@threads for n_grid in eachindex(cutoff_out)
            cutoff_out[n_grid]=getCutOff(obj.stateLonLat[n_grid,:],obj.obsLonLat,validStaNum)
        end
    else
        for n_grid in eachindex(cutoff_out)
            cutoff_out[n_grid]=getCutOff(obj.stateLonLat[n_grid,:],obj.obsLonLat,validStaNum)
        end
    end
    return cutoff_out
end
function getLocCoef(obj::LocalizerStr,n_grid,space_cutoff;time_cutoff=nothing,nargout=1)
    # % get the coef for localization.
    # % time should be in the format of datenum, set time_cutoff=inf or empty will use all time
    if isnothing(obj.spaceCoef_cutoff)
        spaceCoef_cutoff=space_cutoff
    else
        spaceCoef_cutoff=obj.spaceCoef_cutoff
    end

    # for q_loc ---------------------------------------------------
    # spaceDistance=Haversine().([view(obj.stateLonLat,n_grid,:)],obj.obsLonLat)
    spaceDistance =geoDist(obj.stateLonLat[n_grid,:],obj.obsLonLat)
    # spaceDistance=Haversine().([obj.stateLonLat[n_grid,:]],collect(eachrow(obj.obsLonLat)))
    # q_space=BitArray(undef,size(spaceDistance))
    q_space = spaceDistance .<= 2*space_cutoff;  
    # q_space = spaceDistance #.<= 2*space_cutoff;  
    q_loc=q_space 
    if nargout==1; return q_loc;end

    # for loc_coef------------------------------------
    loc_dist=spaceDistance[q_space]
    if !isempty(loc_dist)
        loc_coef = obj.locFun.(loc_dist/ spaceCoef_cutoff)
    else
        loc_coef=Vector{eltype(spaceDistance)}();
    end
    if nargout==2; return q_loc,loc_coef;end
    
    # for loc_coefMatrix ------------------------------
    if length(loc_coef)==0
        loc_coefMatrix=Array{Float64}(undef, 0, 0)
    elseif length(loc_coef)==1
        loc_coefMatrix=1;
    else
        if !isnothing(obj.coefMatrix)
            loc_coefMatrix=obj.coefMatrix[q_loc,q_loc]
        else
            # distMat=pairwise(Haversine(),obj.obsLonLat[q_loc,:])
            distMat=geoDist(obj.obsLonLat[q_loc,:])
            loc_coefMatrix=obj.locFun.(distMat/spaceCoef_cutoff);
        end
    end
    return q_loc,loc_coef,loc_coefMatrix
end
end # module