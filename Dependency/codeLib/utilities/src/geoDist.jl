using Distances
export geoDist
function geoDist(lonlatOri::AbstractVector,lonlatTarget::AbstractMatrix)
    return Haversine().([lonlatOri],eachrow(lonlatTarget))
end
function geoDist(lonlat::AbstractMatrix)
    return pairwise(Haversine(),eachrow(lonlat));
end