
export ratio_toMass,toMass!,toVol!,toMass,toVol
function ratio_toMass(pollName)
    poll_lower=lowercase(pollName)
    if poll_lower in ["pm25", "pm10","aso4","ano3","bc","oc","anh4","tsoa","fdust","fdust","pno3","pso4","pnh4","sand"]
        ratio=1;
    elseif poll_lower ==  "so2"
        ratio=64/22.4;
    elseif poll_lower ==  "no"
        ratio=30/22.4;
    elseif poll_lower ==  "no2"
        ratio=46/22.4;
    elseif poll_lower ==  "no3"
        ratio=62/22.4;
    elseif poll_lower ==  "n2o5"
        ratio=108/22.4;
    elseif poll_lower ==  "o3"
        ratio=48/22.4;
    elseif poll_lower ==  "hono"
        ratio=47/22.4;
    elseif poll_lower ==  "hno3"
        ratio=63/22.4;
    elseif poll_lower ==  "nh3"
        ratio=17/22.4;
    elseif poll_lower ==  "co"
        ratio=28/22.4/1000;
    else
        ratio=1;
        @warn "the type of pollName : $pollName is not define yet"
    end
    return ratio
end
function toMass!(data,pollName)
    data .*=ratio_toMass(pollName)
    return data
end
function toMass(in,pollName)
    out= in.*ratio_toMass(pollName)
    return out
end
function toVol!(data,pollName)
    data ./=ratio_toMass(pollName)
    return data
end
function toVol(in,pollName)
    out =out./ratio_toMass(pollName)
    return out
end

