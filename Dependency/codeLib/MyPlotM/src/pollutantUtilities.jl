

import PyPlot as plt
export getUnitStr,varNameDisp,addPollYLabel # for pollutants 
# ========================= for pollutants ====================================
function varNameDisp(in)
    in=lowercase.(in)
    if in=="pm25"
        out=raw"PM$_{2.5}$"
    elseif in=="pm10"
        out=raw"PM$_{10}$"
    elseif in=="so2"
        out=raw"SO$_{2}$"
    elseif in=="no2"
        out=raw"NO$_{2}$"
    elseif in=="o3"
        out=raw"O$_{3}$"
    elseif in=="nh3"
        out=raw"NH$_{3}$"
    elseif in=="hno3"
        out=raw"HNO$_{3}$"
    elseif in=="ano3"
        out=raw"ANO$_{3}$"
    elseif in=="aso4"
        out=raw"ASO$_{4}$"
    elseif in=="co"
        out=raw"CO"
    else
        out = uppercase.(in)
    end
    return out
end

function getUnitStr(varName)
    varName=lowercase.(varName)
    if varName in ["pm25","pm10","so2","no2","o3"]
        unitStr=raw"μg m$^{-3}$"
    elseif varName=="co"
        unitStr=raw"mg m$^{-3}$"
    else
        unitStr =""
    end
    return unitStr
end

function addPollYLabel(varName)
    plt.ylabel( "$(varNameDisp(varName)) ($(getUnitStr(varName)))")
end
addPollYlabel=addPollYLabel
