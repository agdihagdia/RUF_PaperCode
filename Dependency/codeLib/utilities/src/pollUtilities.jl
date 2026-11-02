
using Interpolations
export getAQIcali ,con2aqi,aqi2con,pollLevel,aqiRangeAccuracy,aqiLevelAccuracy,getAQIRange
    # get the AQI standard enacted by government, type name can be one of the
    # followings:  AQI SO2D SO2H NO2D NO2H PM10D COD COH O3H O3E PM25D
    
    # updates
    # 2016-3-22 add support for lower version of matlab by no using table,
    #           also, it may have a little speed up
    
        # Copyright: Huangjian Wu   All rights reserved.
        # If you have any problem, please contact me: wuhuangjian@126.com
    
        #AQIstandard=array2table(...
    #         [0,0,0,0,0,0,0,0,0,0,0;50,50,150,40,100,50,2,5,160,100,35;100,150,500,80,200,150,4,10,200,160,75;150,475,650,180,700,250,14,35,300,215,115;200,800,800,280,1200,350,24,60,400,265,150;300,1600,NaN,565,2340,420,36,90,800,800,250;400,2100,NaN,750,3090,500,48,120,1000,NaN,350;500,2620,NaN,940,3840,600,60,150,1200,NaN,500],...
    #         "VariableNames",...
    #         ["AQI","SO2D","SO2H","NO2D","NO2H","PM10D","COD","COH","O3H","O3E","PM25D"]);

function getAQIcali(typeName::String;standard="CHINA")
    if uppercase.(standard) == "CHINA"
        return getAQIcali_china(typeName)
    elseif uppercase.(standard) == "EPA"
        return getAQIcali_epa(typeName)
    else
        @error "valid standards are: CHINA EPA"
    end
end
function getAQIcali_china(typeName::String)
    typeName=uppercase.(typeName)
    if typeName in ["AQI"]
        cali=[0;50;100;150;200;300;400;500];#AQIstandard.AQI;
    elseif typeName in ["SO2D","SO2"]
        cali=[0;50;150;475;800;1600;2100;2620];#AQIstandard.SO2D;
    elseif typeName in ["SO2H"]
        cali=[0;150;500;650;800;1e6;1e7;1e8];#AQIstandard.SO2H;
    elseif typeName in ["NO2D","NO2"]
        cali=[0;40;80;180;280;565;750;940];#AQIstandard.NO2D;
    elseif typeName in ["NO2H"]
        cali=[0;100;200;700;1200;2340;3090;3840];#AQIstandard.NO2H;
    elseif typeName in ["COD","CO"]
        cali=[0;2;4;14;24;36;48;60];#AQIstandard.COD;
    elseif typeName in ["COH"]
        cali=[0;5;10;35;60;90;120;150];#AQIstandard.COH;
    elseif typeName in ["O3H","O3"]
        cali=[0;160;200;300;400;800;1000;1200];#AQIstandard.O3H;
    elseif typeName in ["O3E"]
        cali=[0;100;160;215;265;800;1e6;1e7];#AQIstandard.O3E;
    elseif typeName in ["PM25D","PM25"]
        cali=[0;35;75;115;150;250;350;500];#AQIstandard.PM25D;
    elseif typeName in ["PM10D","PM10"]
        cali=[0;50;150;250;350;420;500;600];#AQIstandard.PM10D;
    else
        @error "valid typeNames are: AQI SO2D SO2H NO2D NO2H PM10D COD COH O3H O3E PM25D"
    end
    return Float32.(cali)
end

function getAQIcali_epa(typeName::String)
    typeName=uppercase.(typeName)
    if typeName in ["AQI"]
        cali=[0;50;100;150;200;300;400;500];#AQIstandard.AQI;
    elseif typeName in ["PM25D","PM25"]
        cali=[0;12;35.5;55.5;150.5;250.5;350.5;500.5];#AQIstandard.PM25D;
    elseif typeName in ["PM10D","PM10"]
        cali=[0;55;155;255;355;425;505;605];#AQIstandard.PM10D;
    elseif typeName in ["SO2H","SO2"]
        cali=1000 .*toMass!([0;0.035;0.145;0.255;0.305;0.605;0.805;1.005],"so2");#AQIstandard.SO2D;
    elseif typeName in ["NO2D","NO2"]
        cali=1000 .*toMass!([0;0.054;0.101;0.361;0.65;1.25;1.65;2.049],"no2");#AQIstandard.NO2D;
    elseif typeName in ["COD","CO"]
        cali=1000 .*toMass!([0 ;4.5; 9.5; 12.5 ;15.5 ;30.5 ;40.5 ;50.5],"co");#AQIstandard.COD;
    elseif typeName in ["O3H","O3"]
        cali=1000 .*toMass!([0;0.125/2;0.125;0.165;0.205;0.405;0.505;0.605],"o3");#AQIstandard.O3H;
    elseif typeName in ["O3E"]
        cali=1000 .*toMass!([0;0.06;0.076;0.096;0.116;0.375;NaN;NaN],"o3");#AQIstandard.O3E;
    else
        @error "valid typeNames are: AQI SO2D SO2H NO2D NO2H PM10D COD COH O3H O3E PM25D"
    end
    return Float32.(cali)
end



    # convert concentration to api. support for multiple type (each column a type)
    # option.standard='epa';
        # dependent:
        # getAQIcali
        # typeName should be one of the case in getAQIcali
        # which is 'SO2D','SO2H','NO2D','NO2H','PM10D','COD','COH','O3H','O3E','PM25D'
        
        # sample usage:
        # aqi=getAQI('pm25d',[50 60])
        
        # Copyright: Huangjian Wu   All rights reserved.
        # If you have any problem, please contact me: wuhuangjian@126.com
function con2aqi(value::AbstractArray,typeName::String;standard="CHINA")
    
    aqiCali=getAQIcali("AQI";standard=standard)
    concenCali=getAQIcali(typeName;standard=standard)
    itp = extrapolate(interpolate((concenCali,), aqiCali, Gridded(Linear())),NaN)
    aqi=itp.(value)
    return aqi
end
function aqi2con(value::AbstractArray,typeName::String;standard="CHINA")
    aqiCali=getAQIcali("AQI";standard=standard)
    concenCali=getAQIcali(typeName;standard=standard)
    itp = extrapolate(interpolate((aqiCali,), concenCali, Gridded(Linear())),NaN)
    concen=itp.(value)
    return concen
end
function pollLevel(value,typeName::String;standard="CHINA")
    concenCali=getAQIcali(typeName;standard=standard)
    # f(x)=sum(x .> concenCali)
    # level=f.(value)
    level=sum(value .>= cPoint for cPoint in concenCali)
    return level
end
function aqiRangeAccuracy(obsCon,predCon,pollName)
    obsAQI=ceil.(con2aqi(obsCon,pollName))
    lowerAQI,upperAQI=getAQIRange(predCon,pollName)
    accuracy=mean((obsAQI .>= lowerAQI) .& (obsAQI .<= upperAQI))
    return accuracy
end
function aqiLevelAccuracy(obsCon,predCon,pollName;calLowUpRatio=false)
    # q_valid=.!(isnan.(obsCon) .| isnan.(predCon))
    q_valid=.!isnan.(obsCon) 
    obsCon=obsCon[q_valid]
    predCon=predCon[q_valid]
    obsLevel=pollLevel(obsCon,pollName)
    lowerAQI,upperAQI=getAQIRange(predCon,pollName)
    predUpperLevel=pollLevel(upperAQI,"AQI")
    predLowerLevel=pollLevel(lowerAQI,"AQI")
    accuracy=mean(predLowerLevel .<= obsLevel .<= predUpperLevel)
    if calLowUpRatio
        lowerRatio=mean(predUpperLevel .< obsLevel )
        upperRatio=mean(predLowerLevel .> obsLevel )
        return accuracy,lowerRatio,upperRatio
    else
        return accuracy
    end
end
function getAQIRange(inCon,pollName)
    inAQI=con2aqi(inCon,pollName)
    upperAQI=similar(inAQI)
    lowerAQI=similar(inAQI)
    q_good= inAQI .<= 50
    upperAQI[q_good] = inAQI[q_good] .* 1.1
    lowerAQI[q_good] = inAQI[q_good] .* 0.9
    upperAQI[.!q_good] = inAQI[.!q_good] .* 1.2
    lowerAQI[.!q_good] = inAQI[.!q_good] .* 0.8
    lowerAQI[lowerAQI .<0 ] .= 0
    return ceil.(lowerAQI),ceil.(upperAQI)
end