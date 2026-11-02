using GradDataM
export combinePollFromGrd,getSubVar

function combinePollFromGrd(grdIn,pollName;addSand=true,version="default")
    local pollData
    subVars,subRatio=getSubVar(pollName;addSand=addSand,version=version);
    for n_subVar in eachindex(subVars)
        subVarName=subVars[n_subVar];
        if n_subVar==1
            pollData=GradDataM.read(grdIn,subVarName).*subRatio[n_subVar];
        else
            pollData=pollData.+GradDataM.read(grdIn,subVarName)*subRatio[n_subVar];
        end
    end
    return pollData
end
function combinePollFromGrd(grdIn,pollName,nZ;addSand=true,version="default")
    local pollData
    subVars,subRatio=getSubVar(pollName;addSand=addSand,version=version);
    for n_subVar in eachindex(subVars)
        subVarName=subVars[n_subVar];
        # println(subVarName)
        if n_subVar==1
            pollData=read(grdIn,subVarName,nZ)*subRatio[n_subVar];
        else
            pollData=pollData+read(grdIn,subVarName,nZ)*subRatio[n_subVar];
        end
    end
    return pollData
end
function getSubVar(varName;addSand=true,version="default")
    if uppercase(varName) == "PM25"
        if version == "default"
            if addSand
                subVars=["pm25","aso4","anh4","ano3","bc","oc","soa1","soa2","soa3","soa4","soa5","soa6","aer0201","aer0202"];
            else
                subVars=["pm25","aso4","anh4","ano3","bc","oc","soa1","soa2","soa3","soa4","soa5","soa6"];
            end
            subRatio=ones(size(subVars));subRatio[6]=1.4;
        elseif version == "kong"
            subVars=["pm25","aso4","anh4","ano3","bc","oc","tsoa","fdust"];
            subRatio=ones(size(subVars));subRatio[6]=1.4;
        end
    elseif uppercase(varName) ==  "PM10"
        if version == "default"
            if addSand
                subVars=["pm25","pm10","aso4","anh4","ano3","bc","oc","soa1","soa2","soa3","soa4","soa5","soa6","aer0201","aer0202","aer0101","aer0102"];
            else
                subVars=["pm25","pm10","aso4","anh4","ano3","bc","oc","soa1","soa2","soa3","soa4","soa5","soa6"];
            end
            subRatio=ones(size(subVars));subRatio[7]=1.6;
        elseif version == "kong"
            subVars=["pm25","PM10","aso4","anh4","ano3","bc","oc","tsoa","fdust","cdust"];
            subRatio=ones(size(subVars));subRatio[7]=1.4;
        end
    elseif uppercase(varName) == "PMF"
        subVars=["PM25"];subRatio=ones(size(subVars));
    elseif uppercase(varName) == "PMC"
        subVars=["PM10"];subRatio=ones(size(subVars));
    elseif uppercase(varName) == "OLE"
        subVars=["OLET";"OLEI"];subRatio=ones(size(subVars));
    elseif uppercase(varName) == "SOA"
        subVars=["SOA1";"SOA2";"SOA3";"SOA4";"SOA5";"SOA6"];subRatio=ones(size(subVars));
    elseif uppercase(varName) == "AER"
        subVars=["aer0101";"aer0102";"aer0103";"aer0104";"aer0201";"aer0202";"aer0203";"aer0204"];subRatio=ones(size(subVars));
    elseif uppercase(varName) == "SAND"
        subVars=["pm25";"pm10";"aer0101";"aer0102";"aer0103";"aer0104";"aer0201";"aer0202";"aer0203";"aer0204"];subRatio=ones(size(subVars));
    else
        subVars=[varName];subRatio=1;
    end
    return subVars,subRatio
end