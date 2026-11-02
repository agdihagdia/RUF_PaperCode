export boundMap
function boundMap(indata,lowerLimit,consistLevel,upperLimit)
    # % map unbounded indata to bounded outdata with limits specified and keep
    # % value and first derivative at the consistLevel
    
    q_lowerHalf=indata .< consistLevel;
    
    outdata=fill(NaN,size(indata));
    # % map the lower part
    outdata[q_lowerHalf]=consistLevel .+ 2 .* (lowerLimit .- consistLevel) ./ pi .*atan.((indata[q_lowerHalf] .- consistLevel) .*pi ./2 ./ (lowerLimit .- consistLevel)) ;
    # % map the upper part
    outdata[.~q_lowerHalf]=consistLevel .+ 2 .* (upperLimit .- consistLevel) ./ pi .*atan.((indata[.~q_lowerHalf] .- consistLevel) .* pi ./ 2 ./(upperLimit .-consistLevel)) ;
    
    return outdata
end
    