export generateSpread
function generateSpread(beforeDA::AbstractArray{T};meanRatio=0.3,eleRatio=0.1,member=9) where T
    meanRatio=convert(T,meanRatio)
    eleRatio=convert(T,eleRatio)
    spread= ( meanRatio*mean(beforeDA) .+ eleRatio .* beforeDA ) .* 
                reshape(convert.(T,range(-1,1;length=member)),1,1,member)
    return spread
end