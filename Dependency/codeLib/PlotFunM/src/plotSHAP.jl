using Random,Statistics
export plotBeesWarm
"""
xData=Float64.(reshape(data_shap.shap_effect,:,length(xVar)))
feature_value=Float64.(reshape(data_shap.feature_value,:,length(xVar)))
plt.close("all")
plt.figure(dpi=150)
plotBeesWarm(xData, xVar, feature_value;cMap="rainbow")
plt.gcf()
"""
function plotBeesWarm(xData, varName_all, colorData; spreadScale = 1,xBar_all=nothing,num_xSplice = 200,cMap="jet",addCBar=true)

    if isnothing(xBar_all)
        xBar_all = collect(range(minimum(xData) - 1e-20, stop=maximum(xData) + 1e-20, length=num_xSplice+1))
    end

    yData = getYData(xData, xBar_all;spreadScale=spreadScale)
    cData = (colorData .- mean(colorData;dims=1)) ./ std(colorData;dims=1)
    ax=plt.gca()
    hs=plt.scatter(xData,yData,c=cData,s=0.2,cmap=cMap)
    plt.xlabel("SHAP values (Impact on model output)")
    ax.set_yticks(1:length(varName_all))

    ax.set_yticklabels(varName_all)
    plt.clim([-1,1])
    if addCBar
    hbar = plt.colorbar()    
    hbar.set_ticks([-1, 1])
    hbar.set_ticklabels(["Low", "High"])
    else
        hbar=nothing
    end
    return hs,hbar
end

function getYData(xData, xBar_all; spreadScale=1)
    num_var = size(xData, 2)
    num_xSplice = length(xBar_all) - 1
    yData = fill(NaN, size(xData))
    count_splice = fill(NaN, num_xSplice, num_var)

    for n_splice in 1:num_xSplice
        count_splice[n_splice, :] = sum(xBar_all[n_splice] .<= xData .< xBar_all[n_splice + 1], dims=1)
    end

    spread_all = count_splice ./ maximum(count_splice, dims=1) * spreadScale
    spread_all[spread_all .> 1] .= 1

    for n_splice in 1:num_xSplice
        q_currSplice = xBar_all[n_splice] .<= xData .< xBar_all[n_splice + 1]
        for n_var in 1:num_var
            q_currVar = q_currSplice[:, n_var]
            num_currPoint = sum(q_currVar)
            if num_currPoint > 0
                if num_currPoint==1
                    yData[q_currVar, n_var] .= n_var 
                else
                    yData[q_currVar, n_var] = n_var .+ spread_all[n_splice, n_var] .* shuffle(range(-0.5,0.5,num_currPoint))
                 end
            end
        end
    end

    return yData
end

