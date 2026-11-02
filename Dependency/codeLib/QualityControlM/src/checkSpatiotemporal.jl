
using Distributions
export checkSpatiotemporal
function checkSpatiotemporal(inData::AbstractArray{T},stalonlat;edgeLen_std_time=12*3,edgeLen_misfit_space=24*15,edgeLen_std_space=24*3,probabilityLimit=1e-7,analysis=false,disp=false,ind2ref_all=axes(inData,3),
        numRef = 10 # varing cutoff that ensure 10 sites to reference
        ) where T
    # dimension of inData should be (Time,varName,site)

    # prepare q2ref_all
    q2ref_all=falses(size(inData,3));q2ref_all[ind2ref_all] .= true;

    # check temporal 
    if disp;    print("temporal...");   end
    q_time,fit_time,std_time,pdf_time=mymapslices(x->
        checkTemporal(x;analysis=true,edgeLen_std=edgeLen_std_time)
        ,inData    );
    misfit_time= fit_time - inData
    # check space
    if disp;    println("spatial...");    end
    q_space=falses(size(inData))
    fit_space=fill(convert(T,NaN),size(inData))
    std_space=fill(convert(T,NaN),size(inData))
    pdf_space=fill(convert(T,NaN),size(inData))
    # Threads.@threads for n_var in axes(inData,2)
    for n_var in axes(inData,2)
        data_aVar=view(inData,:,n_var,:)
        q_validSta=vec(sum(.!isnan.(data_aVar),dims=1)) .> 0 # valid stations
        q_validRef=q_validSta .& q2ref_all
        localizer=localizerM(stalonlat[q_validSta,:],stalonlat[q_validRef,:])
        data_aVar_valid=view(inData,:,n_var,q_validSta)
        q_space[:,n_var,q_validSta],fit_space[:,n_var,q_validSta],std_space[:,n_var,q_validSta],pdf_space[:,n_var,q_validSta] =
        checkSpatial(data_aVar_valid,localizer;edgeLen_estimate=edgeLen_misfit_space,edgeLen_std=edgeLen_std_space,analysis=true,ind2ref_all=findall(q2ref_all[q_validSta]),
            numRef = numRef ) 
        # q_space[:,n_var,:],misfit_space[:,n_var,:],std_space[:,n_var,:],pdf_space[:,n_var,:] =
        # checkSpatial(view(inData,:,n_var,:),localizer;edgeLen_estimate=edgeLen_misfit_space,edgeLen_std=edgeLen_std_space,analysis=true,
        #     numRef = numRef ) 
    end
    misfit_space=fit_space - inData
    μ, Σ = [0, 0], [1 .4 ;.4 1]# the correlation is 0.4 by a sample
    pdf_st=similar(inData)
    Threads.@threads   for n_var in axes(pdf_st,2)
        normalizedUnionMisfit=Array{T}(undef,2,size(inData,1),size(inData,3))
        @views normalizedUnionMisfit[1,:,:] = -abs.(misfit_time[:,n_var,:])./std_time[:,n_var,:]
        @views normalizedUnionMisfit[2,:,:] = -abs.(misfit_space[:,n_var,:])./std_space[:,n_var,:]
        pdf_st[:,n_var,:]=pdf_2dMvNormal(μ, Σ,normalizedUnionMisfit)
    end
    q_st = pdf_st .< probabilityLimit

    # output
    if analysis 
        return q_st,fit_time,fit_space,std_time,std_space,pdf_st
    else
        return q_st
    end
end