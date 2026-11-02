using NCDatasets,DataFrames,Dates,TimeZones,Interpolations
export calxco2,preprocess_OCO2,calxco2_3d
function calxco2(con_lev,press_lev)
    # Deficiencies: estimate the boundary pressure by linear interpolation, do not consider the humidity for dry density, revise latter
    diff_press=abs.(diff(press_lev))
    d_press=Array{Float32}(undef,length(press_lev))
    d_press[1]=diff_press[1];d_press[end]=diff_press[end]
    d_press[2:end-1]=(diff_press[1:end-1] + diff_press[2:end])./2
    ratio=d_press ./ sum(d_press)
    xco2=sum(con_lev .* ratio )
    return xco2
end
function calxco2_3d(con_lev,press_lev)
    # Deficiencies: estimate the boundary pressure by linear interpolation, do not consider the humidity for dry density, revise latter
    diff_press=abs.(diff(press_lev,dims=3))
    d_press=Array{Float32}(undef,size(press_lev))
    d_press[:,:,1]=diff_press[:,:,1];d_press[:,:,end]=diff_press[:,:,end]
    d_press[:,:,2:end-1]=(diff_press[:,:,1:end-1] .+ diff_press[:,:,2:end])./2
    ratio=d_press ./ sum(d_press;dims=3)
    xco2=sum(con_lev .* ratio ;dims=3)
    return xco2
end
# function calxco2(press_lev,con_lev)
#     # Deficiencies: estimate the boundary pressure by linear interpolation, do not consider the humidity for dry density, revise latter
#     itp_linear = LinearInterpolation([1:length(press_lev);], press_lev,extrapolation_bc=Line()) 
#     press_LevBound=itp_linear.([0.5:length(press_lev)+0.5;])
#     d_press=diff(press_LevBound)
#     ratio=d_press ./ sum(d_press)
#     xco2=sum(con_lev .* ratio )
#     return xco2
# end
# function calxco2(press_lev,con_lev,humi_lev)
#     function get_dryDense(humi)
#         const g=9.80665
#         const m_dry=28.9634
#         c=(1-humi)/(g*m_dry)
#         return c
#     end
#     itp_linear = LinearInterpolation([1:length(press_lev);], press_lev,extrapolation_bc=Line()) 
#     d_press=diff(itp_linear.([0.5:length(press_lev)+0.5;]))

#     dryDense_level=get_dryDense.(humi_lev)
#     xco2=sum(dryDense_level .* con_lev .* d_press)./sum(dryDense_level .* d_press)
#     return xco2
# end
function preprocess_OCO2(fileName,Times,uniqueLon,unqueLat)
    # read from oco2 file -----------------------
    ds = Dataset(fileName)
    df_obs_ori=DataFrame(time=ds["time"][:],xco2=ds["xco2"][:],lon=ds["longitude"][:],lat=ds["latitude"][:],uncertainty=ds["xco2_uncertainty"][:],flag=ds["xco2_quality_flag"][:])
    df_obs_ori.time=ZonedDateTime.(df_obs_ori.time,tz"UTC") 
    q_select = ( Times[1] .<= df_obs_ori.time .< Times[end]+Hour(1) ) .& ( 0 .==df_obs_ori.flag )
    df_obs_ori=df_obs_ori[q_select,:]

    # calculate unique tlonlat ---------------------------
    df_obs_out=similar(df_obs_ori,0)
    indNearestFun(xBase,x) =xBase[argmin(abs.(xBase.-x))]
    lon_nerest=(x->indNearestFun(uniqueLon,x)).(df_obs_ori.lon)
    lat_nerest=(x->indNearestFun(unqueLat,x)).(df_obs_ori.lat)
    tlonlat_in=zip(floor.(df_obs_ori.time,Hour),lon_nerest,lat_nerest)

    # prepare time lon lat --------------------------------
    tlonlat_out=unique(tlonlat_in)
    ind_inInOut=indexin(tlonlat_in,tlonlat_out)
    df_obs_out=similar(df_obs_ori,length(tlonlat_out))
    df_obs_out.time=[x[1] for x in tlonlat_out]
    df_obs_out.lon=[x[2] for x in tlonlat_out]
    df_obs_out.lat=[x[3] for x in tlonlat_out]

    # supper observation ----------------------
    df_obs_out.xco2=zeros(size(tlonlat_out,1))
    df_obs_out.uncertainty=zeros(size(tlonlat_out,1))
    numSample_inOut=zeros(size(tlonlat_out,1))
    for ind_in =1:length(tlonlat_in)
        ind_out=ind_inInOut[ind_in]
        # if 0==rem(ind_in,round(length(tlonlat_out)/100))
        #   println(ind_in,"  ",ind_out)
        # end
        df_obs_out.xco2[ind_out]+=df_obs_ori.xco2[ind_in]
        df_obs_out.uncertainty[ind_out]+=df_obs_ori.uncertainty[ind_in]
        numSample_inOut[ind_out]+=1
    end
    df_obs_out.xco2=df_obs_out.xco2./numSample_inOut
    df_obs_out.uncertainty=df_obs_out.uncertainty./numSample_inOut

    return df_obs_out
end
