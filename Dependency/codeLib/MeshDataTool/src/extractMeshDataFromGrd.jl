using GradDataM,MeshDataM,NaqpmsUtilities,Dates,TimeZones
export extractMeshDataFromGrd
function extractMeshDataFromGrd(grdIn,datFileNameFun,meshNameOut,nZ_all,lead_all,var_all,tStart_all;varType="direct",transUnit=false)
    # datFileNameFun has two input arguments (tStart,tTarget)
    sizeLonLat=(grdIn.lonlatsize...,)
    dataAll=fill(NaN,(sizeLonLat...,length(nZ_all),length(lead_all),length(var_all),length(tStart_all)))
    for n_start in eachindex(tStart_all)
        tStart=tStart_all[n_start]
        Times=tStart .+ Hour.(lead_all)
        Threads.@threads for n_lead in eachindex(Times)
        # for n_lead=1:length(Times)
            tTarget=Times[n_lead]
            grd_curr=copy(grdIn)
            grd_curr.datFileName=datFileNameFun(tStart,tTarget)
            for n_var in eachindex(var_all)
                var_curr=var_all[n_var]
                # read from grd --------------------
                if varType == "direct"
                    data_curr=read(grd_curr,var_curr,nZ_all)
                elseif varType == "combine"
                    data_curr=combinePollFromGrd(grd_curr,var_curr,nZ_all)
                else
                    @error "varType not found: $varType"
                    data_curr=Array{Float32,2}(NaN,(sizeLonLat...,length(nZ_all)))
                end
                # unit ---------------
                if transUnit == true && (uppercase.(var_curr) in ["SO2","NO2","CO","O3"])
                    toMass!(data_curr,var_curr)
                end
                dataAll[:,:,:,n_lead,n_var,n_start]=data_curr
            end
        end
    end
    meshObj=meshDATA(meshNameOut,sizeLonLat,nZ_all,lead_all,var_all,tStart_all,dataAll)
    return meshObj
end