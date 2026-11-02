using GradDataM,XLSX,DataFrames
export getEmitByID
function getEmitByID(grd_emit,idMat,dx,dy;var_all=grd_emit.varName_all,file_out=splitext(grd_emit.datFileName)[1] * ".xlsx")
    uniqueID_all=sort(unique(idMat[:]))
    var_curr=var_all[1]
    function getEmit(emitM,idM,dx,dy,id_curr)
        ratio= (60*60*24*365)/(1e6*1e6)# convert ug m-2 s-1 to t/year
        q_id=idM .== id_curr
        return ratio * sum(emitM[q_id] .* dx[q_id] .* dy[q_id])
    end
    statistics_curr=DataFrame(ID=uniqueID_all)
    for var_curr in var_all
        emit_curr=sum(read(grd_emit,var_curr),dims=3)[:,:,1]
        emitByID_curr=[ getEmit(emit_curr,idMat,dx,dy,id_curr) for id_curr in uniqueID_all]
        statistics_curr[!,var_curr]=emitByID_curr
    end
    statistics_curr=Float64.(statistics_curr[:,:])
    XLSX.writetable(file_out, statistics_curr, overwrite=true, sheetname="emitByID (t/year)")
    return statistics_curr
end