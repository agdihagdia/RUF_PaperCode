
function massBalanceFilter!(Flux_x,Flux_y,Flux_z,dx,dy,dz,denIn,denOut=denIn)
    # absGridFulx=abs.(Flux_x[1:end-1,:,:]) .+ abs.(Flux_x[2:end,:,:]) .+
    #             abs.(Flux_y[:,1:end-1,:]) .+ abs.(Flux_y[:,2:end,:]) .+
    #             abs.(Flux_z[:,:,1:end-1]) .+ abs.(Flux_z[:,:,2:end])
    # massChangeByFlux=-diff(Flux_x,dims=1) .- diff(Flux_y,dims=2) .- diff(Flux_z,dims=3)
    # massChangeTarget=(denOut-denIn) .* dx .* dy .* dz
    # mass2adjust=massChangeTarget .- massChangeByFlux
    for i in axes(denIn,1)
        for j in axes(denIn,2)
            for k in axes(denIn,3)
                massChangeTarget=(denOut[i,j,k]-denIn[i,j,k]) .* dx[i,j,k] .* dy[i,j,k] .* dz[i,j,k]
                massChangeByFlux = Flux_x[i,j,k] .- Flux_x[i+1,j,k] .+
                                   Flux_y[i,j,k] .- Flux_y[i,j+1,k] .+
                                   Flux_z[i,j,k] .- Flux_z[i,j,k+1]
                mass2adjust=massChangeTarget .- massChangeByFlux
                # adjust the flux to make sure the mass balance. only adjust the upper flux in i+1,j+1,k+1
                absFulx2change= abs(Flux_x[i+1,j,k]) + abs(Flux_y[i,j+1,k]) + abs(Flux_z[i,j,k+1])
                Flux_x[i+1,j,k]= Flux_x[i+1,j,k] + mass2adjust*abs(Flux_x[i+1,j,k])/absFulx2change
                Flux_y[i,j+1,k]= Flux_y[i,j+1,k] + mass2adjust*abs(Flux_y[i,j+1,k])/absFulx2change
                Flux_z[i,j,k+1]= Flux_z[i,j,k+1] + mass2adjust*abs(Flux_z[i,j,k+1])/absFulx2change
            end
        end
    end
    return Flux_x,Flux_y,Flux_z
end

function advec1d_byFlux( Q0::AbstractVector{T}, U, DEN0, DEN1, DXX, DD0, DT;
    IMXMN = zeros(T,size(Q0)), FLUX = zeros(T,size(Q0)), VCMAX = zeros(T,size(Q0)),VCMIN = zeros(T,size(Q0)),
    QN = zeros(T,size(Q0))  ) where T
    # Q0 is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
    # DXX is grid spacing, DD0 is the original density, DT is the time step
    sx = 1 + 1; ex = length(Q0) - 1;
    ZR0 = 0.0;
    LSHP = 0.0;

    # % identify local max and min, specigy mxing ratio limit at new time

    for I = sx:ex
        IMXMN[I] = 0.0;

        if (Q0[I] > max(Q0[I - 1], Q0[I + 1]) || Q0[I] < min(Q0[I - 1], Q0[I + 1]))
            IMXMN[I] = 1;
        end

        CK1 = Q0[I];
        if (U[I] < ZR0); CK1 = Q0[I + 1]; end
        if (U[I - 1] > ZR0); CK1 = Q0[I - 1]; end
        DEN1[I] = max(DEN1[I], 1.E-20);
        DXX[I] = max(DXX[I], 1.E-20);
        VCMAX[I] = max(Q0[I], CK1);
        VCMIN[I] = min(Q0[I], CK1);
    end

    # % update mixing ratios and limit Flux going up where u>0
    if (U[sx - 1] > ZR0); FLUX[sx - 1] = Q0[sx - 1] * U[sx - 1] * DT * DD0[sx - 1]; end
    for I = sx:ex
        if (U[I] < ZR0); continue; end

        if (U[I - 1] < ZR0)
            FLUX[I] = Q0[I] * U[I] * DT * DD0[I];# % outflow only cell
        else
            X1 = DT * U[I] / DXX[I];
            CF1 = CFCALC(X1, Q0[I], Q0[I - 1], Q0[I + 1], IMXMN[I - 1], IMXMN[I + 1], LSHP);
            FLUXI = DD0[I] * U[I] * DT * CF1 / DXX[I];
            QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUXI + FLUX[I - 1] / DXX[I]) / DEN1[I]));
            FLUX[I] = DXX[I] * (Q0[I] * DEN0[I] - QN[I] * DEN1[I]) + FLUX[I - 1];
        end
    end

    # % update mixing ratios and limit Flux going down where u<0
    if (U[ex] < ZR0); FLUX[ex] = Q0[ex + 1] * U[ex] * DT * DD0[ex]; end
    for I = ex:-1:sx
        if (U[I - 1] > ZR0)
            if (U[I] < ZR0)#% inflow-only cell
                QN[I] = (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUX[I - 1] / DXX[I]) / DEN1[I];
            end

        else
            X1 = DT * abs(U[I - 1]) / DXX[I];
            CF1 = CFCALC(X1, Q0[I], Q0[I + 1], Q0[I - 1], IMXMN[I + 1], IMXMN[I - 1], LSHP);

            if (U[I] > ZR0); CF1 = Q0[I]; end#% outflow only cell
            FLUXex1 = DD0[I - 1] * U[I - 1] * DT * CF1 / DXX[I];
            QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUXex1) / DEN1[I]));
            FLUX[I - 1] = DXX[I] * (QN[I] * DEN1[I] - Q0[I] * DEN0[I]) + FLUX[I];
        end

    end
    return QN
end
