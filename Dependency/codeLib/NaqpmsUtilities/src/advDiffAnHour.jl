
# module TransportM
# export advDiffAnHour
using Polyester
export diffusion_hori_Smagorinsky,adv_hori_wuhj,adv_hori_basic,advDiffAnHour,adv_3d_basic,adv_3d_wuhj
export advection,advection_mass,adv_3d_onlineDen,advec1d,myAdvec1d
# functions for diffusion ******************************************

function boundaryMax(data)
    absLimit=prctile(abs.(data[:]),80)*1;
    q=abs.(data)>absLimit;
    data[q]=sign(data[q]).*absLimit;
    data[[ 1 end ],:]=0;data[:,[1 end]]=0;
    return data
end
function ∇²(f,dx,dy)
    ∇²f = zero(f) # initialize to zero since we don't touch the boundaries
    for y = 2:size(f,2)-1, x = 2:size(f,1)-1 # compute ∇²f in interior of f
        ∇²f[x,y] = (f[x-1,y] + f[x+1,y]- 4f[x,y])/dx[x,y]^2 + (f[x,y-1] + f[x,y+1] - 2f[x,y])/dy[x,y]^2
    end
    return ∇²f/4
 end
 function ∇(f,dx,dy)
     dfdx = zero(f) # initialize to zero since we don't touch the boundaries
     dfdy=zero(f)
     for y = 2:size(f,2)-1, x = 2:size(f,1)-1 # compute ∇²f in interior of f
         dfdx[x,y] =  (f[x+1,y]-f[x-1,y])/ (2*dx[x,y])
         dfdy[x,y] =  (f[x,y+1] -f[x,y-1])/(2*dy[x,y])
     end
     return dfdx,dfdy
  end
function diff_oneStep(in::AbstractArray{T},u,v,dx,dy,dt)where T
    # % method in statistic paper
    ddIn=∇²(in,dy,dx);
    # % add boundary limit
    # %ddIn=boundaryMax(ddIn);
    dudx,dudy=∇(u,dx,dy)
    dvdx,dvdy=∇(v,dx,dy)
    K=@. convert(T,(0.28)^2)*dx*dy*sqrt( (dudx- dvdy)^2 + (dudy+dvdx)^2 );
    out=in.+K.*ddIn.*dt;
    return out,ddIn
end

function diffusion_hori_Smagorinsky(in::AbstractArray{T},u,v,dx,dy,dt) where T
    minStepRatio=convert(T,0.2);
    out,ddIn=diff_oneStep(in,u,v,dx,dy,dt);
    if 0==mean(abs.(ddIn.*dx.*dy)); return in;end# not diffusion return input; advoid 0 in the denominator
    stepRatio=mean(abs.(out[:]-in[:]))/mean(abs.(ddIn.*dx.*dy));#% step change / max change in neighbour grid.
    if stepRatio<minStepRatio#% done
        return out;
    else
        t_stepS=ceil(stepRatio/minStepRatio);
        previous=in;
        for n_t=1:t_stepS
            out,~=diff_oneStep(previous,u,v,dx,dy,dt/t_stepS);
            previous=out;
        end
    end
    return out
end
struct advecStr
    IMXMN
    FLUX
    VCMAX
    VCMIN
    QN
end
# struct adv
function reset!(advecVar::advecStr)
    advecVar.IMXMN .= 0
    advecVar.FLUX .= 0
    advecVar.VCMAX .= 0
    advecVar.VCMIN .= 0
    advecVar.QN .= 0
    return advecVar
end
function init_advec(Q0::AbstractVector{T})  where T
    IMXMN,FLUX,VCMAX,VCMIN,QN=zeros(T,size(Q0)),zeros(T,size(Q0)),zeros(T,size(Q0)),zeros(T,size(Q0)),zeros(T,size(Q0))
    return advecStr(IMXMN,FLUX,VCMAX,VCMIN,QN)
end

function advec1d( Q0::AbstractVector{T}, U, DEN0, DEN1, DXX, DD0, DT;
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

        if (U[I - 1] < ZR0) # <- * ->
            FLUX[I] = Q0[I] * U[I] * DT * DD0[I];# % outflow only cell
        else # -> * ->
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
            if (U[I] < ZR0) # -> * <-
                QN[I] = (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUX[I - 1] / DXX[I]) / DEN1[I];
                # println(I," ",QN[I]," ",Q0[I]," ",FLUX[I]," ",FLUX[I-1]," ",DEN1[I])
            end
        else # <- * ==
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

function advec1d_z( Q0::AbstractVector{T}, U, DEN0, DEN1, DXX, DD0, DT;
    IMXMN = zeros(T,size(Q0)), FLUX = zeros(T,size(Q0)), VCMAX = zeros(T,size(Q0)),VCMIN = zeros(T,size(Q0)),
    QN = zeros(T,size(Q0))  ) where T
    # Q0 is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
    # DXX is grid spacing, DD0 is the original density, DT is the time step
    sx = 1; ex = length(Q0) ;
    ZR0 = 0.0;
    LSHP = 0.0;

    # % identify local max and min, specigy mxing ratio limit at new time
    for I = sx:ex
        IMXMN[I] = 0.0;
        if 1==I
            IMXMN[I] = 1;
            CK1 = Q0[I];
        elseif I==length(Q0)
            IMXMN[I] = 1;
            CK1 = Q0[I];
            if (U[I - 1] > ZR0); CK1 = Q0[I - 1]; end
        else
            if (Q0[I] > max(Q0[I - 1], Q0[I + 1]) || Q0[I] < min(Q0[I - 1], Q0[I + 1]))
                IMXMN[I] = 1;
            end
            CK1 = Q0[I];
            if (U[I] < ZR0); CK1 = Q0[I + 1]; end
            if (U[I - 1] > ZR0); CK1 = Q0[I - 1]; end
        end
        DEN1[I] = max(DEN1[I], 1.E-20);
        DXX[I] = max(DXX[I], 1.E-20);
        VCMAX[I] = max(Q0[I], CK1);
        VCMIN[I] = min(Q0[I], CK1);
    end

    # % update mixing ratios and limit Flux going up where u>0
    if (U[sx] > ZR0); FLUX[sx] = Q0[sx] * U[sx ] * DT * DD0[sx]; end
    for I = sx:ex
        if (U[I] < ZR0); continue; end

        if I==1 || (U[I - 1] < ZR0)
            FLUX[I] = Q0[I] * U[I] * DT * DD0[I];# % outflow only cell
        else

            X1 = DT * U[I] / DXX[I];
            if I == length(Q0)
                CF1 = CFCALC(X1, Q0[I], Q0[I - 1], Q0[I], IMXMN[I - 1], IMXMN[I], LSHP);
            else
                CF1 = CFCALC(X1, Q0[I], Q0[I - 1], Q0[I + 1], IMXMN[I - 1], IMXMN[I + 1], LSHP);
            end
            FLUXI = DD0[I] * U[I] * DT * CF1 / DXX[I];
            QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUXI + FLUX[I - 1] / DXX[I]) / DEN1[I]));
            FLUX[I] = DXX[I] * (Q0[I] * DEN0[I] - QN[I] * DEN1[I]) + FLUX[I - 1];
        end
    end

    # % update mixing ratios and limit Flux going down where u<0
    if (U[ex] < ZR0); FLUX[ex] = Q0[ex] * U[ex] * DT * DD0[ex]; end
    for I = ex:-1:sx
        if I==1 || (U[I - 1] > ZR0)
            if (U[I] < ZR0)#% inflow-only cell
                if I ==1
                    QN[I] = (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUX[I] / DXX[I]) / DEN1[I];
                else
                    QN[I] = (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUX[I - 1] / DXX[I]) / DEN1[I];
                end
            end
        else
            if I==length(Q0)
                X1 = DT * abs(U[I - 1]) / DXX[I];
                CF1 = CFCALC(X1, Q0[I], Q0[I], Q0[I - 1], IMXMN[I], IMXMN[I - 1], LSHP);
    
                if (U[I] > ZR0); CF1 = Q0[I]; end#% outflow only cell
                FLUXex1 = DD0[I - 1] * U[I - 1] * DT * CF1 / DXX[I];
                QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUXex1) / DEN1[I]));
                FLUX[I - 1] = DXX[I] * (QN[I] * DEN1[I] - Q0[I] * DEN0[I]) + FLUX[I];
            elseif I == 1 
                X1 = DT * abs(U[I]) / DXX[I];
                CF1 = CFCALC(X1, Q0[I], Q0[I + 1], Q0[I], IMXMN[I + 1], IMXMN[I], LSHP);
                if (U[I] > ZR0); CF1 = Q0[I]; end#% outflow only cell
                FLUXex1 = DD0[I] * U[I] * DT * CF1 / DXX[I];
                QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUXex1) / DEN1[I]));
            else
                X1 = DT * abs(U[I - 1]) / DXX[I];
                CF1 = CFCALC(X1, Q0[I], Q0[I + 1], Q0[I - 1], IMXMN[I + 1], IMXMN[I - 1], LSHP);

                if (U[I] > ZR0); CF1 = Q0[I]; end#% outflow only cell
                FLUXex1 = DD0[I - 1] * U[I - 1] * DT * CF1 / DXX[I];
                QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUX[I] / DXX[I] + FLUXex1) / DEN1[I]));
                FLUX[I - 1] = DXX[I] * (QN[I] * DEN1[I] - Q0[I] * DEN0[I]) + FLUX[I];
            end
        end

    end
    return QN
end

"""
function myAdvec1d( Qin::AbstractVector{T}, U_pre,U_next, DEN0, DEN1, ∆x, DD0, DT;
    q_extrame = falses(size(Qin)), Flux_all = zeros(T,size(Qin)+1), Qmax = similar(Qin),Qmin = similar(Qin), Qout = similar(Qin), 
    fullSharpening = false) where T
    # Qin is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
    # DXX is grid spacing, DD0 is the original density, DT is the time step

    Flux_pre=view(Flux_all,1:length(Flux_all)-1)
    Flux_next=view(Flux_all,2:length(Flux_all))


    # % identify local max and min, specigy mxing ratio limit at new time
    Qmin .= Qin;Qmax .= Qin;
    for i in eachindex(Qin)   
        # set Q_left and Q_right
        if i==1
            Q_left,Q_right=Qin[i],Qin[i+1]
        elseif i==length(Qin)
            Q_left,Q_right=Qin[i-1],Qin[i]
        else
            Q_left,Q_right=Qin[i-1],Qin[i+1]
        end
        # set q_extrame
        if  Qin[i] > max(Q_left, Q_right) || Qin[i] < min(Q_left, Q_right)
            q_extrame[i]=true
        end    
        # set Qmax and Qmin
        if U_pre[i] > 0;        Qmin[i],Qmax[i] = min(Q_left,Qmin[i]),max(Q_left,Qmax[i]);    end #isUpwind_left
        if U_next[i] < 0;       Qmin[i],Qmax[i] = min(Q_right,Qmin[i]),max(Q_right,Qmax[i]);  end #isUpwind_right
    end

    #  update mixing ratios and  Flux for u>0
    # the updated flux will be used in the next section
    Flux_pre[1] = Qin[1] * U_pre[1] * DT * DD0[1];
    for i in eachindex(Qin)
        if U_next[i] > 0 # === * ->
            if U_pre[i] < 0 # <- * ->
                Flux_next[i] = Qin[i] * U_next[i] * DT * DD0[i];# % outflow only cell
            else # -> * ->
                X1 = DT * U_next[i] / ∆x[i]; # raio of transport length to grid spacing
                if i==1;                Q_upwind=Qin[1];Q_downwind=Qin[2];isExtrame_upwind=false;isExtrame_downwind=q_extrame[2]; # handel left boundary
                elseif i==length(Qin);   Q_upwind=Qin[end-1];Q_downwind=Qin[end];isExtrame_upwind=q_extrame[end-1];isExtrame_downwind=false;# handel right boundary
                else Q_upwind=Qin[i-1];     Q_downwind=Qin[i+1];isExtrame_upwind=q_extrame[i-1];isExtrame_downwind=q_extrame[i+1];
                end
                Q_crossFace_downwind = myCFCALC(X1, Qin[i], Q_upwind, Q_downwind, isExtrame_upwind, isExtrame_downwind, fullSharpening);
                Flux_next_tmp = DD0[i] * U_next[i] * DT * Q_crossFace_downwind / ∆x[i];
                Qout[i] = max(Qmin[i], min(Qmax[i], (Qin[i] * DEN0[i] - Flux_next_tmp + Flux_pre[i] / ∆x[i]) / DEN1[i]));
                Flux_next[i] = ∆x[i] * (Qin[i] * DEN0[i] - Qout[i] * DEN1[i]) + Flux_pre[i];
            end
        end
    end

    #  update mixing ratios and Flux for u < 0
    Flux_next[end] = Qin[end] * U_next[end] * DT * DD0[end];
    for i in reverse(eachindex(Qin))
        if U_pre[i] > 0 # -> * ==
            if U_next[i] <= 0# -> * <- inflow only cell
                Qout[i] = (Qin[i] * DEN0[i] - Flux_next[i] / ∆x[i] + Flux_pre[i] / ∆x[i]) / DEN1[i];
                # Qout[i] = Qin[i] 
            end
        else # <- * ==
            X1 = DT * abs(U_pre[i]) / ∆x[i];
            if U_next[i] > 0; # <- * ->  outflow only cell
                Q_crossFace_downwind = Qin[i]; 
            else # <- * <- 
                if i==1;                Q_upwind=Qin[2];Q_downwind=Qin[1];isExtrame_upwind=q_extrame[2];isExtrame_downwind=false; # handel left boundary
                elseif i==length(Qin);   Q_upwind=Qin[end];Q_downwind=Qin[end-1];isExtrame_upwind=false;isExtrame_downwind=q_extrame[end-1];# handel right boundary
                else Q_upwind=Qin[i+1];     Q_downwind=Qin[i-1];isExtrame_upwind=q_extrame[i+1];isExtrame_downwind=q_extrame[i-1];
                end
                Q_crossFace_downwind = myCFCALC(X1, Qin[i], Q_upwind, Q_downwind, isExtrame_upwind, isExtrame_downwind, fullSharpening);
            end
            if i ==1; DD0_curr=DD0[1];else DD0_curr=DD0[i-1];end
            Flux_pre_tmp = DD0_curr * U_pre[i] * DT * Q_crossFace_downwind / ∆x[i];
            Qout[i] = max(Qmin[i], min(Qmax[i], (Qin[i] * DEN0[i] -  Flux_next[i] / ∆x[i] + Flux_pre_tmp) / DEN1[i]));
            Flux_pre[i] = ∆x[i] * (Qout[i] * DEN1[i] - Qin[i] * DEN0[i]) + Flux_next[i];
        end
    end
    return Qout
end
"""

# 计算下风界面的混合比Qf
# CFCALC Cross Face CALCulation
# CF1 = CFCALC(X1, Q0[I], Q0[I + 1], Q0[I - 1], q_locMaxOrMin[I + 1], q_locMaxOrMin[I - 1], LSHP);
function myCFCALC(X1, Q_curr::T, Q_upwind, Q_downwind, isExtrema_upwind, isExtrema_downwind, fullSharpening=false) where T
    # X1 is the ratio of transport length to grid spacing
    # % function to calculate mixing ratio in the fuild moving into neighnor
    # % cell durint DT; X1=Courant No. set LSHP=1 for full sharpening
    local ALFA::T, CF::T
    if fullSharpening # sharpening factor
        ALFA=5.0; 
    else        
        if isExtrema_upwind ; ALFA = 1.2578 + 0.58502 * X1; # 上风向是极值点
        elseif isExtrema_downwind ; ALFA = 1.7763 - 0.5 * X1;  # 下风向是极值点
        else ALFA=1.0;
        end 
    end
    if X1 < 0.5 # 传输不足半个网格
        ALFA = min(ALFA, .98 * 6 ./ (1 - 2 .* X1));
        CF = Q_curr * (1 + ALFA / 6 .* (2 .* X1 - 1.)) + Q_downwind * ALFA / 12 .* (4 - 5 .* X1) + Q_upwind * ALFA / 12 .* (X1 - 2.);
    else # 传输超过半个网格
        CF = Q_curr * (1 - ALFA * (1 ./ X1 + 2 .* X1 -3.) / 6.) + Q_downwind * ALFA * (1 ./ X1 - X1) / 12 + Q_upwind * ALFA * (1 ./ X1 + 5 .* X1 - 6.) / 12.;
    end

    # % limit outflow mixing ratio to reasonable mixing ratio
    out = min(max(CF, min(Q_curr, Q_downwind)), max(Q_curr, Q_downwind));
    return out
end

# 计算界面混合比Qf？
# CF1 = CFCALC(X1, Q0[I], Q0[I + 1], Q0[I - 1], q_locMaxOrMin[I + 1], q_locMaxOrMin[I - 1], LSHP);
function CFCALC(X1, VCUP1, VCUP2, VCDW1, IMXUP2, IMXDN1, LSHP)
    # % function to calculate mixing ratio in the fuild moving into neighnor
    # % cell durint DT; X1=Courant No. set LSHP=1 for full sharpening
    
    if LSHP == 1
        ALFA = 5;
    else
        ALFA = 1;
        if IMXDN1 > 0; ALFA = 1.7763 - 0.5 * X1; end # 上一个网格是极值点
        if IMXUP2 > 0; ALFA = 1.2578 + 0.58502 * X1; end # 下一个网格是极值点
    end
    
    if X1 < 0.5
        ALFA = min(ALFA, .98 * 6 ./ (1 - 2 .* X1));
        CF = VCUP1 * (1 + ALFA / 6 .* (2 .* X1 - 1.)) + VCDW1 * ALFA / 12 .* (4 - 5 .* X1) + VCUP2 * ALFA / 12 .* (X1 - 2.);
    else
        CF = VCUP1 * (1 - ALFA * (1 ./ X1 + 2 .* X1 -3.) / 6.) + VCDW1 * ALFA * (1 ./ X1 - X1) / 12 + VCUP2 * ALFA * (1 ./ X1 + 5 .* X1 - 6.) / 12.;
    end

    # % limit outflow mixing ratio to reasonable mixing ratio
    out = min(max(CF, min(VCUP1, VCDW1)), max(VCUP1, VCDW1));
    return out
end

function adv_hori_wuhj(c_in::AbstractArray{T}, u, v, ∆x, ∆y, dt;par=false,num_tStep=1) where T
    # % the boundary will not change, however, it is necessary for calculating the inner area
    # %function adv_hori( myid, c_out, u, v, ∆x, ∆y, sx, ex, sy, ey ,nx,ny,dt)
    c_out=copy(c_in);

    DEN1 = zeros(T,size(c_out)); DEN1y = zeros(T,size(c_out));

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_out)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_out)); end

    d0 = ones(T,size(c_out));

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = minimum(abs.(vcat(vec(∆x ./ u), vec(∆y ./ v))));
        num_tStep = ceil(dt / dt0);
    end
    dtt = dt / num_tStep;
    # if num_tStep> 1;println("num_tStep= ",num_tStep);end
    if !(num_tStep <50)
        @warn "too much num_tStep!($num_tStep) dtt=$dtt"
    end

    if par
        advecVar_all_dim1=[init_advec(c_in[:,1]) for i in 1:Threads.nthreads()]
        advecVar_all_dim2=[init_advec(c_in[1,:]) for i in 1:Threads.nthreads()]
    else
        advecVar_dim1=init_advec(c_in[:,1])
        advecVar_dim2=init_advec(c_in[1,:])
    end

    for itt = 1:num_tStep
        # % for i direction
        Q0 = c_out;
        u1D = u;
        DXX = ∆x;
        DEN0 = d0;
        @views DEN1[2:end, :] = 1 .- dtt ./ DXX[2:end, :] .* (d0[2:end, :] .* u[2:end, :] .- d0[2:end, :] .* u[1:end - 1, :]);
        DEN1[1, :] .= 0;
        DD0 = DEN0;
        QN=c_out;
        if par == true
            Threads.@threads for ny = 2:(size(QN,2)-1)
                advecVar_curr=reset!(advecVar_all_dim1[Threads.threadid()])
                @views QN[:,ny] =advec1d( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt;                
                    IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
            end
        else
            for ny = 2:(size(QN,2)-1)
                advecVar_curr=advecVar_dim1
                @views QN[:,ny] =advec1d( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt;
                    IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
            end
        end
        c_out = QN;
        # % for j direction
        Q0y = c_out;
        u1Dy = v;
        DXXy = ∆y;
        DEN0y = DEN1;
        DD0y = d0;
        @views DEN1y[:, 1] = DEN0y[:, 1];
        @views  DEN1y[:, 2:end] = DEN0y[:, 2:end] .- dtt ./ ∆y[:, 2:end] .* (d0[:, 2:end] .* v[:, 2:end] .- d0[:, 2:end] .* v[:, 1:end - 1]);
        QNy=c_out;
        if par == true
            Threads.@threads  for nx = 2:(size(QNy,1)-1)
                advecVar_curr=reset!(advecVar_all_dim2[Threads.threadid()])
                @views QNy[nx,:] =advec1d( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt;
                IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
            end
        else
            for nx = 2:(size(QNy,1)-1)
                advecVar_curr=advecVar_dim2
                @views QNy[nx,:] =advec1d( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt;
                    IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
            end
        end
        c_out = QNy;
        # % preserve the boundary
        @views c_out[[1 end],:]=c_in[[1 end],:];
        @views c_out[:,[1 end]]=c_in[:,[1 end]];
        # %c_out=QNy.*DEN1y;
    end
    return c_out
end
function adv_hori_basic(c_in::AbstractArray{T}, u, v, ∆x, ∆y, dt;par=false,num_tStep=1) where T
    # % the boundary will not change, however, it is necessary for calculating the inner area
    # %function adv_hori( myid, c_out, u, v, ∆x, ∆y, sx, ex, sy, ey ,nx,ny,dt)
    c_out=copy(c_in);

    DEN1 = zeros(T,size(c_out)); DEN1y = zeros(T,size(c_out));

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_out)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_out)); end

    d0 = ones(T,size(c_out));

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        # dt0 = minimum(abs.([∆x[:] ./ u[:]; ∆y[:] ./ v[:]]));
        dt0 = minimum(abs.(vcat(vec(∆x ./ u), vec(∆y ./ v))));
        num_tStep = ceil(dt / dt0);
    end
    dtt = dt / num_tStep;

    advecVar_dim1=init_advec(c_in[:,1])
    advecVar_dim2=init_advec(c_in[1,:])
    for itt = 1:num_tStep
        # % for i direction
        Q0 = c_out;
        u1D = u;
        DXX = ∆x;
        DEN0 = d0;

        @views DEN1[2:end, :] = 1 .- dtt ./ DXX[2:end, :] .* (d0[2:end, :] .* u[2:end, :] .- d0[2:end, :] .* u[1:end - 1, :]);
        DEN1[1, :] .= 0;
        DD0 = DEN0;
        QN=c_out;
        
        for ny = 2:(size(QN,2)-1)
            @views QN[:,ny] =advec1d( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt;            
                IMXMN = advecVar_dim1.IMXMN, FLUX = advecVar_dim1.FLUX, VCMAX = advecVar_dim1.VCMAX,VCMIN = advecVar_dim1.VCMIN,   QN = advecVar_dim1.QN);
        end
        
        c_out = QN;
        # % for j direction
        Q0y = c_out;
        u1Dy = v;
        DXXy = ∆y;
        DEN0y = DEN1;
        DD0y = d0;
        @views DEN1y[:, 1] = DEN0y[:, 1];
        @views DEN1y[:, 2:end] = DEN0y[:, 2:end] .- dtt ./ ∆y[:, 2:end] .* (d0[:, 2:end] .* v[:, 2:end] .- d0[:, 2:end] .* v[:, 1:end - 1]);
        QNy=c_out;
        for nx = 2:(size(QNy,1)-1)
            @views QNy[nx,:] =advec1d( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt;
                    IMXMN = advecVar_dim2.IMXMN, FLUX = advecVar_dim2.FLUX, VCMAX = advecVar_dim2.VCMAX,VCMIN = advecVar_dim2.VCMIN,   QN = advecVar_dim2.QN);
        end
        c_out = QNy;
        # % preserve the boundary
        @views c_out[[1 end],:] .= c_in[[1 end],:];
        @views c_out[:,[1 end]] .= c_in[:,[1 end]];
        # %c_out=QNy.*DEN1y;
    end
    return c_out
end
"""
ρ = P/(R×T)

"""
function adv_3d_onlineDen(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=-1,massBalance=false) where T
    # input C_in is the mass concentration
    # set num_tStep=0 to calculate the num_tStep automatically by wind speed and grid size
    # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
    # denOut=copy(denIn);
    Q=similar(c_in);
    c_out=copy(c_in)
    if size(u)==size(c_out); u=cat(u[[1],:,:],u;dims=1); end
    if size(v)==size(c_out); v=cat(v[:,[1],:],v;dims=2); end
    if size(w)==size(c_out); w=cat(w,w[:,:,[1]];dims=3); end
    @views u_pre,u_next=u[1:end-1,:,:],u[2:end,:,:]
    @views v_pre,v_next=v[:,1:end-1,:],v[:,2:end,:]
    @views w_pre,w_next=w[:,:,1:end-1],w[:,:,2:end]

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_in)); end
    if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = min(minimum(abs,∆x ./ u_next), minimum(abs,∆y ./ v_next),minimum(abs,∆z ./ w_pre));
        num_tStep = ceil(dt / dt0);
    end
    dtt = convert(T,dt / num_tStep);
    if !(num_tStep <50)
        @warn "too much num_tStep!($num_tStep) dtt=$dtt"
    end

    # allocation for the advection -----------------------
    advStr_x=init_myAdvec(c_in[:,1,1])
    advStr_y=init_myAdvec(c_in[1,:,1])
    advStr_z=init_myAdvec(c_in[1,1,:])
  
    den_minimum=minimum(denIn)/1e3
    # do the advection -----------------------
    denOut=copy(denIn)
    D0=denOut; D1=similar(denOut);D2=similar(denOut);D3=similar(denOut);
    if !massBalance; Q .= c_in; end
    for itt = 1:num_tStep
        # deal with densities
        if massBalance; Q .= c_out./denOut; end

        # update densities -------------------------------
        D0 .= denOut
        # for i direction  
        for ny in axes(D1,2)
            massFlux_x=Array{T,1}(undef,size(D1,1)+1)
            for nz in axes(D1,3)
                @views D1[:,ny,nz]=updateDen1d!(D0[:,ny,nz],D1[:,ny,nz],D0[:,ny,nz],u_pre[:,ny,nz],u_next[:,ny,nz],∆x[:,ny,nz],dtt,den_minimum; massFlux=massFlux_x)
            end
        end
        # for j direction
        for nx in axes(D2,1)
            massFlux_y=Array{T,1}(undef,size(D2,2)+1)
            for nz in axes(D2,3)
                @views D2[nx,:,nz]=updateDen1d!(D1[nx,:,nz],D2[nx,:,nz],D0[nx,:,nz],v_pre[nx,:,nz],v_next[nx,:,nz],∆y[nx,:,nz],dtt,den_minimum; massFlux=massFlux_y)
            end
        end
        # for k direction
        for nx in axes(D3,1)
            massFulx_z=Array{T,1}(undef,size(D3,3)+1)
            for ny in axes(D3,2)
                @views D3[nx,ny,:]=updateDen1d!(D2[nx,ny,:],D3[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt,den_minimum; massFlux=massFulx_z)
            end
        end
        denOut .= D3
        # update mixing ratio -------------------------------
        # for i direction
        for ny in axes(Q,2), nz in axes(Q,3)
            @views Q[:,ny,nz] =myAdvec1d( Q[:,ny,nz] , u_pre[:,ny,nz] , u_next[:,ny,nz] , D0[:,ny,nz] , D1[:,ny,nz] , ∆x[:,ny,nz] , D0[:,ny,nz] , dtt;   
                q_extrame = advStr_x.q_extrame, Flux_all = advStr_x.Flux_all, Qmax = advStr_x.Qmax,Qmin = advStr_x.Qmin,   Qout = advStr_x.Qout);
        end
        # for j direction
        for nx in axes(Q,1), nz in axes(Q,3)
            @views Q[nx,:,nz] =myAdvec1d( Q[nx,:,nz] , v_pre[nx,:,nz] , v_next[nx,:,nz] , D1[nx,:,nz] , D2[nx,:,nz] , ∆y[nx,:,nz] , D0[nx,:,nz] , dtt;
                    q_extrame = advStr_y.q_extrame, Flux_all = advStr_y.Flux_all, Qmax = advStr_y.Qmax,Qmin = advStr_y.Qmin,   Qout = advStr_y.Qout);
        end
        # for k direction
        for nx in axes(Q,1), ny in axes(Q,2)
            @views Q[nx,ny,:] =myAdvec1d( Q[nx,ny,:] , w_pre[nx,ny,:] , w_next[nx,ny,:] , D2[nx,ny,:] , D3[nx,ny,:] , ∆z[nx,ny,:] , D0[nx,ny,:] , dtt;
                q_extrame = advStr_z.q_extrame, Flux_all = advStr_z.Flux_all, Qmax = advStr_z.Qmax,Qmin = advStr_z.Qmin,   Qout = advStr_z.Qout);
        end
        if massBalance; c_out .= Q .* denOut; end
        # c_out .= Q.*denOut;
    end
    # c_out = Q.*denIn;
    if !massBalance; c_out .= Q;    end
    return c_out
end

function adv_3d_basic(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=1) where T
    # input C_in is the mass concentration
    # set num_tStep=0 to calculate the num_tStep automatically by wind speed and grid size
    # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
    # denOut=copy(denIn);
    Q=similar(c_in);
    c_out=copy(c_in)
    if size(u)==size(c_out); u=cat(u[[1],:,:],u;dims=1); end
    if size(v)==size(c_out); v=cat(v[:,[1],:],v;dims=2); end
    if size(w)==size(c_out); w=cat(w,w[:,:,[1]];dims=3); end
    @views u_pre,u_next=u[1:end-1,:,:],u[2:end,:,:]
    @views v_pre,v_next=v[:,1:end-1,:],v[:,2:end,:]
    @views w_pre,w_next=w[:,:,1:end-1],w[:,:,2:end]

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_in)); end
    if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = min(minimum(abs,∆x ./ u_next), minimum(abs,∆y ./ v_next),minimum(abs,∆z ./ w_pre));
        num_tStep = ceil(dt / dt0);
    end
    dtt = convert(T,dt / num_tStep);

    # allocation for the advection -----------------------
    advStr_x=init_myAdvec(c_in[:,1,1])
    advStr_y=init_myAdvec(c_in[1,:,1])
    advStr_z=init_myAdvec(c_in[1,1,:])

    # deal for densities -----------------------    
    den_minimum=minimum(denIn)/1e3
    D0=denIn; D1=similar(denIn);D2=similar(denIn);D3=similar(denIn);
    # for i direction  
    for ny in axes(D1,2)
        massFlux_x=Array{T,1}(undef,size(D1,1)+1)
        for nz in axes(D1,3)
            @views D1[:,ny,nz]=updateDen1d!(D0[:,ny,nz],D1[:,ny,nz],D0[:,ny,nz],u_pre[:,ny,nz],u_next[:,ny,nz],∆x[:,ny,nz],dtt,den_minimum; massFlux=massFlux_x)
        end
    end
    # for j direction
    for nx in axes(D2,1)
        massFlux_y=Array{T,1}(undef,size(D2,2)+1)
        for nz in axes(D2,3)
            @views D2[nx,:,nz]=updateDen1d!(D1[nx,:,nz],D2[nx,:,nz],D0[nx,:,nz],v_pre[nx,:,nz],v_next[nx,:,nz],∆y[nx,:,nz],dtt,den_minimum; massFlux=massFlux_y)
        end
    end
    # for k direction
    for nx in axes(D3,1)
        massFulx_z=Array{T,1}(undef,size(D3,3)+1)
        for ny in axes(D3,2)
            @views D3[nx,ny,:]=updateDen1d!(D2[nx,ny,:],D3[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt,den_minimum; massFlux=massFulx_z)
        end
    end
    denOut = D3
    # do the advection -----------------------
    for itt = 1:num_tStep
        Q .= c_out./denIn
        # for i direction
        for ny in axes(Q,2), nz in axes(Q,3)
            @views Q[:,ny,nz] =myAdvec1d( Q[:,ny,nz] , u_pre[:,ny,nz] , u_next[:,ny,nz] , D0[:,ny,nz] , D1[:,ny,nz] , ∆x[:,ny,nz] , D0[:,ny,nz] , dtt;   
            q_extrame = advStr_x.q_extrame, Flux_all = advStr_x.Flux_all, Qmax = advStr_x.Qmax,Qmin = advStr_x.Qmin,   Qout = advStr_x.Qout);
        end
        # for j direction
        for nx in axes(Q,1), nz in axes(Q,3)
            @views Q[nx,:,nz] =myAdvec1d( Q[nx,:,nz] , v_pre[nx,:,nz] , v_next[nx,:,nz] , D1[nx,:,nz] , D2[nx,:,nz] , ∆y[nx,:,nz] , D0[nx,:,nz] , dtt;
                    q_extrame = advStr_y.q_extrame, Flux_all = advStr_y.Flux_all, Qmax = advStr_y.Qmax,Qmin = advStr_y.Qmin,   Qout = advStr_y.Qout);
        end
        # for k direction
        for nx in axes(Q,1), ny in axes(Q,2)
            @views Q[nx,ny,:] =myAdvec1d( Q[nx,ny,:] , w_pre[nx,ny,:] , w_next[nx,ny,:] , D2[nx,ny,:] , D3[nx,ny,:] , ∆z[nx,ny,:] , D0[nx,ny,:] , dtt;
            q_extrame = advStr_z.q_extrame, Flux_all = advStr_z.Flux_all, Qmax = advStr_z.Qmax,Qmin = advStr_z.Qmin,   Qout = advStr_z.Qout);
        end
    c_out .= Q.*denOut;
    end
    # c_out = Q.*denIn;
    return c_out
end

function myAdvec1d( Qin::AbstractVector{T}, U_pre,U_next, DEN0, DEN1, ∆x, DD0, DT;LSHP = 0.0,
    q_extrame = Array{Bool}(undef,size(Qin)), Flux_all = Array{T}(undef,length(Qin)+1), Qmax = similar(Qin),Qmin = similar(Qin), Qout = similar(Qin), 
    fullSharpening = false) where T
    # Qin is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
    # DXX is grid spacing, DD0 is the original density, DT is the time step
    # Q_bound::T=0
    @views Flux_pre,Flux_next=Flux_all[1:end-1],Flux_all[2:end]

    # % identify local max and min, specigy mxing ratio limit at new time
    Qmin .= Qin;Qmax .= Qin;
    @inbounds for i in eachindex(Qin)   
        # set Q_left and Q_right
        if i==1
            Q_left,Q_right=Qin[i],Qin[i+1] # *boundary*
        elseif i==length(Qin)
            Q_left,Q_right=Qin[i-1],Qin[i] # *boundary*
        else
            Q_left,Q_right=Qin[i-1],Qin[i+1]
        end
        # set q_extrame
        if  Qin[i] > max(Q_left, Q_right) || Qin[i] < min(Q_left, Q_right)
            q_extrame[i]=true
        else
            q_extrame[i]=false
        end
        # set Qmax and Qmin
        if U_pre[i] > 0;        Qmin[i],Qmax[i] = min(Q_left,Qmin[i]),max(Q_left,Qmax[i]);    end #isUpwind_left
        if U_next[i] < 0;       Qmin[i],Qmax[i] = min(Q_right,Qmin[i]),max(Q_right,Qmax[i]);  end #isUpwind_right   
    end

    # % update mixing ratios and limit Flux going up where u>0   = * ->
    Flux_pre[1]=  Qin[1] *  U_pre[1] * DT * DD0[1] 
    @inbounds for i in eachindex(Qin)
        if U_next[i] < 0; continue;end
        if U_pre[i] < 0 # <- * ->
            Flux_next[i] = Qin[i] * U_next[i] * DT * DD0[i]; 
        else # -> * ->
            # set Q_left and Q_right
            # if i==length(Qin)
            #     Q_upwind,Q_downwind=Qin[i],Qin[i] # *boundary*
            #     isExtrame_upwind,isExtrame_downwind=q_extrame[i],false
            # else
            #     Q_upwind,Q_downwind=Qin[i],Qin[i+1]
            #     isExtrame_upwind,isExtrame_downwind=q_extrame[i],q_extrame[i+1]
            # end
            if i == 1
                Q_upwind,Q_downwind=Qin[i],Qin[i+1] # *boundary*
                isExtrame_upwind,isExtrame_downwind=false,q_extrame[i+1]
            elseif i == length(Qin)
                Q_upwind,Q_downwind=Qin[i-1],Qin[i]
                isExtrame_upwind,isExtrame_downwind=q_extrame[i-1],q_extrame[i]
            else
                Q_upwind,Q_downwind=Qin[i-1],Qin[i+1]
                isExtrame_upwind,isExtrame_downwind=q_extrame[i-1],q_extrame[i+1]
            end


            X1 = DT * U_next[i] / ∆x[i]; # raio of transport length to grid spacing
            CF1 = myCFCALC(X1, Qin[i], Q_upwind, Q_downwind, isExtrame_upwind, isExtrame_downwind, fullSharpening);
            # CF1 = CFCALC(X1,  Qin[i],Q_right, Q_left, isExtrame_right, isExtrame_left, LSHP);
            FLUXI = DD0[i] * U_next[i] * DT * CF1 / ∆x[i];

            Qout[i] = max(Qmin[i], min(Qmax[i], (Qin[i] * DEN0[i] - FLUXI + Flux_pre[i] / ∆x[i])/ DEN1[i]));
            Flux_next[i] = ∆x[i] * (Qin[i] * DEN0[i] - Qout[i] * DEN1[i]) + Flux_pre[i];
            # println(i," ",Qin[i]," ",Qout[i]," ",Flux_next[i]," ",Flux_pre[i]," ",FLUXI," ",Qmin[i]," ",Qmax[i])
        end
    end

    # % update mixing ratios and limit Flux going down where u<0
    if U_next[end] < 0; Flux_next[end]=  Qin[end] *  U_next[end] * DT * DD0[end]; end
    @inbounds for i in reverse(eachindex(Qin))
        if U_pre[i] > 0
            if U_next[i] < 0 # -> * <-
                Qout[i] = max(Qmin[i], min(Qmax[i], (Qin[i] * DEN0[i] - Flux_next[i] / ∆x[i] + Flux_pre[i] / ∆x[i]) / DEN1[i]));
                # println(i," ",Qout[i]," ",Qin[i]," ",Flux_next[i]," ",Flux_pre[i]," ",DEN1[i])
            end
        else # <- * ==
            # if i==1
            #     Q_upwind,Q_downwind=Qin[i],Qin[i]
            #     isExtrame_upwind,isExtrame_downwind=q_extrame[i],false
            #     DD0_pre=DD0[i]
            # else
            #     Q_upwind,Q_downwind=Qin[i],Qin[i-1]
            #     isExtrame_upwind,isExtrame_downwind=q_extrame[i],q_extrame[i-1]
            #     DD0_pre=DD0[i-1]
            # end
            
            if i==1
                Q_upwind,Q_downwind=Qin[i+1],Qin[i]
                isExtrame_upwind,isExtrame_downwind=q_extrame[i+1],q_extrame[i]
                DD0_pre=DD0[i]
            elseif i ==length(Qin)
                Q_upwind,Q_downwind=Qin[i],Qin[i-1]
                isExtrame_upwind,isExtrame_downwind=q_extrame[i],q_extrame[i-1]
                DD0_pre=DD0[i-1]
            else
                Q_upwind,Q_downwind=Qin[i+1],Qin[i-1]
                isExtrame_upwind,isExtrame_downwind=q_extrame[i+1],q_extrame[i-1]
                DD0_pre=DD0[i-1]
            end
            X1 = DT * abs(U_pre[i]) / ∆x[i];
            CF1 = myCFCALC(X1, Qin[i], Q_upwind, Q_downwind, isExtrame_upwind, isExtrame_downwind, fullSharpening);
            # CF1 = myCFCALC(X1, Qin[i], Q_left, Q_right, isExtrame_left, isExtrame_right, fullSharpening);
            # CF1 = CFCALC(X1,  Qin[i],Q_right, Q_left, isExtrame_right, isExtrame_left, LSHP);
            
            if U_next[i] > 0; CF1 = Qin[i]; end#% outflow only cell
            FLUXex1 = DD0_pre * U_pre[i] * DT * CF1 / ∆x[i];
            Qout[i] = max(Qmin[i], min(Qmax[i], (Qin[i] * DEN0[i] - Flux_next[i] / ∆x[i] + FLUXex1) / DEN1[i]));
            Flux_pre[i] = ∆x[i] * (Qout[i] * DEN1[i] - Qin[i] * DEN0[i]) + Flux_next[i];
        end
    end
    # println(round.(Flux_pre))
    # println(round.(Qin))
    # println(round.(Qout))
    # println(round.(Flux_next))
    return Qout
end

# function myAdvec1d( Qin::AbstractVector{T}, U_pre,U_next, DEN0, DEN1, ∆x, DD0, DT;LSHP = 0.0,
#     q_extrame = Array{Bool}(undef,size(Qin)), Flux_all = Array{T}(undef,length(Qin)+1), Qmax = similar(Qin),Qmin = similar(Qin), Qout = similar(Qin), 
#     fullSharpening = false) where T
#     # Qin is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
#     # DXX is grid spacing, DD0 is the original density, DT is the time step
#     # Q_bound::T=0
#     @views Flux_pre,Flux_next=Flux_all[1:end-1],Flux_all[2:end]

#     # % identify local max and min, specigy mxing ratio limit at new time
#     Qmin .= Qin;Qmax .= Qin;
#     @inbounds for i in eachindex(Qin)   
#         # set Q_left and Q_right
#         if i==1
#             Q_left,Q_right=0,Qin[i+1] # *boundary*
#         elseif i==length(Qin)
#             Q_left,Q_right=Qin[i-1],0 # *boundary*
#         else
#             Q_left,Q_right=Qin[i-1],Qin[i+1]
#         end
#         # set q_extrame
#         if  Qin[i] > max(Q_left, Q_right) || Qin[i] < min(Q_left, Q_right)
#             q_extrame[i]=true
#         else
#             q_extrame[i]=false
#         end
#         # set Qmax and Qmin
#         if U_pre[i] > 0;        Qmin[i],Qmax[i] = min(Q_left,Qmin[i]),max(Q_left,Qmax[i]);    end #isUpwind_left
#         if U_next[i] < 0;       Qmin[i],Qmax[i] = min(Q_right,Qmin[i]),max(Q_right,Qmax[i]);  end #isUpwind_right   
#     end

#     Flux_all .= NaN
#     #  update Flux using myCFCALC for -> * -> and <- * <-
#     @inbounds for i in eachindex(Qin)
#         # set Q_left and Q_right
#         if i==1
#             Q_left,Q_right=0,Qin[i+1] # *boundary*
#             isExtrame_left,isExtrame_right=false,q_extrame[i+1]
#         elseif i==length(Qin)
#             Q_left,Q_right=Qin[i-1],0 # *boundary*
#             isExtrame_left,isExtrame_right=q_extrame[i-1],false
#         else
#             Q_left,Q_right=Qin[i-1],Qin[i+1]
#             isExtrame_left,isExtrame_right=q_extrame[i-1],q_extrame[i+1]
#         end
#         # update Flux
#         if U_pre[i] > 0 
#             if U_next[i] > 0 # -> * ->
#                 X1 = DT * U_next[i] / ∆x[i]; # raio of transport length to grid spacing
#                 Q_crossFace = myCFCALC(X1, Qin[i], Q_left, Q_right, isExtrame_left, isExtrame_right, fullSharpening);
#                 Flux_next[i] = DD0[i] * U_next[i] * DT * Q_crossFace ;
#             end
#         else # <- * ==
#             if U_next[i] < 0 # <- * <-
#                 X1 = DT * U_pre[i] / ∆x[i]; # raio of transport length to grid spacing
#                 Q_crossFace = myCFCALC(X1, Qin[i],Q_right, Q_left,isExtrame_right,  isExtrame_left, fullSharpening);
#                 Flux_pre[i] = DD0[i] * U_pre[i] * DT * Q_crossFace;
#             end
#         end
#     end

#     # update remaining flux
#     if isnan(Flux_pre[1]); 
#         if U_pre[1] > 0; Flux_pre[1] = 0; # *boundary*
#         else;            Flux_pre[1] = Qin[1] * U_pre[1] * DT * DD0[1]; 
#         end
#     end
#     if isnan(Flux_next[end]); 
#         if U_next[end] < 0; Flux_next[end] = 0; # *boundary*
#         else;               Flux_next[end] = Qin[end] * U_next[end] * DT * DD0[end]; 
#         end
#     end
#     @inbounds for i in eachindex(Qin)[1:end-1]
#         if isnan(Flux_next[i])
#             if U_next[i] > 0 # === * ->
#                 Flux_next[i] = Qin[i] * U_next[i] * DT * DD0[i];
#             else # === * <-
#                 Flux_next[i] = Qin[i+1] * U_next[i] * DT * DD0[i+1];
#             end
#         end
#     end

#     # limit Flux and update Qout
#     @inbounds for i in eachindex(Qin)
#         Q_out_tmp =  (Qin[i] * DEN0[i] - Flux_next[i] / ∆x[i] + Flux_pre[i] / ∆x[i]) / DEN1[i];
#         if Q_out_tmp > Qmax[i]
#             Q_out_tmp = Qmax[i]
#             Flux_next[i] = ∆x[i] * (Qin[i] * DEN0[i] -Q_out_tmp * DEN1[i]) + Flux_pre[i];
#         elseif Q_out_tmp < Qmin[i]
#             Q_out_tmp = Qmin[i]
#             Flux_next[i] = ∆x[i] * (Qin[i] * DEN0[i] -Q_out_tmp * DEN1[i]) + Flux_pre[i];
#         end
#         Qout[i] = Q_out_tmp
#     end
#     return Qout
# end

function advec1d_mass( Qin::AbstractVector{T}, U_pre,U_next, DEN0, DEN1, ∆x, massFlux_all, DT;
    q_extrame = Array{Bool}(undef,size(Qin)), Flux_all = Array{T}(undef,length(Qin)+1), Qmax = similar(Qin),Qmin = similar(Qin), Qout = similar(Qin), 
    fullSharpening = false) where T
    # Qin is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
    # DXX is grid spacing, DD0 is the original density, DT is the time step
    # Q_bound::T=0
    @views Flux_pre,Flux_next=Flux_all[1:end-1],Flux_all[2:end]
    @views massFlux_all_pre,massFlux_all_next=massFlux_all[1:end-1],massFlux_all[2:end]

    # % identify local max and min, specigy mxing ratio limit at new time
    Qmin .= Qin;Qmax .= Qin;
    @inbounds for i in eachindex(Qin)   
        # set Q_left and Q_right
        if i==1
            Q_left,Q_right=0,Qin[i+1] # *boundary*
        elseif i==length(Qin)
            Q_left,Q_right=Qin[i-1],0 # *boundary*
        else
            Q_left,Q_right=Qin[i-1],Qin[i+1]
        end
        # set q_extrame
        if  Qin[i] > max(Q_left, Q_right) || Qin[i] < min(Q_left, Q_right)
            q_extrame[i]=true
        else
            q_extrame[i]=false
        end     
        # set Qmax and Qmin
        if U_pre[i] > 0;        Qmin[i],Qmax[i] = min(Q_left,Qmin[i]),max(Q_left,Qmax[i]);    end #isUpwind_left
        if U_next[i] < 0;       Qmin[i],Qmax[i] = min(Q_right,Qmin[i]),max(Q_right,Qmax[i]);  end #isUpwind_right   
    end

    Flux_all .= NaN
    #  update Flux using myCFCALC for -> * -> and <- * <-
    @inbounds for i in eachindex(Qin)
        # set Q_left and Q_right
        if i==1
            Q_left,Q_right=0,Qin[i+1] # *boundary*
            isExtrame_left,isExtrame_right=false,q_extrame[i+1]
        elseif i==length(Qin)
            Q_left,Q_right=Qin[i-1],0 # *boundary*
            isExtrame_left,isExtrame_right=q_extrame[i-1],false
        else
            Q_left,Q_right=Qin[i-1],Qin[i+1]
            isExtrame_left,isExtrame_right=q_extrame[i-1],q_extrame[i+1]
        end
        # update Flux
        if U_pre[i] > 0 
            if U_next[i] > 0 # -> * ->
                X1 = DT * U_next[i] / ∆x[i]; # raio of transport length to grid spacing
                Q_crossFace = myCFCALC(X1, Qin[i], Q_left, Q_right, isExtrame_left, isExtrame_right, fullSharpening);
                Flux_next[i] = massFlux_all_next[i] * Q_crossFace ;
            end
        else # <- * ==
            if U_next[i] < 0 # <- * <-
                X1 = DT * U_pre[i] / ∆x[i]; # raio of transport length to grid spacing
                Q_crossFace = myCFCALC(X1, Qin[i],Q_right, Q_left,isExtrame_right,  isExtrame_left, fullSharpening);
                Flux_pre[i] = massFlux_all_pre[i] * Q_crossFace;
            end
        end
    end

    # update remaining flux
    if isnan(Flux_pre[1]); 
        if U_pre[1] > 0; Flux_pre[1] = 0; # *boundary*
        else;            Flux_pre[1] = massFlux_all_pre[1] * Qin[1]; 
        end
    end
    if isnan(Flux_next[end]); 
        if U_next[end] < 0; Flux_next[end] = 0; # *boundary*
        else;               Flux_next[end] =massFlux_all_next[end] * Qin[end]; 
        end
    end
    @inbounds for i in eachindex(Qin)[1:end-1]
        if isnan(Flux_next[i])
            if U_next[i] > 0 # === * ->
                Flux_next[i] = Qin[i] * massFlux_all_next[i];
            else # === * <-
                Flux_next[i] = Qin[i+1] * massFlux_all_next[i];
            end
        end
    end

    # limit Flux and update Qout
    @inbounds for i in eachindex(Qin)
        Q_out_tmp =  (Qin[i] * DEN0[i] - Flux_next[i] / ∆x[i] + Flux_pre[i] / ∆x[i]) / DEN1[i];
        if Q_out_tmp > Qmax[i]
            Q_out_tmp = Qmax[i]
            Flux_next[i] = ∆x[i] * (Qin[i] * DEN0[i] -Q_out_tmp * DEN1[i]) + Flux_pre[i];
        elseif Q_out_tmp < Qmin[i]
            Q_out_tmp = Qmin[i]
            Flux_next[i] = ∆x[i] * (Qin[i] * DEN0[i] -Q_out_tmp * DEN1[i]) + Flux_pre[i];
        end
        Qout[i] = Q_out_tmp
    end
    return Qout
end

struct myAvecStr
    q_extrame
    Flux_all
    Qmax
    Qmin
    Qout
end
function init_myAdvec(Qin::AbstractVector{T})  where T    
    q_extrame,Flux_all=falses(length(Qin)),zeros(T,length(Qin)+1)
    Qmax,Qmin,Qout=zeros(T,length(Qin)),zeros(T,length(Qin)),zeros(T,length(Qin))
    return myAvecStr(q_extrame,Flux_all,Qmax,Qmin,Qout)
end

function adv_3d_wuhj(c_in::AbstractArray{T}, u, v,w, ∆x, ∆y,∆z, dt;par=false,num_tStep=-1) where T
    # input C_in is the mass concentration
    # set num_tStep=0 to calculate the num_tStep automatically by wind speed and grid size
    # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
    # parallel version--------------------------------------
    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_in)); end
    if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end
    function expendZ(in)
        out=cat(in[:,:,[1]],in,in[:,:,[end]],dims=3)
        return out
    end
    c_in, u, v,w, ∆x, ∆y,∆z = expendZ.([c_in, u, v,w, ∆x, ∆y,∆z])

    c_out=copy(c_in);
    
    DEN1 = zeros(T,size(c_in)); DEN1y = zeros(T,size(c_in)); DEN1z = zeros(T,size(c_in));

    d0 = ones(T,size(c_out));

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = min(minimum(abs.(∆x ./ u)), minimum(abs.(∆y ./ v)),minimum(abs.(∆z ./ w)));
        num_tStep = ceil(dt / dt0);
    end
    dtt=convert(T,dt / num_tStep);
    if !(num_tStep <50)
        @warn "too much num_tStep!($num_tStep) dtt=$dtt"
    end
    
    if par
        advecVar_all_dim1=[init_advec(c_in[:,1,1]) for i in 1:Threads.nthreads()]
        advecVar_all_dim2=[init_advec(c_in[1,:,1]) for i in 1:Threads.nthreads()]
        advecVar_all_dim3=[init_advec(c_in[1,1,:]) for i in 1:Threads.nthreads()]
    else
        advecVar_dim1=init_advec(c_in[:,1,1])
        advecVar_dim2=init_advec(c_in[1,:,1])
        advecVar_dim3=init_advec(c_in[1,1,:])
    end

    for itt = 1:num_tStep
        # --------------------------------update mixing ratio -------------------------------
        # for i direction
        Q0 = c_out;
        u1D = u;
        DXX = ∆x;
        DEN0 = d0;
        @views DEN1[2:end, :,:] = 1 .- dtt ./ DXX[2:end, :,:] .* (d0[2:end, :, :] .* u[2:end, :, :] .- d0[2:end, :, :] .* u[1:end - 1, :, :]);
        DEN1[1, :,:] .= 0;
        DD0 = DEN0;
        QN=c_out;
        if par == true
            Threads.@threads for ny = 2:(size(QN,2)-1)
                advecVar_curr=reset!(advecVar_all_dim1[Threads.threadid()])
                for nz in axes(QN,3)
                    @views QN[:,ny,nz] =advec1d( Q0[:,ny,nz] , u1D[:,ny,nz] , DEN0[:,ny,nz] , DEN1[:,ny,nz] , DXX[:,ny,nz] , DD0[:,ny,nz] , dtt;
                        IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                end
            end
        else
            for ny = 2:(size(QN,2)-1)
                advecVar_curr=advecVar_dim1
                for nz in axes(QN,3)
                    @views QN[:,ny,nz] =advec1d( Q0[:,ny,nz] , u1D[:,ny,nz] , DEN0[:,ny,nz] , DEN1[:,ny,nz] , DXX[:,ny,nz] , DD0[:,ny,nz] , dtt; 
                        IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                end
            end
        end
        c_out = QN;

        # for j direction        
        Q0y = c_out;
        u1Dy = v;
        DXXy = ∆y;
        DEN0y = DEN1;
        DD0y = d0;
        @views DEN1y[:, 1,:] = DEN0y[:, 1,:];
        # @views  DEN1y[:, 2:end] = DEN0y[:, 2:end] .- dtt ./ ∆y[:, 2:end] .* (d0[:, 2:end] .* v[:, 2:end] .- d0[:, 2:end] .* v[:, 1:end - 1]);
        @views  DEN1y[:, 2:end,:] = DEN0y[:, 2:end,:] .- dtt ./ ∆y[:, 2:end,:] .* (d0[:, 2:end,:] .* v[:, 2:end,:] .- d0[:, 2:end,:] .* v[:, 1:end - 1,:]);
        QNy=c_out;
        if par == true
            Threads.@threads for nx = 2:(size(QNy,1)-1)
                advecVar_curr=reset!(advecVar_all_dim2[Threads.threadid()])
                for nz in axes(QNy,3)
                    @views QNy[nx,:,nz] =advec1d( Q0y[nx,:,nz] , u1Dy[nx,:,nz] , DEN0y[nx,:,nz] , DEN1y[nx,:,nz] , DXXy[nx,:,nz] , DD0y[nx,:,nz] , dtt;                
                        IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                end
            end
        else
            for nx = 2:(size(QNy,1)-1)
                advecVar_curr=advecVar_dim2
                for nz in axes(QNy,3)
                    @views QNy[nx,:,nz] =advec1d( Q0y[nx,:,nz] , u1Dy[nx,:,nz] , DEN0y[nx,:,nz] , DEN1y[nx,:,nz] , DXXy[nx,:,nz] , DD0y[nx,:,nz] , dtt;                
                        IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                end
            end
        end
        c_out = QNy;
        # for k direction 
        Q0z = c_out;
        u1Dz = w;
        DXXz = ∆z;
        DEN0z = DEN1y;
        DD0z = d0;
        @views DEN1z[:, :, 1] = DEN0z[:, :, 1];
        @views DEN1z[:, :, 2:end] = DEN0z[:, :, 2:end] .- dtt ./ ∆z[:, :, 2:end] .* (d0[:, :, 2:end] .* w[:, :, 2:end] .- d0[:, :, 2:end] .* w[:, :, 1:end - 1]);
        QNz=c_out;
        if par == true
            # println("new advz")
            Threads.@threads for nx = 2:(size(QNz,1)-1)
                advecVar_curr=reset!(advecVar_all_dim3[Threads.threadid()])
                for ny = 2:(size(QNz,2)-1)
                    @views QNz[nx,ny,:] =advec1d( Q0z[nx,ny,:] , u1Dz[nx,ny,:] , DEN0z[nx,ny,:] , DEN1z[nx,ny,:] , DXXz[nx,ny,:] , DD0z[nx,ny,:] , dtt;                
                        IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                end
            end
        else
            for nx  = 2:(size(QNz,1)-1)
                advecVar_curr=advecVar_dim3
                for ny = 2:(size(QNz,2)-1)
                    @views QNz[nx,ny,:] =advec1d( Q0z[nx,ny,:] , u1Dz[nx,ny,:] , DEN0z[nx,ny,:] , DEN1z[nx,ny,:] , DXXz[nx,ny,:] , DD0z[nx,ny,:] , dtt;                
                        IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                end
            end
        end
        c_out = QNz;
        
        # % preserve the boundary
        @views c_out[[1 end],:,:]=c_in[[1 end],:,:];
        @views c_out[:,[1 end],:]=c_in[:,[1 end],:];
    end
    c_out=c_out[:,:,2:end-1]
    return c_out
end
function advection(c_in::AbstractArray{T}, u, v,w, ∆x, ∆y,∆z, dt;par=false,num_tStep=-1,denIn=ones(T,size(c_in)),denOut=ones(T,size(c_in)),massBalance=false,massInOut=false,constDen=true) where T
    # input C_in default to be volumn concentration like ppb, set massInOut to use mass concentration
    # set num_tStep=0 to calculate the num_tStep automatically by wind speed and grid size
    # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
    if !par     # non parallel version
        return adv_3d_onlineDen(c_in, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=num_tStep,massBalance=massBalance) 
        # adv_3d_onlineDen(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=1) 
        # return adv_3d_basic(c_in, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=num_tStep) 
    end
    # parallel version--------------------------------------
    Q=similar(c_in);
    c_out=copy(c_in);
    if size(u)==size(c_out); u=cat(u[[1],:,:],u;dims=1); end
    if size(v)==size(c_out); v=cat(v[:,[1],:],v;dims=2); end
    # if size(u)==size(c_out); u=cat(u,u[[end],:,:];dims=1); end
    # if size(v)==size(c_out); v=cat(v,v[:,[end],:];dims=2); end
    if size(w)==size(c_out); w=cat(w,w[:,:,[end]];dims=3); end
    @views u_pre,u_next=u[1:end-1,:,:],u[2:end,:,:];
    @views v_pre,v_next=v[:,1:end-1,:],v[:,2:end,:];
    @views w_pre,w_next=w[:,:,1:end-1],w[:,:,2:end];

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_in)); end
    if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = min(minimum(abs.(∆x ./ u_next)), minimum(abs.(∆y ./ v_next)),minimum(abs.(∆z ./ w_next)));
        num_tStep = ceil(dt / dt0);
    end
    dtt=convert(T,dt / num_tStep);
    if !(num_tStep <50)
        @warn "too much num_tStep!($num_tStep) dtt=$dtt"
    end

    advStr_all_dim1=[init_myAdvec(c_in[:,1,1]) for i in 1:Threads.nthreads()]
    advStr_all_dim2=[init_myAdvec(c_in[1,:,1]) for i in 1:Threads.nthreads()]
    advStr_all_dim3=[init_myAdvec(c_in[1,1,:]) for i in 1:Threads.nthreads()]

    # denStr_all_dim1=[Array{eltype(denIn),1}(undef,size(denIn,1)+1)  for i in 1:Threads.nthreads()]
    # denStr_all_dim2=[Array{eltype(denIn),1}(undef,size(denIn,2)+1)  for i in 1:Threads.nthreads()]
    # denStr_all_dim3=[Array{eltype(denIn),1}(undef,size(denIn,3)+1)  for i in 1:Threads.nthreads()]


    den_minimum=minimum(denIn)/1e3
    denCurr = copy(denIn);
    den2Change_aStep=(denOut .- denIn )./num_tStep
    D0=copy(denCurr); D1=copy(D0);D2=copy(D0);D3=copy(D0);
    # Q .= c_in./denIn;
    Q = massInOut ? c_in./denIn : c_in
    for itt = 1:num_tStep
        # deal with densities
        # if massInOut; Q = c_out./denCurr; end
        
        # ---------------------------------- update densities -------------------------------
        D0 .= constDen ? denIn : denCurr
        # d_D1=dtt .* (u_pre .* denIn - u_next .* denIn) ./ ∆x
        # D1 = D0 + d_D1
        # d_D2=dtt .* (v_pre .* denIn - v_next .* denIn) ./ ∆y
        # D2 = D1 + d_D2
        # d_D3=dtt .* (w_pre .* denIn - w_next .* denIn) ./ ∆z
        # D3= D2 + d_D3

        # den2adjust=den2Change_aStep .- (d_D1 .+ d_D2 .+ d_D3)
        # allAbsDenChange=abs.(d_D1) .+ abs.(d_D2) .+ abs.(d_D3)
        # q_noChange=allAbsDenChange .< den_minimum
        # allAbsDenChange[q_noChange] .= den_minimum
        # ratio_d1= abs.(d_D1) ./ allAbsDenChange
        # ratio_d2= abs.(d_D2) ./ allAbsDenChange
        # ratio_d3= abs.(d_D3) ./ allAbsDenChange
        # ratio_d1[q_noChange] .= 1/3;ratio_d2[q_noChange] .= 1/3;ratio_d3[q_noChange] .= 1/3;
        # # print(minimum(abs.(d_D1) .+ abs.(d_D2) .+ abs.(d_D3))," ")
        # d_D1 +=  den2adjust .* ratio_d1
        # d_D2 +=  den2adjust .* ratio_d2
        # d_D3 +=  den2adjust .* ratio_d3
        # D1= D0 .+ d_D1
        # D2= D1 .+ d_D2
        # D3= D2 .+ d_D3
        # D1 = D0 + dtt .* (u_pre .* D0 - u_next .* D0) ./ ∆x
        # D2 = D1 + dtt .* (v_pre .* D0 - v_next .* D0) ./ ∆y
        # D3 = D2 + dtt .* (w_pre .* D0 - w_next .* D0) ./ ∆z

        Threads.@threads         for ny in axes(D1,2)
           @views D1[:,ny,:]=D0[:,ny,:] + dtt .* (u_pre[:,ny,:] .* denIn[:,ny,:] - u_next[:,ny,:] .* denIn[:,ny,:]) ./ ∆x[:,ny,:]
            # massFlux_x=Array{T,1}(undef,size(D1,1)+1)
            # for nz in axes(D1,3)#updateDen1d_middelDen!     updateDen1d!
            #     @views D1[:,ny,nz]=updateDen1d_middelDen!(D0[:,ny,nz],D1[:,ny,nz],D0[:,ny,nz],u_pre[:,ny,nz],u_next[:,ny,nz],∆x[:,ny,nz],dtt,den_minimum; massFlux=massFlux_x)
            # end
        end
        # for j direction
        Threads.@threads          for nx  in axes(D2,1)
            @views D2[nx,:,:]=D1[nx,:,:] + dtt .* (v_pre[nx,:,:] .* denIn[nx,:,:] - v_next[nx,:,:] .* denIn[nx,:,:]) ./ ∆y[nx,:,:]
            # massFlux_y=Array{T,1}(undef,size(D2,2)+1)
            # for nz in axes(D2,3)
            #     # @views D2[nx,:,nz]=D1[nx,:,nz] .+ dtt .* (v_pre[nx,:,nz] .* D1[nx,:,nz] - v_next[nx,:,nz] .* D1[nx,:,nz]) ./ ∆y[nx,:,nz]
            #     @views D2[nx,:,nz]=updateDen1d_middelDen!(D1[nx,:,nz],D2[nx,:,nz],D0[nx,:,nz],v_pre[nx,:,nz],v_next[nx,:,nz],∆y[nx,:,nz],dtt,den_minimum; massFlux=massFlux_y)
            # end
        end
        # # D3= D0 +  den2Change_aStep
        # # for k direction
        # D3=D2 + dtt .* (w_pre .* denIn - w_next .* denIn) ./ ∆z
        
        # D3 = copy(D0)
        Threads.@threads   for nx in axes(D3,1)
            @views  D3[nx,:,:] = D2[nx,:,:] + dtt .* (w_pre[nx,:,:] .* denIn[nx,:,:] - w_next[nx,:,:] .* denIn[nx,:,:]) ./ ∆z[nx,:,:]
            # massFulx_z=Array{T,1}(undef,size(D3,3)+1)
            # for ny in axes(D3,2)
            #     # @views D3[nx,ny,:]=D2[nx,ny,:] .+ dtt .* (w_pre[nx,ny,:] .* D2[nx,ny,:] - w_next[nx,ny,:] .* D2[nx,ny,:]) ./ ∆z[nx,ny,:]
            #     @views D3[nx,ny,:]=updateDen1d_middelDen!(D2[nx,ny,:],D3[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt,den_minimum; massFlux=massFulx_z)
            # end
        end

        # D1 = D0 + dtt .* (u_pre .* denIn - u_next .* denIn) ./ ∆x
        # D2 = D1 + dtt .* (v_pre .* denIn - v_next .* denIn) ./ ∆y
        # D3 = D2 + dtt .* (w_pre .* denIn - w_next .* denIn) ./ ∆z

        # println("minimum Den ",minimum(D0)," ",minimum(D1)," ",minimum(D2)," ",minimum(D3))

        # # adjust density ----------------------        
        # Threads.@threads  for nx  in axes(den2Change_aStep,1)
        #     for ny in axes(den2Change_aStep,2)
        #         d_D1= D1[nx,ny,:]  - D0[nx,ny,:]
        #         d_D2= D2[nx,ny,:]  - D1[nx,ny,:]
        #         d_D3= D3[nx,ny,:]  - D2[nx,ny,:]
        #         den2adjust=den2Change_aStep[nx,ny,:] .- (d_D1 .+ d_D2 .+ d_D3)
        #         allAbsDenChange=abs.(d_D1) .+ abs.(d_D2) .+ abs.(d_D3)
        #         for (i,aChange) in enumerate(allAbsDenChange)
        #             if aChange < den_minimum
        #                 d_D1[i] += den2adjust[i] / 3
        #                 d_D2[i] += den2adjust[i] / 3
        #                 d_D3[i] += den2adjust[i] / 3
        #             else
        #                 d_D1[i] += den2adjust[i] * abs(d_D1[i]) / allAbsDenChange[i]
        #                 d_D2[i] += den2adjust[i] * abs(d_D2[i]) / allAbsDenChange[i]
        #                 d_D3[i] += den2adjust[i] * abs(d_D3[i]) / allAbsDenChange[i]
        #             end
        #         end
        #         D1[nx,ny,:] = D0[nx,ny,:] .+ d_D1
        #         D2[nx,ny,:] = D1[nx,ny,:] .+ d_D2
        #         D3[nx,ny,:] = D2[nx,ny,:] .+ d_D3
        #     end
        # end

        denCurr .=  D3
        # denCurr =  copy(D3)
        # --------------------------------update mixing ratio -------------------------------
        # for i direction
        Threads.@threads for ny in axes(Q,2)
            # advStr_curr=init_myAdvec(Q[:,ny,1])
            advStr_curr=advStr_all_dim1[Threads.threadid()]
            for nz in axes(Q,3)
                @views Q[:,ny,nz] =myAdvec1d( Q[:,ny,nz] , u_pre[:,ny,nz] , u_next[:,ny,nz] , D0[:,ny,nz] , D1[:,ny,nz] , ∆x[:,ny,nz] , D0[:,ny,nz] , dtt;            
                    q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
                # if minimum(Q[:,ny,nz]) <0 ;println(ny," ",nz," ", minimum(Q[:,ny,nz]));end
            end
        end
        # for j direction
        Threads.@threads  for nx in axes(Q,1)
            # advStr_curr=init_myAdvec(Q[nx,:,1])
            advStr_curr=advStr_all_dim2[Threads.threadid()]
            for nz in axes(Q,3)
                @views Q[nx,:,nz] =myAdvec1d( Q[nx,:,nz] , v_pre[nx,:,nz] , v_next[nx,:,nz] , D1[nx,:,nz] , D2[nx,:,nz] , ∆y[nx,:,nz] , D0[nx,:,nz] , dtt;
                    q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
            end
        end
        # for k direction 
        Threads.@threads  for nx in axes(Q,1)
            # advStr_curr=init_myAdvec(Q[nx,1,:])
            advStr_curr=advStr_all_dim3[Threads.threadid()]
            for ny in axes(Q,2)
                @views Q[nx,ny,:] =myAdvec1d( Q[nx,ny,:] , w_pre[nx,ny,:] , w_next[nx,ny,:] , D2[nx,ny,:] , D3[nx,ny,:] , ∆z[nx,ny,:] , D0[nx,ny,:] , dtt;
                    q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);     # fullSharpening = true, 
            end
        end
        # if massInOut; c_out = Q .* denCurr; end
    end
    c_out= massInOut ? Q .* denOut : Q
    # if !massInOut; c_out = Q;    end
    return c_out
end
# ================== update density only once and offline =============
# function advection(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;par=false,num_tStep=1,massBalance=true) where T
#     # input C_in is the mass concentration
#     # set num_tStep=0 to calculate the num_tStep automatically by wind speed and grid size
#     # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
#     if !par     # non parallel version
#         return adv_3d_onlineDen(c_in, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=num_tStep) 
#         # adv_3d_onlineDen(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=1) 
#         # return adv_3d_basic(c_in, denIn , u, v,w, ∆x, ∆y,∆z, dt;num_tStep=num_tStep) 
#     end
#     # parallel version--------------------------------------
#     c_out=copy(c_in);
#     Q=similar(c_out);
#     if size(u)==size(c_out); u=cat(u[[1],:,:],u;dims=1); end
#     if size(v)==size(c_out); v=cat(v[:,[1],:],v;dims=2); end
#     if size(w)==size(c_out); w=cat(w,w[:,:,[1]];dims=3); end
#     @views u_pre,u_next=u[1:end-1,:,:],u[2:end,:,:];
#     @views v_pre,v_next=v[:,1:end-1,:],v[:,2:end,:];
#     @views w_pre,w_next=w[:,:,1:end-1],w[:,:,2:end];


#     if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_out)); end
#     if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_out)); end
#     if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end

#     # divide dt into ddt
#     if !(num_tStep >0) # num_tStep <= 0 or isnan
#         dt0 = min(minimum(abs.(∆x ./ u_next)), minimum(abs.(∆y ./ v_next)),minimum(abs.(∆z ./ w_pre)));
#         # for x in [u_next,v_next,w_pre,∆x,∆y,∆z]; 
#         #     println(minimum(x),mean(x),maximum(x));  
#         # end
#         num_tStep = ceil(dt / dt0);
#     end
#     dtt=convert(T,dt / num_tStep);
#     if !(num_tStep <50)
#         @warn "too much num_tStep!($num_tStep) dtt=$dtt"
#     end
    
#     # update densities ========================
#     den_minimum=minimum(denIn)/1e3
    
#     # D0=copy(denIn); D1=copy(denIn);D2=copy(denIn);D3=copy(denIn);
#     D0=denIn; D1=similar(denIn);D2=similar(denIn);D3=similar(denIn);
#     @batch per=thread  for ny in axes(D1,2)
#         massFlux_x=Array{T,1}(undef,size(D1,1)+1)
#         for nz in axes(D1,3)
#             @views D1[:,ny,nz]=updateDen1d!(D0[:,ny,nz],D1[:,ny,nz],D0[:,ny,nz],u_pre[:,ny,nz],u_next[:,ny,nz],∆x[:,ny,nz],dtt,den_minimum; massFlux=massFlux_x)
#         end
#     end
#     # for j direction
#     @batch per=thread  for nx in axes(D2,1)
#         massFlux_y=Array{T,1}(undef,size(D2,2)+1)
#         for nz in axes(D2,3)
#             @views D2[nx,:,nz]=updateDen1d!(D1[nx,:,nz],D2[nx,:,nz],D0[nx,:,nz],v_pre[nx,:,nz],v_next[nx,:,nz],∆y[nx,:,nz],dtt,den_minimum; massFlux=massFlux_y)
#         end
#     end
#     # for k direction
#     @batch per=thread  for nx in axes(D3,1)
#         massFulx_z=Array{T,1}(undef,size(D3,3)+1)
#         for ny in axes(D3,2)
#             @views D3[nx,ny,:]=updateDen1d!(D2[nx,ny,:],D3[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt,den_minimum; massFlux=massFulx_z)
#         end
#     end
    
#     # D1 .= D0 .- (denIn .* u_next .- denIn .* u_pre) .* dtt ./ ∆x ;
#     # D2 .= D1 .- (denIn .* v_next .- denIn .* v_pre) .* dtt ./ ∆y ;
#     # D3 .= D2 .- (denIn .* w_next .- denIn .* w_pre) .* dtt ./ ∆z ;
#     denOut = D3;

#     # Q .= c_in./denIn;
#     if !massBalance; Q .= c_in; end
#     for itt = 1:num_tStep
#         # deal with densities
#         if massBalance; Q .= c_out./denIn; end
#         # for i direction
#         @batch per=thread  for ny in axes(Q,2)
#             advStr_curr=init_myAdvec(Q[:,ny,1])
#             for nz in axes(Q,3)
#                 @views Q[:,ny,nz] =myAdvec1d( Q[:,ny,nz] , u_pre[:,ny,nz] , u_next[:,ny,nz] , D0[:,ny,nz] , D1[:,ny,nz] , ∆x[:,ny,nz] , D0[:,ny,nz] , dtt;            
#                     q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
#             end
#         end
#         # for j direction
#         @batch per=thread  for nx in axes(Q,1)
#             advStr_curr=init_myAdvec(Q[nx,:,1])
#             for nz in axes(Q,3)
#                 @views Q[nx,:,nz] =myAdvec1d( Q[nx,:,nz] , v_pre[nx,:,nz] , v_next[nx,:,nz] , D1[nx,:,nz] , D2[nx,:,nz] , ∆y[nx,:,nz] , D0[nx,:,nz] , dtt;
#                     q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
#             end
#         end
#         # for k direction       
#         @batch per=thread  for nx in axes(Q,1)
#             advStr_curr=init_myAdvec(Q[nx,1,:])
#             for ny in axes(Q,2)
#                 @views Q[nx,ny,:] =myAdvec1d( Q[nx,ny,:] , w_pre[nx,ny,:] , w_next[nx,ny,:] , D2[nx,ny,:] , D3[nx,ny,:] , ∆z[nx,ny,:] , D0[nx,ny,:] , dtt;
#                     q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
#             end
#         end
#         if massBalance; c_out .= Q.*denOut; end
#     end
#     if !massBalance; c_out .= Q;    end
#     return c_out
# end
function advection_mass(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;par=false,num_tStep=1) where T
    # input C_in is the mass concentration
    # set num_tStep=0 to calculate the num_tStep automatically by wind speed and grid size
    # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
    # parallel version--------------------------------------
    c_out=copy(c_in);
    Q=similar(c_out);
    if size(u)==size(c_out); u=cat(u[[1],:,:],u;dims=1); end
    if size(v)==size(c_out); v=cat(v[:,[1],:],v;dims=2); end
    if size(w)==size(c_out); w=cat(w,w[:,:,[1]];dims=3); end
    massFlux_x=similar(u);massFlux_y=similar(v);massFlux_z=similar(w);
    @views u_pre,u_next=u[1:end-1,:,:],u[2:end,:,:];
    @views v_pre,v_next=v[:,1:end-1,:],v[:,2:end,:];
    @views w_pre,w_next=w[:,:,1:end-1],w[:,:,2:end];


    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_out)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_out)); end

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = min(minimum(abs.(∆x ./ u_next)), minimum(abs.(∆y ./ v_next)),minimum(abs.(∆z ./ w_pre)));
        num_tStep = ceil(dt / dt0);
    end
    dtt=convert(T,dt / num_tStep);
    # println("num_tStep=",num_tStep," dtt=",dtt)
    
    # update densities ========================
    den_minimum=minimum(denIn)/1e3
    
    # D0=copy(denIn); D1=copy(denIn);D2=copy(denIn);D3=copy(denIn);
    D0=denIn; D1=similar(denIn);D2=similar(denIn);D3=similar(denIn);
    @batch per=thread  for ny in axes(D1,2)
        for nz in axes(D1,3)
            @views updateDenFlux1d!(D1[:,ny,nz],massFlux_x[:,ny,nz],D0[:,ny,nz],D0[:,ny,nz],u_pre[:,ny,nz],u_next[:,ny,nz],∆x[:,ny,nz],dtt,den_minimum)
        end
    end
    # for j direction
    @batch per=thread  for nx in axes(D2,1)
        for nz in axes(D2,3)
            @views updateDenFlux1d!(D2[nx,:,nz],massFlux_y[nx,:,nz],D0[nx,:,nz],D0[nx,:,nz],v_pre[nx,:,nz],v_next[nx,:,nz],∆y[nx,:,nz],dtt,den_minimum)
        end
    end
    # for k direction
    @batch per=thread  for nx in axes(D3,1)
        for ny in axes(D3,2)
            @views updateDenFlux1d!(D3[nx,ny,:],massFlux_z[nx,ny,:],D0[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt,den_minimum)
        end
    end

    denOut = D3;

    # Q .= c_in./denIn;
    for itt = 1:num_tStep
        # deal with densities
        Q .= c_out./denIn;
        # for i direction
        @batch per=thread  for ny in axes(Q,2)
            advStr_curr=init_myAdvec(Q[:,ny,1])
            for nz in axes(Q,3)
                @views Q[:,ny,nz] =advec1d_mass( Q[:,ny,nz] , u_pre[:,ny,nz] , u_next[:,ny,nz] , D0[:,ny,nz] , D1[:,ny,nz] , ∆x[:,ny,nz] , massFlux_x[:,ny,nz] , dtt;            
                    q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
            end
        end
        # for j direction
        @batch per=thread  for nx in axes(Q,1)
            advStr_curr=init_myAdvec(Q[nx,:,1])
            for nz in axes(Q,3)
                @views Q[nx,:,nz] =advec1d_mass( Q[nx,:,nz] , v_pre[nx,:,nz] , v_next[nx,:,nz] , D1[nx,:,nz] , D2[nx,:,nz] , ∆y[nx,:,nz] , massFlux_y[nx,:,nz] , dtt;
                    q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
            end
        end
        # for k direction
        @batch per=thread  for nx in axes(Q,1)
            advStr_curr=init_myAdvec(Q[nx,1,:])
            for ny in axes(Q,2)
                @views Q[nx,ny,:] =advec1d_mass( Q[nx,ny,:] , w_pre[nx,ny,:] , w_next[nx,ny,:] , D2[nx,ny,:] , D3[nx,ny,:] , ∆z[nx,ny,:] , massFlux_z[nx,ny,:] , dtt;
                    q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
            end
        end
        c_out .= Q.*denOut;
    end
    return c_out
end
function adv_z_par!(Q_yz,w_pre_yz,w_next_yz,D2_yz,D3_yz,∆z_yz,D0_yz,dtt,advStr_curr)
    local ind_pre,ind_next
    for ny in axes(Q_yz,1)
        for nz in axes(D3_yz,2)
            if w_pre_yz[ny,nz]  > 0   &&  nz > 1;           ind_pre=nz-1;else;ind_pre=nz; end
            if w_next_yz[ny,nz] < 0   &&  nz < size(D3_yz,2);  ind_next=nz+1;else;ind_next=nz; end
            D3_yz[ny,nz]=D2_yz[ny,nz] - (D0_yz[ny,ind_next] * w_next_yz[ny,nz] - D0_yz[ny,ind_pre] * w_pre_yz[ny,nz]) * dtt / ∆z_yz[ny,nz]
        end
        @views Q_yz[ny,:] =myAdvec1d( Q_yz[ny,:] , w_pre_yz[ny,:] , w_next_yz[ny,:] , D2_yz[ny,:] , D3_yz[ny,:] , ∆z_yz[ny,:] , D0_yz[ny,:] , dtt;
            q_extrame = advStr_curr.q_extrame, Flux_all = advStr_curr.Flux_all, Qmax = advStr_curr.Qmax,Qmin = advStr_curr.Qmin,   Qout = advStr_curr.Qout);
    end
    return Q_yz,D3_yz
end

# update the density while ensure positive density and mass balance ===============================
function updateDenFlux1d!(Den_after, massFlux,Den_before::AbstractArray{T},D0,u_pre,u_next,gridSpace,dtt,den_minimum) where T
    @views massFlux_pre,massFlux_next=massFlux[1:end-1],massFlux[2:end]
    # update the mass flux ===============================
    massFlux_pre[1]=D0[1] * u_pre[1] * dtt
    massFlux_next[end]=D0[end] * u_next[end] * dtt
    @inbounds for i in eachindex(Den_after)[1:end-1]
        if u_next[i] < 0 
            massFlux_next[i]=D0[i+1] * u_next[i] * dtt 
        else
            massFlux_next[i]=D0[i] * u_next[i] * dtt 
        end
    end
    # limit the flux ===============================
    @inbounds for i in eachindex(Den_after)
        den_tmp=Den_before[i] + (massFlux_pre[i] - massFlux_next[i]) / gridSpace[i]
        if den_tmp < den_minimum
            # Den_after[i]=den_minimum
            Flux2adjust= (den_minimum -den_tmp)*gridSpace[i]
            massFlux_next[i] -= Flux2adjust
        # else
        #     Den_after[i]=den_tmp
        end
    end
    Den_after .= Den_before .+ (massFlux_pre - massFlux_next) ./ gridSpace
    
    return Den_after,massFlux
end
function inferWindFromDen!(Den_before::AbstractArray{T},Den_after,D0,u_pre,u_next,gridSpace,dtt;
    massFlux=Array{T,1}(undef,length(Den_next)+1)) where T    
    @views Flux_pre,Flux_next=massFlux[1:end-1],massFlux[2:end]
    # update the mass flux ===============================
    Flux_pre[1]=D0[1] * u_pre[1] * dtt
    # Flux_next[end]=D0[end] * u_next[end] * dtt
    @inbounds for i in eachindex(u_next)[1:end-1]
        Flux_next[i]=Flux_pre[i] - (Den_after[i] - Den_before[i]) * gridSpace[i]
        if Flux_next[i] > 0 
            u_next=Flux_next[i] / D0[i] / dtt
        else
            u_next=Flux_next[i] / D0[i+1] / dtt
        end
    end
    Flux_end=Flux_pre[end] - (Den_after[end] - Den_before[end]) * gridSpace[end]
    u_next[end]=Flux_end / D0[i] / dtt

    # update the wind speed =============================== 
    return u_pre
        # return u_pre,u_next
end
function getWFromDensity(denIn,denOut, u::AbstractArray{T}, v, ∆x, ∆y,∆z, dtt)   where T 
    
    if size(u)==size(denIn); u=cat(u[[1],:,:],u;dims=1); end
    if size(v)==size(denIn); v=cat(v[:,[1],:],v;dims=2); end
    w=zeros(T,size(denIn,1),size(denIn,2),size(denIn,3)+1)
    @views u_pre,u_next=u[1:end-1,:,:],u[2:end,:,:]
    @views v_pre,v_next=v[:,1:end-1,:],v[:,2:end,:]
    @views w_pre,w_next=w[:,:,1:end-1],w[:,:,2:end]

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_in)); end
    if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end

    den_minimum=minimum(denIn)/1e3
    D0=denIn; D1=similar(denIn);D2=similar(denIn);D3=denOut;
    # for i direction  
    for ny in axes(D1,2)
        massFlux_x=Array{T,1}(undef,size(D1,1)+1)
        for nz in axes(D1,3)
            @views D1[:,ny,nz]=updateDen1d!(D0[:,ny,nz],D1[:,ny,nz],D0[:,ny,nz],u_pre[:,ny,nz],u_next[:,ny,nz],∆x[:,ny,nz],dtt,den_minimum; massFlux=massFlux_x)
        end
    end
    # for j direction
    for nx in axes(D2,1)
        massFlux_y=Array{T,1}(undef,size(D2,2)+1)
        for nz in axes(D2,3)
            @views D2[nx,:,nz]=updateDen1d!(D1[nx,:,nz],D2[nx,:,nz],D0[nx,:,nz],v_pre[nx,:,nz],v_next[nx,:,nz],∆y[nx,:,nz],dtt,den_minimum; massFlux=massFlux_y)
        end
    end
    # for k direction
    for nx in axes(D3,1)
        massFulx_z=Array{T,1}(undef,size(D3,3)+1)
        for ny in axes(D3,2)
            @views w_pre[nx,ny,:]=inferWindFromDen!(D2[nx,ny,:],D3[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt; massFlux=massFulx_z) 
            # @views w_pre[nx,ny,:],w_next[nx,ny,:]=inferWindFromDen!(D2[nx,ny,:],D3[nx,ny,:],D0[nx,ny,:],w_pre[nx,ny,:],w_next[nx,ny,:],∆z[nx,ny,:],dtt; massFlux=massFulx_z) 
        end
    end
    return w_pre,w_next
end
function updateDen1d_middelDen!(Den_before::AbstractArray{T},Den_after,D0,u_pre,u_next,gridSpace,dtt,den_minimum;
    massFlux=Array{T,1}(undef,length(Den_before)+1)) where T
    @views Flux_pre,Flux_next=massFlux[1:end-1],massFlux[2:end]
    # update the mass flux ===============================
    @inbounds Flux_pre[1]=D0[1] * u_pre[1] * dtt
    @inbounds Flux_next[end]=D0[end] * u_next[end] * dtt
    @inbounds for i in eachindex(Den_after)[1:end-1]
        Flux_next[i]=(D0[i]+D0[i+1])/2 * u_next[i] * dtt 
    end
    @inbounds for i in eachindex(Den_after)
        Den_after[i]=Den_before[i] + (Flux_pre[i] - Flux_next[i]) / gridSpace[i]
    end
    # Flux_pre[1]=D0[1] * u_pre[1] * dtt
    # Flux_next[end]=D0[end] * u_next[end] * dtt
    # @views Flux_pre[2:end]=(D0[1:end-1] .+ D0[2:end]) ./2  .* u_pre[2:end] * dtt
    # Den_after .=  Den_before .+ (Flux_pre .- Flux_next) ./ gridSpace

    return Den_after
end
function updateDen1d!(Den_before::AbstractArray{T},Den_after,D0,u_pre,u_next,gridSpace,dtt,den_minimum;
    massFlux=Array{T,1}(undef,length(Den_before)+1)) where T
    @views Flux_pre,Flux_next=massFlux[1:end-1],massFlux[2:end]
    # update the mass flux ===============================
    Flux_pre[1]=D0[1] * u_pre[1] * dtt
    Flux_next[end]=D0[end] * u_next[end] * dtt
    @inbounds for i in eachindex(Den_after)[1:end-1]
        if u_next[i] < 0 
            Flux_next[i]=D0[i+1] * u_next[i] * dtt 
        else
            Flux_next[i]=D0[i] * u_next[i] * dtt 
        end
    end
    # D1= D0 .+ dtt ./ ∆x .* (D0 .* u_pre .- D0 .* u_next)
    # limit the flux ===============================
    @inbounds for i in eachindex(Den_after)
        den_tmp=Den_before[i] + (Flux_pre[i] - Flux_next[i]) / gridSpace[i]
        if den_tmp < den_minimum
            Den_after[i]=den_minimum
            Flux2adjust= (den_minimum -den_tmp)*gridSpace[i]
            Flux_next[i] =Flux_next[i] - Flux2adjust
        else
            Den_after[i]=den_tmp
        end
    end

    # update the density =============================== 
    # Den_after .= Den_before .+ dtt ./ gridSpace .* (D0 .* u_pre .- D0 .* u_next) 
    return Den_after
end
function updateDensities!(D0,D1,D2,D3,u_pre,u_next,v_pre,v_next,w_pre,w_next,∆x,∆y,∆z,dtt)
    
    # D1 .= D0 .- (denOut .* u_next .- denOut .* u_pre) .* dtt ./ ∆x ;
    # D2 .= D1 .- (denOut .* v_next .- denOut .* v_pre) .* dtt ./ ∆y ;
    # D3 .= D2 .- (denOut .* w_next .- denOut .* w_pre) .* dtt ./ ∆z ;

    # for D1
    local ind_pre,ind_next
    @inbounds for i in axes(D1,1), j in axes(D1,2), k in axes(D1,3)
        if u_pre[i,j,k]  > 0   &&  i > 1;           ind_pre=i-1;else;ind_pre=i; end
        if u_next[i,j,k] < 0   &&  i < size(D1,1);  ind_next=i+1;else;ind_next=i; end
        D1[i,j,k] = D0[i,j,k] - (D0[ind_next,j,k] * u_next[i,j,k] - D0[ind_pre,j,k] * u_pre[i,j,k]) * dtt / ∆x[i,j,k]
    end
    # for D2
    @inbounds for i in axes(D2,1), j in axes(D2,2), k in axes(D2,3)
        if v_pre[i,j,k]  > 0   &&  j > 1;           ind_pre=j-1;else;ind_pre=j; end
        if v_next[i,j,k] < 0   &&  j < size(D2,2);  ind_next=j+1;else;ind_next=j; end
        D2[i,j,k] = D1[i,j,k] - (D0[i,ind_next,k] * v_next[i,j,k] - D0[i,ind_pre,k] * v_pre[i,j,k]) * dtt / ∆y[i,j,k]
    end
    # for D3
    @inbounds for i in axes(D3,1), j in axes(D3,2), k in axes(D3,3)
        if w_pre[i,j,k]  > 0   &&  k > 1;           ind_pre=k-1;else;ind_pre=k; end
        if w_next[i,j,k] < 0   &&  k < size(D3,3);  ind_next=k+1;else;ind_next=k; end
        D3[i,j,k] = D2[i,j,k] - (D0[i,j,ind_next] * w_next[i,j,k] - D0[i,j,ind_pre] * w_pre[i,j,k]) * dtt / ∆z[i,j,k]
    end
    return D1,D2,D3
end

function advDiffAnHour(in::AbstractArray{T}, u, v, ∆x, ∆y;par=false,num_tStep=NaN,disp=false) where T

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(in)); end
    dt=60*60;
    transport2grid_max=1;
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        # numSub=ceil(max(u.*dt./dx...,v.*dt./dy...)/transport2grid_max);
        # num_tStep=ceil(max(maximum(u.*dt./dx),maximum(v.*dt./dy))/transport2grid_max);
        
        dt0 = transport2grid_max * minimum(abs.(vcat(vec(∆x ./ u), vec(∆y ./ v))));# seconds needed to traval one grid
        num_tStep = ceil(dt / dt0);
    end
    if disp; println("num_tStep=",num_tStep); end
    # println("num_tStep=",num_tStep)
    # out=adv_hori_wuhj( in, u, v, dx, dy,dt;par=par,num_tStep=num_tStep); 
    out=copy(in);
    dt_sub=dt/num_tStep;
    for _=1:num_tStep
        out=adv_hori_wuhj( out, u, v, ∆x, ∆y,dt_sub;par=par,num_tStep=1); 
        # out=adv_hori_basic( out, u, v, dx, dy,dt_sub;par=par,divideT=false); 
        out=diffusion_hori_Smagorinsky(out,u,v,∆x,∆y,dt_sub);
    end
    return out
end
function advDiffAnHour(c_in::AbstractArray{T}, denIn , u, v,w, ∆x, ∆y,∆z, dt;par=true,num_tStep=-1,massBalance=false,disp=false) where T

    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_in)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_in)); end
    if 1 == length(∆z); ∆z = ∆z .+ zeros(T,size(c_in)); end
    dt=60*60;
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = min(minimum(abs.(∆x ./ u)), minimum(abs.(∆y ./ v)),minimum(abs.(∆z ./ w)));
        num_tStep = ceil(dt / dt0);
    end
    if disp; println("num_tStep=",num_tStep); end
    c_out=copy(c_in);
    dt_sub=dt/num_tStep;
    for _=1:num_tStep
        # advection ------------------------
        c_out=advection(c_out, denIn , u, v,w, ∆x, ∆y,∆z, dt_sub;par=par,num_tStep=1,massBalance=massBalance)
        # diffusion ------------------------
        if par
            @batch per=thread  for nz in axes(c_out,3)
                @views c_out[:,:,nz]=diffusion_hori_Smagorinsky(c_out[:,:,nz],u[:,:,nz],v[:,:,nz],∆x[:,:,nz],∆y[:,:,nz],dt_sub);
                # @views c_out[:,:,nz]=diff_oneStep(c_out[:,:,nz],u[:,:,nz],v[:,:,nz],∆x[:,:,nz],∆y[:,:,nz],dt_sub);
            end
        else
            for nz in axes(c_out,3)
                @views c_out[:,:,nz]=diffusion_hori_Smagorinsky(c_out[:,:,nz],u[:,:,nz],v[:,:,nz],∆x[:,:,nz],∆y[:,:,nz],dt_sub);
            end
        end
    end
    return c_out
end