
# module TransportM
# export advDiffAnHour

export diffusion_hori_Smagorinsky,adv_hori_wuhj,adv_hori_basic,advDiffAnHour,adv_3d_basic
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
function diff_oneStep(in,u,v,dx,dy,dt)
    # % method in statistic paper
    ddIn=∇²(in,dy,dx);
    # % add boundary limit
    # %ddIn=boundaryMax(ddIn);
    dudx,dudy=∇(u,dx,dy)
    dvdx,dvdy=∇(v,dx,dy)
    K=@. (0.28)^2*dx*dy*sqrt( (dudx- dvdy)^2 + (dudy+dvdx)^2 );
    out=in.+K.*ddIn.*dt;
    return out,ddIn
end

function diffusion_hori_Smagorinsky(in,u,v,dx,dy,dt)
    minStepRatio=0.2;
    out,ddIn=diff_oneStep(in,u,v,dx,dy,dt);
    if 0==mean(abs.(ddIn.*dx.*dy)); return;end#% advoid 0 in the denominator
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
# struct advecStr
#     IMXMN::Vector{Float32}
#     FLUX::Vector{Float32}
#     VCMAX::Vector{Float32}
#     VCMIN::Vector{Float32}
#     QN::Vector{Float32}
# end
struct advecStr
    IMXMN
    FLUX
    VCMAX
    VCMIN
    QN
end
# struct advecStr
#     IMXMN::Vector{<:Real}
#     FLUX::Vector{<:Real}
#     VCMAX::Vector{<:Real}
#     VCMIN::Vector{<:Real}
#     QN::Vector{<:Real}
# end
function reset!(advecVar::advecStr)
    advecVar.IMXMN .= 0
    advecVar.FLUX .= 0
    advecVar.VCMAX .= 0
    advecVar.VCMIN .= 0
    advecVar.QN .= 0
    return advecVar
end
function init_advec(Q0::AbstractVector{T})  where T
    sx = 1 + 1;  ex = length(Q0) - 1; 
    
    IMXMN,FLUX,VCMAX,VCMIN,QN=zeros(T,size(Q0)),zeros(T,size(Q0)),zeros(T,size(Q0)),zeros(T,size(Q0)),zeros(T,size(Q0))
    return advecStr(IMXMN,FLUX,VCMAX,VCMIN,QN)
end
# Q0=rand(10)
# advecVar=init_advec(Q0)

# @time IMXMN,FLUX,VCMAX,VCMIN,QN=advecVar.IMXMN,advecVar.FLUX,advecVar.VCMAX,advecVar.VCMIN,advecVar.QN
# advecVar...
# advecVar.QN
# functions for advection ******************************************
# function advec1d( Q0::AbstractVector{T}, U, DEN0, DEN1, DXX, DD0, DT;advecVar=init_advec(Q0))  where T
#     sx = 1 + 1;  ex = length(Q0) - 1; 
    # IMXMN,FLUX,VCMAX,VCMIN,QN=Ref{AbstractVector{T}}(),Ref{AbstractVector{T}}(),Ref{AbstractVector{T}}(),Ref{AbstractVector{T}}(),Ref{AbstractVector{T}}()
    # IMXMN[],FLUX[],VCMAX[],VCMIN[],QN[]=advecVar.IMXMN,advecVar.FLUX,advecVar.VCMAX,advecVar.VCMIN,advecVar.QN
    # IMXMN,FLUX,VCMAX,VCMIN,QN=advecVar.IMXMN,advecVar.FLUX,advecVar.VCMAX,advecVar.VCMIN,advecVar.QN
    # IMXMN,FLUX,VCMAX,VCMIN,QN=Ref([advecVar.IMXMN,advecVar.FLUX,advecVar.VCMAX,advecVar.VCMIN,advecVar.QN])...
    # IMXMN=Ref(advecVar.IMXMN)
    # FLUX=Ref(advecVar.FLUX)
    # VCMAX=Ref(advecVar.VCMAX)
    # VCMIN=Ref(advecVar.VCMIN)
    # QN=Ref(advecVar.QN)
    
    
    # function advec1d( Q0, U, DEN0, DEN1, DXX, DD0, DT)
    # sx = 1 + 1;  ex = length(Q0) - 1; 
    # IMXMN = zeros(ex - sx + 3, 1); FLUX = zeros(ex - sx + 3, 1); VCMAX = zeros(ex - sx + 3, 1); VCMIN = zeros(ex - sx + 3, 1);
    # QN = zeros(size(Q0));

    function advec1d( Q0::AbstractVector{T}, U, DEN0, DEN1, DXX, DD0, DT;
        IMXMN = zeros(T,size(Q0)), FLUX = zeros(T,size(Q0)), VCMAX = zeros(T,size(Q0)),VCMIN = zeros(T,size(Q0)),
        QN = zeros(T,size(Q0))
    ) where T
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

function myAdvec1d( Q0::AbstractVector{T}, U_next,U_pre, DEN0, DEN1, DXX, DD0, DT;
    q_locMaxOrMin = falses(size(Q0)), Flux_all = zeros(T,size(Q0)+1), Q_max = similar(Q0),Q_min = similar(Q0),
    QN = zeros(T,size(Q0))
) where T
# Q0 is mixing ratio, U is wind speed, DEN0 is density at previous time step, DEN1 is density at current time step, 
# DXX is grid spacing, DD0 is the original density, DT is the time step
ZR0 = 0.0;
LSHP = 0.0;
Flux_pre=view(FLUX_all,1:lenght(Flux_all)-1)
Flux_next=view(FLUX_all,2:lenght(Flux_all))

# % identify local max and min, specigy mxing ratio limit at new time
q_locMaxOrMin[2:end-1]= (Q0[2:end-1] .> max.(Q0[1:end-2], Q0[3:end]) .| Q0[2:end-1] .< min.(Q0[1:end-2], Q0[3:end]))
Q_transported=similar(Q0)
for i in eachindex(Q0)    
    if U_next[i] < 0 && i<length(Q0)
        Q_upwind_curr=Q0[i + 1]
    elseif U_pre[i] > 0 && i>1
        Q_upwind_curr=Q0[I - 1]
    else
        Q_upwind_curr=Q0[i]
    end
    Q_transported[i]=Q_upwind_curr
    Q_max[i] = max(Q0[i], Q_upwind_curr);
    Q_min[i] = min(Q0[i], Q_upwind_curr);
end

# % update mixing ratios and limit Flux going up where u_pre > 0
if U_pre[1] >0; Flux_pre[1] = Q0[1] * U_pre[1] * DT * DD0[1]; end
for i in eachindex(Q0)
    if U_pre[i] >0
        if U_next[i] < 0
            FLUX[I] = Q0[I] * U[I] * DT * DD0[I];# % outflow only cell
        else
            X1 = DT * U[I] / DXX[I];
            CF1 = CFCALC(X1, Q0[I], Q0[I - 1], Q0[I + 1], IMXMN[I - 1], IMXMN[I + 1], LSHP);
            FLUXI = DD0[I] * U[I] * DT * CF1 / DXX[I];
            QN[I] = max(VCMIN[I], min(VCMAX[I], (Q0[I] * DEN0[I] - FLUXI + FLUX[I - 1] / DXX[I]) / DEN1[I]));
            FLUX[I] = DXX[I] * (Q0[I] * DEN0[I] - QN[I] * DEN1[I]) + FLUX[I - 1];
        end
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

function advec1d!( Q0::AbstractVector{T}, U, DEN0, DEN1, DXX, DD0, DT;pre::advecStr=init_advec(Q0)) where T
sx = 1 + 1; ex = length(Q0) - 1;
ZR0 = 0.0;
LSHP = 0.0;

# % identify local max and min, specigy mxing ratio limit at new time

for I = sx:ex
    pre.IMXMN[I] = 0.0;

    if (Q0[I] > max(Q0[I - 1], Q0[I + 1]) || Q0[I] < min(Q0[I - 1], Q0[I + 1]))
        pre.IMXMN[I] = 1;
    end

    CK1 = Q0[I];
    if (U[I] < ZR0); CK1 = Q0[I + 1]; end
    if (U[I - 1] > ZR0); CK1 = Q0[I - 1]; end
    DEN1[I] = max(DEN1[I], 1.E-20);
    DXX[I] = max(DXX[I], 1.E-20);
    pre.VCMAX[I] = max(Q0[I], CK1);
    pre.VCMIN[I] = min(Q0[I], CK1);
end

# % update mixing ratios and limit pre.FLUX going up where u>0
if (U[sx - 1] > ZR0); pre.FLUX[sx - 1] = Q0[sx - 1] * U[sx - 1] * DT * DD0[sx - 1]; end

for I = sx:ex
    if (U[I] < ZR0); continue; end

    if (U[I - 1] < ZR0)
        pre.FLUX[I] = Q0[I] * U[I] * DT * DD0[I];# % outflow only cell
    else
        X1 = DT * U[I] / DXX[I];
        CF1 = CFCALC(X1, Q0[I], Q0[I - 1], Q0[I + 1], pre.IMXMN[I - 1], pre.IMXMN[I + 1], LSHP);
        FLUXI = DD0[I] * U[I] * DT * CF1 / DXX[I];
        pre.QN[I] = max(pre.VCMIN[I], min(pre.VCMAX[I], (Q0[I] * DEN0[I] - FLUXI + pre.FLUX[I - 1] / DXX[I]) / DEN1[I]));
        pre.FLUX[I] = DXX[I] * (Q0[I] * DEN0[I] - pre.QN[I] * DEN1[I]) + pre.FLUX[I - 1];
    end

end

# % update mixing ratios and limit pre.FLUX going down where u<0
if (U[ex] < ZR0); pre.FLUX[ex] = Q0[ex + 1] * U[ex] * DT * DD0[ex]; end
for I = ex:-1:sx
    if (U[I - 1] > ZR0)
        if (U[I] < ZR0)#% inflow-only cell
            pre.QN[I] = (Q0[I] * DEN0[I] - pre.FLUX[I] / DXX[I] + pre.FLUX[I - 1] / DXX[I]) / DEN1[I];
        end

    else
        X1 = DT * abs(U[I - 1]) / DXX[I];
        CF1 = CFCALC(X1, Q0[I], Q0[I + 1], Q0[I - 1], pre.IMXMN[I + 1], pre.IMXMN[I - 1], LSHP);

        if (U[I] > ZR0); CF1 = Q0[I]; end#% outflow only cell
        FLUXex1 = DD0[I - 1] * U[I - 1] * DT * CF1 / DXX[I];
        pre.QN[I] = max(pre.VCMIN[I], min(pre.VCMAX[I], (Q0[I] * DEN0[I] - pre.FLUX[I] / DXX[I] + FLUXex1) / DEN1[I]));
        pre.FLUX[I - 1] = DXX[I] * (pre.QN[I] * DEN1[I] - Q0[I] * DEN0[I]) + pre.FLUX[I];
    end

end
return pre.QN
end

function CFCALC(X1, VCUP1, VCUP2, VCDW1, IMXUP2, IMXDN1, LSHP)
    # % function to calculate mixing ratio in the fuild moving into neighnor
    # % cell durint DT; X1=Courant No. set LSHP=1 for full sharpening
    
    if LSHP == 1
        ALFA = 5;
    else
        ALFA = 1;
        if IMXDN1 > 0; ALFA = 1.7763 - 0.5 * X1; end
        if IMXUP2 > 0; ALFA = 1.2578 + 0.58502 * X1; end
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
        # if false
        #     # if par == true
        #     Threads.@threads  for i=2:size(DEN1,1)
        #         @views DEN1[i, :] = 1 .- dtt ./ DXX[i, :] .* (d0[i, :] .* u[i, :] .- d0[i, :] .* u[i-1, :]);
        #     end
        # else
            @views DEN1[2:end, :] = 1 .- dtt ./ DXX[2:end, :] .* (d0[2:end, :] .* u[2:end, :] .- d0[2:end, :] .* u[1:end - 1, :]);
        # end
        DEN1[1, :] .= 0;
        DD0 = DEN0;
        QN=c_out;
        if par == true
            Threads.@threads for ny = 2:(size(QN,2)-1)
                # advecVar_curr=view(advecVar_all_dim1,Threads.threadid())
                advecVar_curr=reset!(advecVar_all_dim1[Threads.threadid()])
                # advecVar_curr=advecVar_all_dim1[Threads.threadid()]
                # advecVar_curr=advecVar_dim1
                @views QN[:,ny] =advec1d( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt;                
                    IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
            end
        else
            for ny = 2:(size(QN,2)-1)
                advecVar_curr=advecVar_dim1
                @views QN[:,ny] =advec1d( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt;
                    IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                # @views    QN[:,ny] =advec1d!( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt;pre=advecVar_curr);
            end
        end
        c_out = QN;
        # % for j direction
        Q0y = c_out;
        u1Dy = v;
        DXXy = ∆y;
        DEN0y = DEN1;
        # %DEN0y[2:end, :] = 1 - dtt ./ DXX[2:end, :] .* (d0[2:end, :] .* u[2:end, :] - d0[2:end, :] .* u[1:end - 1, :]);
        DD0y = d0;
        @views DEN1y[:, 1] = DEN0y[:, 1];
        # if false
        #     # if par == true
        #     Threads.@threads  for j=2:size(DEN1y,2)
        #         @views DEN1y[:, j] = DEN0y[:, j] .- dtt ./ ∆y[:, j] .* (d0[:, j] .* v[:, j] .- d0[:, j] .* v[:, j-1]);
        #     end
        # else
           @views  DEN1y[:, 2:end] = DEN0y[:, 2:end] .- dtt ./ ∆y[:, 2:end] .* (d0[:, 2:end] .* v[:, 2:end] .- d0[:, 2:end] .* v[:, 1:end - 1]);
        # end
        QNy=c_out;
        if par == true
            Threads.@threads  for nx = 2:(size(QNy,1)-1)
                advecVar_curr=reset!(advecVar_all_dim2[Threads.threadid()])
                # advecVar_curr=advecVar_all_dim2[Threads.threadid()]
                # advecVar_curr=advecVar_dim2
                @views QNy[nx,:] =advec1d( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt;
                # QNy[nx,:] =advec1d( view(Q0y,nx,:) , view(u1Dy,nx,:) , view(DEN0y,nx,:) , view(DEN1y,nx,:) , view(DXXy,nx,:)  , view(DD0y,nx,:)  , dtt;
                # IMXMN = pointer(advecVar_curr.IMXMN), FLUX = pointer(advecVar_curr.FLUX), VCMAX = pointer(advecVar_curr.VCMAX),VCMIN = pointer(advecVar_curr.VCMIN),   QN = pointer(advecVar_curr.QN));
                IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
            end
        else
            for nx = 2:(size(QNy,1)-1)
                advecVar_curr=advecVar_dim2
                    # IMXMN = advecVar_curr.IMXMN; FLUX = advecVar_curr.FLUX; VCMAX = advecVar_curr.VCMAX;VCMIN = advecVar_curr.VCMIN;   QN = advecVar_curr.QN;
                @views QNy[nx,:] =advec1d( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt;
                #     IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                # QNy[nx,:] =advec1d( view(Q0y,nx,:) , view(u1Dy,nx,:) , view(DEN0y,nx,:) , view(DEN1y,nx,:) , view(DXXy,nx,:)  , view(DD0y,nx,:) , dtt;
                #     IMXMN = IMXMN, FLUX = FLUX, VCMAX = VCMAX,VCMIN = VCMIN,   QN = QN);
                    IMXMN = advecVar_curr.IMXMN, FLUX = advecVar_curr.FLUX, VCMAX = advecVar_curr.VCMAX,VCMIN = advecVar_curr.VCMIN,   QN = advecVar_curr.QN);
                    
                # @views QNy[nx,:] =advec1d!( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt;    pre = advecVar_curr);
                
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

function adv_3d_basic(c_in::AbstractArray{T}, ρ0 , u, v,w, ∆x, ∆y,∆z, dt;par=false,num_tStep=1) where T
    # input C_in is the mass concentration
    # % the boundary will not change, however, it is necessary for calculating the inner area
    # %function adv_hori( myid, c_out, u, v, ∆x, ∆y, sx, ex, sy, ey ,nx,ny,dt)
    # accroding to the appendix in <Improving PM2.5 Forecasts in China Using an Initial Error Transport Model>
    c_out=copy(c_in);
    Q=c_in./ρ0
    u_next=u;u_pre=zeros(u);u_pre[2:end,:,:]=u[1:end-1,:,:]
    v_next=v;v_pre=zeros(v);v_pre[:,2:end,:]=v[:,1:end-1,:]
    z_next=w;z_pre=zeros(w);z_pre[:,:,2:end]=w[:,:,1:end-1]

    D0=ρ0
    D1=similar(D0);D2=similar(D0);D3=similar(D0);
    if 1 == length(∆x); ∆x = ∆x .+ zeros(T,size(c_out)); end
    if 1 == length(∆y); ∆y = ∆y .+ zeros(T,size(c_out)); end

    # divide dt into ddt
    if !(num_tStep >0) # num_tStep <= 0 or isnan
        dt0 = minimum(abs.(vcat(vec(∆x ./ u), vec(∆y ./ v))));
        num_tStep = ceil(dt / dt0);
    end
    dtt = dt / num_tStep;

    D1 = D0 .- (ρ0 .* u_next .- ρ0 .* u_pre) .* dtt ./ ∆x ;
    D2 = D1 .- (ρ0 .* v_next .- ρ0 .* v_pre) .* dtt ./ ∆y ;
    D3 = D2 .- (ρ0 .* z_next .- ρ0 .* z_pre) .* dtt ./ ∆z ;


    advStr_x=init_advec(c_in[:,1,1])
    advStr_y=init_advec(c_in[1,:,1])
    advStr_z=init_advec(c_in[1,1,:])
    for itt = 1:num_tStep
        # % for i direction        
        for ny in axes(Q,2), nz in axes(Q,3)
            @views Q[:,ny,nz] =advec1d( Q[:,ny,nz] , u[:,ny,nz] , D0[:,ny,nz] , D1[:,ny,nz] , ∆x[:,ny,nz] , ρ0[:,ny,nz] , dtt;            
                IMXMN = advStr_x.IMXMN, FLUX = advStr_x.FLUX, VCMAX = advStr_x.VCMAX,VCMIN = advStr_x.VCMIN,   QN = advStr_x.QN);
        end
        # % for j direction
        for nx in axes(Q,1), nz in axes(Q,3)
            @views Q[nx,:,nz] =advec1d( Q[nx,:,nz] , v[nx,:,nz] , D1[nx,:,nz] , D2[nx,:,nz] , ∆y[nx,:,nz] , ρ0[nx,:,nz] , dtt;            
                IMXMN = advStr_y.IMXMN, FLUX = advStr_y.FLUX, VCMAX = advStr_y.VCMAX,VCMIN = advStr_y.VCMIN,   QN = advStr_y.QN);
        end
        # for k direction
        for nx in axes(Q,1), ny in axes(Q,2)
            @views Q[nx,ny,:] =advec1d( Q[nx,ny,:] , w[nx,ny,:] , D2[nx,ny,:] , D3[nx,ny,:] , ∆z[nx,ny,:] , ρ0[nx,ny,:] , dtt;            
                IMXMN = advStr_z.IMXMN, FLUX = advStr_z.FLUX, VCMAX = advStr_z.VCMAX,VCMIN = advStr_z.VCMIN,   QN = advStr_z.QN);
        end
    end
    # c_out = Q.*ρ0;
    c_out = Q.*D3;
    return c_out
end


function advDiffAnHour(in, u, v, dx, dy;par=false,num_tStep=NaN)

    dt=60*60;
    transport2grid_max=1;
    if isnan(num_tStep)
        # numSub=ceil(max(u.*dt./dx...,v.*dt./dy...)/transport2grid_max);
        num_tStep=ceil(max(maximum(u.*dt./dx),maximum(v.*dt./dy))/transport2grid_max);
    end
    # dt_sub=dt/numSub;
    out=copy(in);
    out=adv_hori_wuhj( out, u, v, dx, dy,dt;par=par,num_tStep=num_tStep); 
    # for n_t=1:numSub
    #     out=adv_hori_wuhj( out, u, v, dx, dy,dt_sub;par=par,divideT=false); 
    #     # out=adv_hori_basic( out, u, v, dx, dy,dt_sub;par=par,divideT=false); 
    #     # out=diffusion_hori_Smagorinsky(out,u,v,dx,dy,dt_sub);
    # end
    return out
end
 # ================================ dump =================================
# function adv_hori_wuhj(c_in, u, v, ∆x, ∆y, dt)
#     # % the boundary will not change, however, it is necessary for calculating the inner area
#     # %function adv_hori( myid, c_out, u, v, ∆x, ∆y, sx, ex, sy, ey ,nx,ny,dt)
#     c_out=copy(c_in);

#     DEN1 = zeros(size(c_out)); DEN1y = zeros(size(c_out));

#     if 1 == length(∆x); ∆x = ∆x .+ zeros(size(c_out)); end
#     if 1 == length(∆y); ∆y = ∆y .+ zeros(size(c_out)); end

#     d0 = ones(size(c_out));

#     dt0 = min(abs.([∆x[:] ./ u[:]; ∆y[:] ./ v[:]])...);
#     istep = ceil(dt / dt0);
#     dtt = dt / istep;

#     for itt = 1:istep
#         # % for i direction
#         Q0 = c_out;
#         u1D = u;
#         DXX = ∆x;
#         DEN0 = d0;
#         DEN1[2:end, :] = 1 .- dtt ./ DXX[2:end, :] .* (d0[2:end, :] .* u[2:end, :] .- d0[2:end, :] .* u[1:end - 1, :]);
#         DEN1[1, :] .= 0;
#         DD0 = DEN0;
#         QN=c_out;
#         for ny = 2:(size(QN,2)-1)
#             QN[:,ny] =advec1d( Q0[:,ny] , u1D[:,ny] , DEN0[:,ny] , DEN1[:,ny] , DXX[:,ny] , DD0[:,ny] , dtt);
#         end
#         c_out = QN;
#         # % for j direction
#         Q0y = c_out;
#         u1Dy = v;
#         DXXy = ∆y;
#         DEN0y = DEN1;
#         # %DEN0y[2:end, :] = 1 - dtt ./ DXX[2:end, :] .* (d0[2:end, :] .* u[2:end, :] - d0[2:end, :] .* u[1:end - 1, :]);
#         DD0y = d0;
#         DEN1y[:, 1] = DEN0y[:, 1];
#         DEN1y[:, 2:end] = DEN0y[:, 2:end] - dtt ./ ∆y[:, 2:end] .* (d0[:, 2:end] .* v[:, 2:end] - d0[:, 2:end] .* v[:, 1:end - 1]);
#         QNy=c_out;
#         for nx = 2:(size(QNy,1)-1)
#             QNy[nx,:] =advec1d( Q0y[nx,:] , u1Dy[nx,:] , DEN0y[nx,:] , DEN1y[nx,:] , DXXy[nx,:] , DD0y[nx,:] , dtt);
#         end
#         c_out = QNy;
#         # % preserve the boundary
#         c_out[[1 end],:]=c_in[[1 end],:];
#         c_out[:,[1 end]]=c_in[:,[1 end]];
#         # %c_out=QNy.*DEN1y;
#     end
#     return c_out
# end
# function advDiffAnHour(in, u, v, dx, dy)

#     dt=60*60;
#     transport2grid_max=1;
#     numSub=ceil(max(u.*dt./dx...,v.*dt./dy...)/transport2grid_max);
#     dt_sub=dt/numSub;
#     out=copy(in);
#     for n_t=1:numSub
#         out=adv_hori_wuhj( out, u, v, dx, dy,dt_sub); 
#         # out=diffusion_hori_Smagorinsky(out,u,v,dx,dy,dt_sub);
#     end
#     return out
# end
# end