using GradDataM
export getFooInfo,read_aLayer,combineFoo,read_aVar,getTermInfo,getInfoFromFoo,combineTerm,getData_aPCorr
export getSMInfo,combineSM,getData_aCore
struct fooInfoStr
    sx_all::Vector
    ex_all::Vector
    sy_all::Vector
    ey_all::Vector
    lonlatSize::Tuple
    z_all::Vector # vector for z of layers
    ivar_all::Vector
    varName_all::Vector # vector for variable name, upper case
    nZ_all::Vector# vector for z of variables
end
struct termInfoStr
    sx_all::Vector
    ex_all::Vector
    sy_all::Vector
    ey_all::Vector
    lonlatSize::Tuple
    z_all::Vector
    p_all::Vector
    g_all::Vector
end

struct smInfoStr
    sx_all::Vector
    ex_all::Vector
    sy_all::Vector
    ey_all::Vector
    lonlatSize::Tuple
    z_all::Vector
    ism_all::Vector
    ivar_all::Vector
end

"""
# codes to view the foo file
f=open(fooNameFun(0)) 
info=Array{Int32}(undef,(12 + num_grid +1),2)
data=Array{Float32}(undef,(12 + num_grid +1),2)
read!(f,info);info=ntoh.(info)
info'
seek(f,0)
read!(f,data);data=ntoh.(data)
data'
close(f)

"""

function calIvarByZ(z_all)
    # return ivar_all starting from 1, increasing when z not increase
    ivar_all=similar(z_all)
    for n in eachindex(z_all)
        if n == 1 
            ivar_all[n] = 1;
        else
            if z_all[n] > z_all[n-1]
                ivar_all[n] = ivar_all[n-1];
            else
                ivar_all[n] = 1 + ivar_all[n-1];
            end
        end
    end
    return ivar_all
end

function getFooInfo(ctlFileName::AbstractString,fooNameFun)    
    grd=gradData(ctlFileName)
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all=getInfoFromFoo(fooNameFun)
    ivar_all=calIvarByZ(z_all)
    return fooInfoStr(sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,ivar_all,grd.varName_all,grd.nZ_all)
end
function getFooInfo(fooNameFun)    
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all=getInfoFromFoo(fooNameFun)
    ivar_all=calIvarByZ(z_all)
    return fooInfoStr(sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,ivar_all,[],[])
end
function getInfoFromFoo(fooNameFun)
    # get all files ----------------------------------
    file_all=[fooNameFun(0)];id_core=1;
    while isfile(fooNameFun(id_core))
        append!(file_all,[fooNameFun(id_core)])
        id_core+=1
    end
    # get all sx, ex, sy, ey ----------------------
    sx_all, ex_all, sy_all, ey_all, numGrid_all= (fill(0,length(file_all)) for i in 1:5)
    Threads.@threads for n_core  in eachindex(file_all)
        open(file_all[n_core]) do f
            info=Array{Int32}(undef,9); read!(f,info);info=ntoh.(info)
            z, sx, ex, sy, ey = info[3:7]; num_grid=(ex - sx + 1)*(ey - sy + 1)
            sx_all[n_core], ex_all[n_core], sy_all[n_core], ey_all[n_core], numGrid_all[n_core]=(sx,ex,sy,ey,num_grid)
        end
    end
    # get all z_all ------------------------
    z_all=[]    
    f= open(fooNameFun(1))
    while !eof(f)
        info=Array{Int32}(undef,9); read!(f,info);info=ntoh.(info)
        z, sx, ex, sy, ey = info[3:7];  num_grid=(ex - sx + 1)*(ey - sy + 1)
        # println(info[1:9]) # for testing
        skip(f,(1+num_grid)*4)
        append!(z_all,z)
    end
    # get others ---------------------
    lonlatSize=(max(ex_all...),max(ey_all...))
    return sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all
end

# fooNameFun(nCore)="/data1/elpt_2023_000305/Xiang/source10/sm/tmp/20220926/makd2.$(string(nCore,pad=3)).2022092619"
# f= open(fooNameFun(0))
# info=Array{Int32}(undef,40); read!(f,info);info=ntoh.(info)
# smInfo=getSMInfo(fooNameFun)
# write(id)myid,k,ism,idm,id,sx,ex,sy,ey

function getSMInfo(fooNameFun)

    # get all files ----------------------------------
    file_all=[fooNameFun(0)];id_core=1;
    while isfile(fooNameFun(id_core))
        append!(file_all,[fooNameFun(id_core)])
        id_core+=1
    end
    # get all sx, ex, sy, ey ----------------------
    sx_all, ex_all, sy_all, ey_all, numGrid_all= (fill(0,length(file_all)) for i in 1:5)
    Threads.@threads for n_core  in eachindex(file_all)
        open(file_all[n_core]) do f
            info=Array{Int32}(undef,12); read!(f,info);info=ntoh.(info)
            sx, ex, sy, ey = info[7:10]; num_grid=(ex - sx + 1)*(ey - sy + 1)
            sx_all[n_core], ex_all[n_core], sy_all[n_core], ey_all[n_core], numGrid_all[n_core]=(sx,ex,sy,ey,num_grid)
        end
    end
    # get all layer(z) process(p) and species(g) ------------------------
    z_all,ism_all,ivar_all=[],[],[] 
    f= open(fooNameFun(0))
    while !eof(f)
        info=Array{Int32}(undef,12); read!(f,info);info=ntoh.(info)
        myid, z, ism, ivar,ifile, sx, ex, sy, ey = info[2:10];  num_grid=(ex - sx + 1)*(ey - sy + 1)
        # println(info[1:10]) # for testing
        skip(f,(1+num_grid)*4)
        append!(z_all,z)
        append!(ism_all,ism)
        append!(ivar_all,ivar)
    end
    # get others ---------------------
    lonlatSize=(max(ex_all...),max(ey_all...))
    return smInfoStr( sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,ism_all,ivar_all)
end
function getTermInfo(fooNameFun)

    # get all files ----------------------------------
    file_all=[fooNameFun(0)];id_core=1;
    while isfile(fooNameFun(id_core))
        append!(file_all,[fooNameFun(id_core)])
        id_core+=1
    end
    # get all sx, ex, sy, ey ----------------------
    sx_all, ex_all, sy_all, ey_all, numGrid_all= (fill(0,length(file_all)) for i in 1:5)
    Threads.@threads for n_core  in eachindex(file_all)
        open(file_all[n_core]) do f
            info=Array{Int32}(undef,11); read!(f,info);info=ntoh.(info)
            sx, ex, sy, ey = info[6:9]; num_grid=(ex - sx + 1)*(ey - sy + 1)
            sx_all[n_core], ex_all[n_core], sy_all[n_core], ey_all[n_core], numGrid_all[n_core]=(sx,ex,sy,ey,num_grid)
        end
    end
    # get all layer(z) process(p) and species(g) ------------------------
    z_all,p_all,g_all=[],[],[] 
    f= open(fooNameFun(0))
    while !eof(f)
        info=Array{Int32}(undef,11); read!(f,info);info=ntoh.(info)
        myid, z, ip, ig, sx, ex, sy, ey = info[2:9];  num_grid=(ex - sx + 1)*(ey - sy + 1)
        # println(info[1:9]) # for testing
        skip(f,(1+num_grid)*4)
        append!(z_all,z)
        append!(p_all,ip)
        append!(g_all,ig)
    end
    # get others ---------------------
    lonlatSize=(max(ex_all...),max(ey_all...))
    return termInfoStr( sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,p_all,g_all)
end

# =============== for restart file or combine without extraction ========
# fooNameFun(nCore)="/data/tangxiao/wuhj/enFore_April/ideal/out/tmp/20180422/rfood3.$(string(nCore,pad=5)).2018042212"
# fileName_out="/data/tangxiao/wuhj/enFore_April/ideal/out/init/rtestd3.2018042212.grd"
# combineFoo(fooNameFun,fileName_out)
function combineFoo(fooNameFun,fileName_out)
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all=getInfoFromFoo(fooNameFun)

    dataAll=Array{Float32}(undef,(lonlatSize...,length(z_all)))
    #read --------------------------
    Threads.@threads for n_core in eachindex(sx_all)
        id_core=n_core-1
        file_curr=fooNameFun(id_core)
        open(file_curr) do f
            sx, ex, sy, ey =sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core]
            nx=(ex - sx + 1); ny=(ey - sy + 1);
            num_grid=nx*ny
            data_aCore=Array{Float32}(undef,(10+num_grid,length(z_all)))
            read!(f,data_aCore);data_aCore=ntoh.(data_aCore)
            dataAll[sx:ex,sy:ey,:]=reshape(data_aCore[10:end-1,:],(nx,ny,length(z_all)))
        end
    end
    #write ----------------------------
    # make folder if no exist
    outFolder=dirname(fileName_out)
    if ! isdir(outFolder);        mkpath(outFolder);     end
    open(fileName_out,"w") do f
        write(f,hton.(Float32.(dataAll)))
    end
end

# function combineFoo(fooNameFun,fileName_out,fooInfo::fooInfoStr)
#     dataAll=Array{Float32}(undef,(fooInfo.lonlatSize...,length(fooInfo.z_all)))
#     #read --------------------------
#     Threads.@threads for n_core in eachindex(fooInfo.sx_all)
#         id_core=n_core-1
#         file_curr=fooNameFun(id_core)
#         open(file_curr) do f
#             sx, ex, sy, ey =fooInfo.sx_all[n_core],fooInfo.ex_all[n_core],fooInfo.sy_all[n_core],fooInfo.ey_all[n_core]
#             nx=(ex - sx + 1); ny=(ey - sy + 1);
#             num_grid=nx*ny
#             data_aCore=Array{Float32}(undef,(10+num_grid,length(fooInfo.z_all)))
#             read!(f,data_aCore);data_aCore=ntoh.(data_aCore)
#             dataAll[sx:ex,sy:ey,:]=reshape(data_aCore[10:end-1,:],(nx,ny,length(fooInfo.z_all)))
#         end
#     end
#     #write ----------------------------
#     # make folder if no exist
#     outFolder=dirname(fileName_out)
#     if ! isdir(outFolder);        mkpath(outFolder);     end
#     open(fileName_out,"w") do f
#         write(f,hton.(Float32.(dataAll)))
#     end
# end

# fileFun_out(processName)="term/$processName/2022092500.grd"
function combineTerm(termNameFun,fileFun_out;processName_all=["proc$(string(n,pad=2))" for n in 1:24],isWriteCtl=true)
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,p_all,g_all=getTermInfo(termNameFun)
    varName_all=["H2SO4(g)","HNO3","HCL","NH3","NO","NO2","NO3","N2O5","HONO","HNO4","O3","O1D","O3P","OH","HO2","H2O2","CO","SO2","CH4","C2H6","CH3O2","ETHP","HCHO",
    "CH3OH","ANOL","CH3OOH","ETHOOH","ALD2","HCOOH","RCOOH","C2O3","PAN","PAR","AONE","MGLY","ETH","OLET","OLEI","TOL","XYL","CRES","TO2","CRO","OPEN","ONIT","ROOH",
    "RO2","ANO2","NAP","XO2","XPAR","ISOP","ISOPRD","ISOPP","ISOPN","ISOPO2","DMS","MSA","DMSO","DMSO2","CH3SO2H","CH3SCH2OO","CH3SO2","CH3SO3","CH3SO2OO","CH3SO2CH2OO",
    "SULFH","TERP","SV1","SV2","SV3","SV4","SV5","SV6","PM25","PM10","BC","OC","H+(AQ)","NA+(AQ)","NH4+(AQ)","CL-(AQ)","SO4--(AQ)","HSO4-(AQ)","NO3-(AQ)","NACL(S)",
    "NA2SO4(S)","NANO3(S)","NH42SO4(S)","NH4NO3(S)","NH4CL(S)","H2SO4(AQ)","NH4HSO4(S)","NAHSO4(S)","(NH4)4H(SO4)2(S)","SOA1","SOA2","SOA3","SOA4","SOA5","SOA6","AH2O"]

    dataAll=Array{Float32}(undef,(lonlatSize...,length(p_all)))
    #read --------------------------
    Threads.@threads for n_core in eachindex(sx_all)
        id_core=n_core-1
        file_curr=termNameFun(id_core)
        open(file_curr) do f
            sx, ex, sy, ey =sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core]
            nx=(ex - sx + 1); ny=(ey - sy + 1);
            num_grid=nx*ny
            data_aCore=Array{Float32}(undef,(12+num_grid,length(p_all)))
            read!(f,data_aCore);data_aCore=ntoh.(data_aCore)
            dataAll[sx:ex,sy:ey,:]=reshape(data_aCore[12:end-1,:],(nx,ny,length(p_all)))
        end
    end
    #write ----------------------------
    for n_proc in unique(p_all)
        fileName_out=fileFun_out(processName_all[n_proc])
        outFolder=dirname(fileName_out)
        if ! isdir(outFolder);        mkpath(outFolder);     end
        q_currProc= n_proc .== p_all
        open(fileName_out,"w") do f
            write(f,hton.(Float32.(dataAll[:,:,q_currProc])))
        end
        if isWriteCtl
            ctlName=splitext(fileName_out)[1]*".ctl"
            varNameS=varName_all[unique(g_all)]
            zS=fill(maximum(z_all),length(varNameS))
            writeCtl(ctlName,lonlatSize,varNameS,zS)
        end
    end
end


function getData_aPCorr(file_curr,inds_currProc,sx, ex, sy, ey,readMethod)
    local dataOut
    nx=(ex - sx + 1); ny=(ey - sy + 1);
    num_grid=nx*ny

    open(file_curr) do f
        if readMethod == "sequence"
            # method 2, read all the layers and then select
            data_aCore=Array{Float32}(undef,(12+num_grid,length(p_all)))
            read!(f,data_aCore);#data_aCore=ntoh.(data_aCore)
            # data_aCore=reshape(data_aCore[12:end-1,inds_currProc],(nx,ny,length(inds_currProc)))
            @views dataOut=data_aCore[12:end-1,inds_currProc]
        else
            # method 1, read layer by layer
            data_aCore=Array{Float32}(undef,(12+num_grid,length(inds_currProc)))
            # data_aLayer=Array{Float32}(undef,12+num_grid)
            for (n_slice,ind_slice) in enumerate(inds_currProc)
                seek(f,4*(12+num_grid)*(ind_slice-1))
                read!(f,view(data_aCore,:,n_slice));
            end
            #data_aCore=ntoh.(data_aCore)
            # data_aCore=reshape(data_aCore[12:end-1,:],(nx,ny,length(inds_currProc)))
            @views dataOut=data_aCore[12:end-1,:]
        end
    end
    return reshape(dataOut,nx,ny,length(inds_currProc))
end

function getData_aCore(file_curr,iSlice_all,sx, ex, sy, ey,readMethod;infoLen=11)
    local dataOut
    nx=(ex - sx + 1); ny=(ey - sy + 1);
    num_grid=nx*ny

    open(file_curr,"r") do f
        if readMethod == "sequence"
            # method 2, read all the layers and then select
            data_aCore=Array{Float32}(undef,(infoLen+num_grid+1,maximum(iSlice_all)))
            read!(f,data_aCore);
            @views dataOut=data_aCore[(infoLen+1):end-1,iSlice_all]
        else
            # method 1, read layer by layer
            data_aCore=Array{Float32}(undef,(infoLen+num_grid+1,length(iSlice_all)))
            for (n_slice,ind_slice) in enumerate(iSlice_all)
                seek(f,4*(infoLen+num_grid+1)*(ind_slice-1))
                read!(f,view(data_aCore,:,n_slice));
            end
            @views dataOut=data_aCore[(infoLen+1):end-1,:]
        end
    end
    # println(dataOut);
    return reshape(dataOut,nx,ny,length(iSlice_all))
end


# fooNameFun(nCore)="/data1/elpt_2023_000305/Xiang/source10/sm/tmp/20220926/makd2.$(string(nCore,pad=3)).2022092619"
# fileFun_out(nZ)="/data1/elpt_2023_000305/Xiang/source10/sm/data/makd2.z$(nZ).2022092619.grd"
# smInfo=getSMInfo(fooNameFun)
# # write(id)myid,k,ism,idm,id,sx,ex,sy,ey
# combineSM(fooNameFun,fileFun_out,smInfo;disp=true,readMethod="slices")


# fileFun_out(processName)="term/$processName/2022092500.grd"
function combineSM(fooNameFun,fileFun_out,smInfoStr;outZ_all=unique(smInfoStr.z_all),outIvar_all=unique(smInfoStr.ivar_all),disp=false,readByVar=true,readMethod="slices")
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,ivar_all,ism_all=smInfoStr.sx_all,smInfoStr.ex_all,smInfoStr.sy_all,smInfoStr.ey_all,smInfoStr.lonlatSize,smInfoStr.z_all,smInfoStr.ivar_all,smInfoStr.ism_all
    fooName_all=collect(fooNameFun(n_core-1) for n_core in eachindex(sx_all))
    for z_curr in unique(outZ_all)
        fileName_out=fileFun_out(z_curr)
        if isfile(fileName_out);rm(fileName_out);end # remove the file if not exist
        if readByVar
            for iVar in outIvar_all
                if disp;println("Extracting z_curr = $(z_curr), iVar = $(iVar)");end
                iSlice_all= findall( (z_all .== z_curr) .& (ivar_all .== iVar) )
                data_aVar=Array{Float32}(undef,(lonlatSize...,length(iSlice_all))) 
                # data_allCore=Array{Any}(undef,length(sx_all))    
                #read --------------------------        
                Threads.@threads for n_core in eachindex(sx_all)
                # for n_core in eachindex(sx_all)
                #     if disp;print(n_core," ");end
                    data_aVar[sx_all[n_core]:ex_all[n_core],sy_all[n_core]:ey_all[n_core],:]=getData_aCore(fooName_all[n_core],iSlice_all,sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core],readMethod;infoLen=12)
                end 
                # write the file in append mode ------------------------
                open(fileName_out,"a") do f
                    write(f,data_aVar)
                end
            end
        else
            iSlice_all= findall( (z_all .== z_curr) .& (in.(ivar_all,[outIvar_all])) )
            if disp;println("Extracting z_curr = $(z_curr)");end
            data_aZ=Array{Float32}(undef,(lonlatSize...,length(iSlice_all))) 
            #read --------------------------        
            Threads.@threads for n_core in eachindex(sx_all)
                # data_allCore[n_core]=getData_aCore(fooName_all[n_core],iSlice_all,sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core],readMethod;infoLen=12)
                data_aZ[sx_all[n_core]:ex_all[n_core],sy_all[n_core]:ey_all[n_core],:]=getData_aCore(fooName_all[n_core],iSlice_all,sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core],readMethod;infoLen=12)
            end 
            # write the file in append mode ------------------------
            if disp;println("Writing $(fileName_out)");end
            open(fileName_out,"w") do f
                write(f,data_aZ)
            end
        end
    end
end

# fileFun_out(processName)="term/$processName/2022092500.grd"
function combineFoo(fooNameFun,fileName_out,fooInfoStr;outZ_all=unique(fooInfoStr.z_all),outIvar_all=unique(fooInfoStr.ivar_all),readByVar=false,disp=false,readMethod="slices")
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,ivar_all=
        fooInfoStr.sx_all,fooInfoStr.ex_all,fooInfoStr.sy_all,fooInfoStr.ey_all,fooInfoStr.lonlatSize,fooInfoStr.z_all,fooInfoStr.ivar_all
    fooName_all=collect(fooNameFun(n_core-1) for n_core in eachindex(sx_all))
    if isfile(fileName_out);rm(fileName_out);end # remove the file if not exist
    if readByVar
        for iVar in outIvar_all
            if disp;println("Extracting iVar = $(iVar)");end
            iSlice_all= findall( (z_all .== z_curr) .& (ivar_all .== iVar) )
            data_aVar=Array{Float32}(undef,(lonlatSize...,length(iSlice_all))) 
            #read --------------------------        
            Threads.@threads for n_core in eachindex(sx_all)
                data_aVar[sx_all[n_core]:ex_all[n_core],sy_all[n_core]:ey_all[n_core],:]=getData_aCore(fooName_all[n_core],iSlice_all,sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core],readMethod;infoLen=9)
            end 
            # write the file in append mode ------------------------
            open(fileName_out,"a") do f
                write(f,data_aVar)
            end
        end
    else
        iSlice_all= findall( in.(z_all,[outZ_all])  .&&  in.(ivar_all,[outIvar_all]))
        data_aZ=Array{Float32}(undef,(lonlatSize...,length(iSlice_all))) 
        #read --------------------------        
        Threads.@threads for n_core in eachindex(sx_all)
            data_aZ[sx_all[n_core]:ex_all[n_core],sy_all[n_core]:ey_all[n_core],:]=getData_aCore(fooName_all[n_core],iSlice_all,sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core],readMethod;infoLen=9)
        end 
        # write
        if disp;println("Writing $(fileName_out)");end
        open(fileName_out,"w") do f
            write(f,data_aZ)
        end
    end
end

# fileFun_out(processName)="term/$processName/2022092500.grd"
function combineTerm(termNameFun,fileFun_out,termInfoStr;processName_all=["proc$(string(n,pad=2))" for n in 1:24],isWriteCtl=true,disp=true,readMethod="slices",
    nproc2ext=unique(termInfoStr.p_all))
    sx_all, ex_all, sy_all, ey_all,lonlatSize,z_all,p_all,g_all=termInfoStr.sx_all,termInfoStr.ex_all,termInfoStr.sy_all,termInfoStr.ey_all,termInfoStr.lonlatSize,termInfoStr.z_all,termInfoStr.p_all,termInfoStr.g_all
    varName_all=["H2SO4(g)","HNO3","HCL","NH3","NO","NO2","NO3","N2O5","HONO","HNO4","O3","O1D","O3P","OH","HO2","H2O2","CO","SO2","CH4","C2H6","CH3O2","ETHP","HCHO",
    "CH3OH","ANOL","CH3OOH","ETHOOH","ALD2","HCOOH","RCOOH","C2O3","PAN","PAR","AONE","MGLY","ETH","OLET","OLEI","TOL","XYL","CRES","TO2","CRO","OPEN","ONIT","ROOH",
    "RO2","ANO2","NAP","XO2","XPAR","ISOP","ISOPRD","ISOPP","ISOPN","ISOPO2","DMS","MSA","DMSO","DMSO2","CH3SO2H","CH3SCH2OO","CH3SO2","CH3SO3","CH3SO2OO","CH3SO2CH2OO",
    "SULFH","TERP","SV1","SV2","SV3","SV4","SV5","SV6","PM25","PM10","BC","OC","H+(AQ)","NA+(AQ)","NH4+(AQ)","CL-(AQ)","SO4--(AQ)","HSO4-(AQ)","NO3-(AQ)","NACL(S)",
    "NA2SO4(S)","NANO3(S)","NH42SO4(S)","NH4NO3(S)","NH4CL(S)","H2SO4(AQ)","NH4HSO4(S)","NAHSO4(S)","(NH4)4H(SO4)2(S)","SOA1","SOA2","SOA3","SOA4","SOA5","SOA6","AH2O"]

    #write ----------------------------
    for n_proc in nproc2ext
        fileName_out=fileFun_out(processName_all[n_proc])
        if disp; println("processing ",fileName_out);end
        outFolder=dirname(fileName_out)
        if ! isdir(outFolder);        mkpath(outFolder);     end
        inds_currProc= findall(n_proc .== p_all)
        data_aP=Array{Float32}(undef,(lonlatSize...,length(inds_currProc)))    
        data_allCore=Array{Any}(undef,length(sx_all))    
        #read --------------------------        
        Threads.@threads for n_core in eachindex(sx_all)
        # for n_core in eachindex(sx_all)
            data_allCore[n_core]=getData_aPCorr(termNameFun(n_core-1),inds_currProc,sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core],readMethod)
            # id_core=n_core-1
            # file_curr=termNameFun(id_core)
            # open(file_curr) do f
            #     sx, ex, sy, ey =sx_all[n_core],ex_all[n_core],sy_all[n_core],ey_all[n_core]
            #     nx=(ex - sx + 1); ny=(ey - sy + 1);
            #     num_grid=nx*ny

            #     if readMethod == "sequence"
            #         # method 2, read all the layers and then select
            #         data_aCore=Array{Float32}(undef,(12+num_grid,length(p_all)))
            #         read!(f,data_aCore);
            #         data_aP[sx:ex,sy:ey,:]=reshape(data_aCore[12:end-1,inds_currProc],(nx,ny,length(inds_currProc)))
            #     else
            #         # method 1, read layer by layer
            #         n_out=1
            #         data_aLayer=Array{Float32}(undef,12+num_grid)
            #         for p in p_all
            #             if p == n_proc
            #                 read!(f,data_aLayer);
            #                 data_aP[sx:ex,sy:ey,n_out]=reshape(data_aLayer[12:end-1],(nx,ny))
            #                 n_out+=1
            #             else # skip if not current process
            #                 skip(f,4*(12+num_grid))
            #             end
            #         end
            #     end
            # end
        end
        # combine ------------------------
        for n_core in eachindex(sx_all)
            data_aP[sx_all[n_core]:ex_all[n_core],sy_all[n_core]:ey_all[n_core],:]=data_allCore[n_core]
        end
        # write ----------------------------
        open(fileName_out,"w") do f
            # write(f,hton.(Float32.(data_aP)))
            write(f,data_aP)
        end
        if isWriteCtl
            ctlName=splitext(fileName_out)[1]*".ctl"
            varNameS=varName_all[unique(g_all)]
            zS=fill(maximum(z_all[inds_currProc]),length(varNameS))
            writeCtl(ctlName,lonlatSize,varNameS,zS)
        end
    end
end

function read_aVar(fooInfo::fooInfoStr,foodNameFun,varName)
    num_z=fooInfo.nZ_all[findfirst(x->x == uppercase(varName),fooInfo.varName_all)]
    dataOut=fill(NaN32,(fooInfo.lonlatSize...,num_z))
    Threads.@threads  for n_core=eachindex(fooInfo.sx_all)
        id_core=n_core-1
        file_curr=foodNameFun(id_core)
        sx, ex, sy, ey =fooInfo.sx_all[n_core],fooInfo.ex_all[n_core],fooInfo.sy_all[n_core],fooInfo.ey_all[n_core]
        open(file_curr) do f
            for z_curr=1:num_z
                seek(f,getDataLoc(fooInfo,n_core,varName,z_curr))
                data_aCore=Array{Float32}(undef,(ex-sx+1,ey-sy+1))
                read!(f,data_aCore);data_aCore=ntoh.(data_aCore)
                dataOut[sx:ex,sy:ey,z_curr]=data_aCore
            end
        end
    end
    return dataOut
end

function read_aLayer(fooInfo::fooInfoStr,foodNameFun,varName,nZ)
    dataOut=fill(NaN32,fooInfo.lonlatSize)
    @sync for n_core=eachindex(fooInfo.sx_all)
        id_core=n_core-1
        file_curr=foodNameFun(id_core)
        @async open(file_curr) do f
            skip(f,getDataLoc(fooInfo,n_core,varName,nZ))
            sx, ex, sy, ey =fooInfo.sx_all[n_core],fooInfo.ex_all[n_core],fooInfo.sy_all[n_core],fooInfo.ey_all[n_core]
            data_aCore=Array{Float32}(undef,(ex-sx+1,ey-sy+1))
            read!(f,data_aCore);data_aCore=ntoh.(data_aCore)
            dataOut[sx:ex,sy:ey]=data_aCore
        end
    end
    return dataOut
end

# function read_aLayer(fooInfo::fooInfoStr,foodNameFun,varName,nZ)
#     dataOut=fill(NaN32,fooInfo.lonlatSize)
#     for n_core=eachindex(fooInfo.sx_all)
#         id_core=n_core-1
#         file_curr=foodNameFun(id_core)
#         open(file_curr) do f
#             skip(f,getDataLoc(fooInfo,n_core,varName,nZ))
#             sx, ex, sy, ey =fooInfo.sx_all[n_core],fooInfo.ex_all[n_core],fooInfo.sy_all[n_core],fooInfo.ey_all[n_core]
#             data_aCore=Array{Float32}(undef,(ex-sx+1,ey-sy+1))
#             read!(f,data_aCore);data_aCore=ntoh.(data_aCore)
#             dataOut[sx:ex,sy:ey]=data_aCore
#         end
#     end
#     return dataOut
# end
function getDataLoc(fooInfo::fooInfoStr,nCorr,varName,nZ)
    infoSize_start=9
    return getLayerLoc(fooInfo,nCorr,varName,nZ)+infoSize_start*4
end
function getLayerLoc(fooInfo::fooInfoStr,nCorr,varName,nZ)
    ind_var=findfirst(x->x==uppercase(varName),fooInfo.varName_all)
    if isnothing(ind_var);@error varName * " not found in fooInfo";end
    num_layer2skip=sum(fooInfo.nZ_all[1:ind_var]) - fooInfo.nZ_all[ind_var] + nZ - 1    
    infoSize_start=9;infoSize_end=1;
    sx, ex, sy, ey =fooInfo.sx_all[nCorr],fooInfo.ex_all[nCorr],fooInfo.sy_all[nCorr],fooInfo.ey_all[nCorr]
    return  ((ex-sx+1)*(ey-sy+1)+infoSize_start+infoSize_end)*num_layer2skip*4
end
# ctlFileName="testd1.ctl"
# foodNameFun(id_core)=raw"tmp/food2." * string(id_core,pad=5) * ".2018022013"

# fooInfo=getFooInfo(ctlFileName,foodNameFun)
# # include("GradDataM.jl")
# varName="CO2"
# nZ=1
# data=read_aLayer(fooInfo,foodNameFun,"co2",1)
# using PyPlot
# pcolor(data')
# show()

# file_curr=foodNameFun(0)
# f=open(file_curr)
# seek(f,0)
# info=Array{Int32}(undef,9); read!(f,info);info=ntoh.(info)
# z, sx, ex, sy, ey = info[3:7]; num_grid=(ex - sx + 1)*(ey - sy + 1)
# skip(f,(1+num_grid)*4)
# info=Array{Float32}(undef,15)
