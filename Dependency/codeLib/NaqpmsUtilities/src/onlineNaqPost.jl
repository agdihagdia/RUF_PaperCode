using Statistics,TimeZones,Dates,GradDataM#,MeshDataM
export onlineNaqPost_var,onlineNaqPost_grd,naqPost_grd,waitTillReady
# ===============sample usage:   run this post in shell by :
# julia -t 10 -e '''include("test_onlineNaqPost_conlumn.jl"); onlineNaqPost_column("/data/tangxiao/wuhj/enFore_April/member04/out/tmp","testd1.ctl","meshDATA2018052618_m04.nc","2018-05-26T00","2018-05-26T23") '''

# function onlineNaqPost_var(fooNameFun,ctlFileName,Times,var_all,fileNameOut;disp=true,seconds2wait=60*10)
#     # ============= configuration =============
#     tStart=Times[1]
#     leads_all=[Hour(t-tStart).value for t in Times]
    
#     # get fooInfo
#     fooInfo=getFooInfo(ctlFileName,nCore->fooNameFun(nCore,tStart))
#     num_core=length(fooInfo.sx_all)
#     ind_all=indexin(uppercase.(var_all),uppercase.(fooInfo.varName_all))
#     nZs=[1:fooInfo.nZ_all[ind_all[1]];]
#     # ============== read data ==============
#     meshObj=meshDATA(fileNameOut,fooInfo.lonlatSize,nZs,leads_all,var_all,tStart)
#     for (n_t,t_curr) in enumerate(Times)
#         if disp;println("-----------reading ",t_curr," ");end
#         # ------------- wait till file ready -------
#         waitTillReady(()->isFileReady(t->fooNameFun(num_core-1,t),t_curr),seconds2wait)
#         sleep(0.1)
#         # --------- read and write
#         data_curr=fill(NaN,(fooInfo.lonlatSize...,length(nZs),1,length(var_all),1))
#         for n_var  in eachindex(var_all)
#             data_curr[:,:,:,:,n_var,:]=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),var_all[n_var])
#         end
#         write(meshObj,data_curr;leadTimes=leads_all[n_t],varNames=var_all)
#     end
# end

"""

sample usage:
# cd to the folder of naqpms
tStart=DateTime(2023,4,10,12)
Times=tStart .+ Hour.(0:244)
# prepare folder and ctlFile
outFolder="out/data"
if ! isdir(outFolder);    mkpath(outFolder);end
cp("testd1.ctl",outFolder)
fooNameFun=
ctlName_out,Times

"""
function onlineNaqPost_grd(fooNameFun,ctlName_out,Times;disp=true,seconds2wait=60*10,deleteFoo=false,ctlName_foo=ctlName_out)
    # ============= configuration =============
    tStart=Times[1]
    # get fooInfo --------------------
    # start post process after the second hour exist, because it is troublesome to check when the data of the first hour is finished
    waitTillReady(()->isfile(fooNameFun(0,Times[2])),seconds2wait) # wait for the second hour of food to be ready
    fooInfo=getFooInfo(ctlName_foo,nCore->fooNameFun(nCore,tStart))
    num_core=length(fooInfo.sx_all)
    
    # prepare output folder ----------------------
    grdOut=setTime!(gradData(ctlName_out),tStart)
    file_tmp=grdOut.datFileName
    outfolder=dirname(file_tmp);    if ! isdir(outfolder);        mkpath(outfolder);     end
    # get the varName to extract -------------
    var_all=grdOut.varName_all
    # ============== read data ==============    
    for t_curr in Times
        setTime!(grdOut,t_curr)
        if disp;println("-----------reading ",t_curr," ");end
        # ------------- wait till file ready -------
        waitTillReady(()->isFileReady(fooNameFun,num_core,t_curr),seconds2wait)
        sleep(0.1) # for robust, wait for file ready
        # --------- read and write
        for var_curr in var_all
            if uppercase.(var_curr) == "CO2_COL"
                co2_curr=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),"co2") 
                press_curr=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),"pres")
                data_aVar=calxco2_3d(co2_curr,press_curr)
            else
                # read directly
                data_aVar=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),var_curr)
            end
            write(grdOut,data_aVar,var_curr)
        end
        # delete the file after Post extraction
        if deleteFoo
            for nCore in 0:num_core
                rm(fooNameFun(nCore,t_curr), force=true)
            end
        end
    end
end

function naqPost_grd(fooNameFun,ctlName_out,Times;disp=true,deleteFoo=false,ctlName_foo=ctlName_out)
    # ============= configuration =============
    tStart=Times[1]
    # get fooInfo --------------------
    # start post process after the second hour exist, because it is troublesome to check when the data of the first hour is finished
    fooInfo=getFooInfo(ctlName_foo,nCore->fooNameFun(nCore,tStart))
    num_core=length(fooInfo.sx_all)
    
    # prepare output folder ----------------------
    grdOut=setTime!(gradData(ctlName_out),tStart)
    file_tmp=grdOut.datFileName
    outfolder=dirname(file_tmp);    if ! isdir(outfolder);        mkpath(outfolder);     end
    # get the varName to extract -------------
    var_all=grdOut.varName_all
    # ============== read data ==============    
    Threads.@threads for t_curr in Times
        grd_curr=copy(grdOut)
        setTime!(grd_curr,t_curr)
        if disp;println("-----------reading ",t_curr," ");end
        # --------- read and write
        for var_curr in var_all
            if uppercase.(var_curr) == "CO2_COL"
                co2_curr=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),"co2") 
                press_curr=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),"pres")
                data_aVar=calxco2_3d(co2_curr,press_curr)
            else
                # read directly
                data_aVar=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),var_curr)
            end
            write(grd_curr,data_aVar,var_curr)
        end
        # delete the file after Post extraction
        if deleteFoo
            for nCore in 0:num_core
                rm(fooNameFun(nCore,t_curr), force=true)
            end
        end
    end
end

# # ================= for direct variables only ===================
# function onlineNaqPost_grd(fooNameFun,ctl_foo,Times,var_all,fileFun_grd;disp=true,seconds2wait=60*10)
#     # ============= configuration =============
#     tStart=Times[1]

#     # get fooInfo --------------------
#     # start post process after the second hour exist, because it is troublesome to check when the data of the first hour is finished
#     waitTillReady(()->isfile(fooNameFun(0,Times[2])),seconds2wait)
#     fooInfo=getFooInfo(ctl_foo,nCore->fooNameFun(nCore,tStart))
#     num_core=length(fooInfo.sx_all)
#     ind_all=indexin(uppercase.(var_all),uppercase.(fooInfo.varName_all))
#     nZs=fooInfo.nZ_all[ind_all];

#     # write ctl for the first hour ----------------------
#     file_tmp=fileFun_grd(tStart)
#     outfolder=dirname(file_tmp);    if ! isdir(outfolder);        mkpath(outfolder);     end
#     ctlName_out=splitext(file_tmp)[1] * ".ctl"
#     writeCtl(ctlName_out,fooInfo.lonlatSize,var_all,nZs)
#     grdOut=gradData(ctlName_out)
#     # ============== read data ==============
#     for n_t =1:length(Times)
#         t_curr=Times[n_t]
#         grdOut.datFileName=fileFun_grd(t_curr)
#         if disp;println("-----------reading ",t_curr," ");end
#         # ------------- wait till file ready -------
#         waitTillReady(()->isFileReady(fooNameFun,num_core,t_curr),seconds2wait)
#         sleep(0.1)
#         # --------- read and write
#         for n_var =1:length(var_all)
#             var_curr=var_all[n_var]
#             data_aVar=read_aVar(fooInfo,nCore->fooNameFun(nCore,t_curr),var_curr)
#             write(grdOut,data_aVar,var_curr)
#         end
#     end
# end

function waitTillReady(isReadyFun,seconds2wait)
    t_start=now()
    while (now()-t_start) < Second(seconds2wait)
        if isReadyFun()
            break
        else
            sleep(1)
        end
    end
    if !isReadyFun()
        @error "foo files not read, break........"
    end

end
function isFileReady(nameFun,num_core,t2check)
    # the file is ready if the next hour file exist or the size equal to the previous hour
    if isfile(nameFun(0,t2check)) && isfile(nameFun(0,t2check+Hour(1)))
        return true # return true if next hour file exist
    else
        num_coreReady=0
        for core = 0: (num_core-1)
            if filesize(nameFun(core,t2check)) == filesize(nameFun(core,t2check-Hour(1)))
                num_coreReady +=1
            end
        end
        if num_coreReady ==   num_core
            return true # return true if file size equal to previous hour
        else
            return false
        end
    end
end
function isFileReady(nameFun,t2check)
    # the file is ready if the next hour file exist or the size equal to the previous hour
    if isfile(nameFun(t2check)) && isfile(nameFun(t2check+Hour(1)))
        return true # return true if next hour file exist
    elseif filesize(nameFun(t2check)) == filesize(nameFun(t2check-Hour(1)))
        # poltential bug there, some core might not finished----------------------
        return true # return true if file size equal to previous hour
    else
        return false
    end
end
# # ============== test ==============
# tmpFolder="/data/tangxiao/wuhj/enFore_April/member04/out/tmp"
# ctlFileName,fileNameOut,tStr_start,tStr_end = ("testd1.ctl","meshDATA2018052618_m04.nc","2018-05-26T00","2018-05-26T23")
# onlineNaqPost_column(tmpFolder,ctlFileName,fileNameOut,tStr_start,tStr_end)

# "/data/tangxiao/wuhj/enFore_April/member04/out/tmp/20180526/food2.00018.2018052619"
