
using DataFrames
export carefulImportPoll
"""
fileName="/mnt/data01/obs_jinan/201301/20130102120000.txt"
df_data=carefulImportPoll(fileName;disp=true)
"""
function carefulImportPoll(filename::AbstractString;disp=false,writeUnrecognizedLine=false)
    # staID is a cell containing strings; Data contains 'longitude','latitude','PM25','PM10','CO','NO2','SO2','O3'
    # discard O38
    # 2021-12-31 updated by wuhj to support string staID
    # return if file no valid ========================
    if !isfile(filename)
        if disp;println("file $filename does not exist");end
        return nothing
    end
    # read all lines and initialize outputs ================
    allLines = readlines(filename)
    num_lines = length(allLines)
    staID_all = fill("", num_lines)
    Data_all = fill(NaN, (num_lines, 8))
    # parse data line by line ========================
    n_line = 1
    while n_line <= num_lines # instead of 'if', 'while' is used in here to support combining multiple lines
        linStr_ori = allLines[n_line]
        lineCell, lineData = line2numbers(linStr_ori)
        staID_curr, otherData,n_line = checkALine(lineCell, lineData, linStr_ori,n_line,num_lines,allLines)
        # save to Buffer
        if isempty(otherData)
            @warn "unrecognized line : $linStr_ori";
            # write the warning message to file with append mode to log.insertDF
            if writeUnrecognizedLine
                open("log.unrecognizedLine", "a") do io
                    write(io, "$filename: $linStr_ori\n")
                end
            end

            # if disp;println("unrecognized line : $linStr_ori");end
        else
            Data_all[n_line, :] = otherData
            staID_all[n_line] = staID_curr
        end
        n_line += 1
    end

    # return if not valid data
    if isempty(Data_all)
        return nothing
    end
    # exclue invalid data ========================
    Data_all[:, 1] = round.(Data_all[:, 1], digits = 3)
    Data_all[:, 2] = round.(Data_all[:, 2], digits = 4)
    q_invalid = hcat(Data_all[:, 1:2] .== -999, Data_all[:, 3:end] .< 0)
    Data_all[q_invalid] .= NaN    
    Data_all=replace(Data_all,NaN=>missing)
    # exclude invalid rows
    q_invalid= ismissing.(Data_all)
    # q_valid=  .!(q_invalid[:,1] .| q_invalid[:,2] .| vec(sum(.!q_invalid[:,3:end], dims=2)) .== 0)
    q_valid=  .!isempty(staID_all) .& (vec(sum(.!q_invalid[:,3:end], dims=2)) .> 0)
    # q_valid=fill(true,size(q_valid))
    if disp;println("removing $(sum(.!q_valid)) rows without valid observations");end

    # construct data frames for output    
    df_data=hcat(DataFrame(StaID=staID_all[q_valid]),DataFrame(Data_all[q_valid,:],Symbol.(["lon","lat","PM25","PM10","CO","NO2","SO2","O3"])))
    # df_data.staID = staID_all[q_valid]
    
    # remove duplicate
    staID_unique=unique(df_data.StaID)
    if length(staID_unique) < length(df_data.StaID)
        @warn "removing $(length(df_data.StaID) - length(staID_unique)) duplicate staID"
        df_data=df_data[indexin(staID_unique,df_data.StaID),:]
    end
    return df_data
    
    # return staID_all, q_valid,Data_all
end
function line2numbers(lineStr)
    C = split(lineStr)
    lineData = tryparse.(Float64, C)
    return C, lineData
end

function checkALine(C, lineData, linStr_ori,n_line,num_lines,allLines)
    # check the line, return empty if fail the check
    if length(lineData) == 12 && lineData[11] == -999
        # case 1 data in the right format, C{10} maybe -999 or o38  *********************
        otherData = lineData[2:9]
        staID_curr = C[1]
    else
        if length(lineData) == 11 && C[2] == "null"
            # case 2 without longitude and latitude  *********************
            staID_curr = C[1]
            otherData = vcat(NaN, NaN, lineData[3:8])
        elseif length(C) < 12
            # case 3 no enough data   *********************
            if n_line + 1 > num_lines # end of file
                staID_curr = ""
                otherData = [] # drop it
            else
                # read the next line, combine them together and check the combined line
                linStr_oriNext = allLines[n_line + 1]
                Cnext, lineDataNext = line2numbers(linStr_oriNext)
                staID_curr, otherData , n_line= checkALine(vcat(C, Cnext), vcat(lineData, lineDataNext), string(linStr_ori, linStr_oriNext),n_line,num_lines,allLines)
                if !isempty(otherData)
                    n_line += 1
                end
            end
        else # too many data
            # don't contain the right number of elements in a line and no recognized
            staID_curr = ""
            otherData = [] # drop it
        end
    end
    return staID_curr, otherData, n_line
end