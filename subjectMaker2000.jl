# the subjectMaker2000 will parse the XDF object for the pupil labs stream and the PsyPy stream
# it will then return two DFs
# one called 'df' which contains all the pupil labs data
# another one called 'times' which contains all the time stamps correpsonding to markers in the experiment
# these two DFs can the be used by trialMaster30k to slice out single trials

using XDF, DataFrames

# column names for the XDF pupil labs object
global gaze_names = [
    "conf",
    "x_pos","y_pos",
    "gaze_point_3d_x", "gaze_point_3d_y", "gaze_point_3d_z",
    "eye_centerright_3d_x","eye_centerright_3d_y", "eye_centerright_3d_z",
    "eye_centerleft_3d_x", "eye_centerleft_3d_y", "eye_centerleft_3d_z",
    "gaze_normalright_x", "gaze_normalright_y", "gaze_normalright_z",
    "gaze_normalleft_x", "gaze_normalleft_y", "gaze_normalleft_z",
    "diameterright_2d", "diameterleft_2d",
    "diameterright_3d", "diameterleft_3d"
    ]


function subjectMaker2000(filename::String)

    # configuration for the answers in the experiment
    correct_T = [1,1,1,1]
    correct_C = [1, 0, 1, 0, 0]
    correct_I = [1,1,0,0,1,1,0,1,1,1,0,0,1,0,0,0]
    correct = [correct_T; correct_I; correct_C]

    #read the XDF file
    streams = read_xdf("../data/" * filename * ".xdf")

    pupil = []
    psypi = []
    times = []

    # get the right data from dict
    for i in 1:length(streams)
        if get(streams[i], "name",2) == "pupil_capture"
            pupil = streams[i]
            println("gawk")
        elseif get(streams[i], "name",2) == "insightful_stream_psycoPy"
            psypi = streams[i]
            println("bok")
        else
            println("beep")
        end
    end


    # add a header to the eye tracking data
    df = DataFrame(pupil["data"], gaze_names)

    #get time stamps for eye tracking data
    df.time = (pupil["time"])

    # build a dataframe with timestamps and the markers from psychopy
    df_psypi = DataFrame(psypi["data"],:auto)
    df_psypi.time = (psypi["time"])

    # function returning a tuple of marker, start time, RTA time, answer time and choice selected for a single trial
    function taskWindow(x)
        time  = df_psypi[df_psypi.x1 .== x, 2][1] # this [1] is just for the rare cases, where a marker is duplicate due to people repeating the training
        table = df_psypi[df_psypi.time .>= time, :][1:3,:]

        (table[1,1], table[1,2], table[2,2], table[3,2], table[3,1])
    end

    # now we will construct a data frame containing all the relevant trials
    # construct empty DataFrame
    times = DataFrame(Marker = [], Start = Float64[], RTA = Float64[], Answer = Float64[], Choice = Int[])


    for i in 1:4
        new = taskWindow(1000+i)
        push!(times, (new))
    end

    # add lines for control condition trials
    for i in 1:5
        new = taskWindow(2000+i)
        push!(times, (new))
    end

    # add lines for impasse condition trials

    for i in 1:16
        new = taskWindow(300+i)
        push!(times, (new))
    end

    sort!(times, order(:Marker))
    times.correct = correct
    times.right = (times.correct .== times.Choice)
    sort!(times, order(:Start))

    plot(times.right)
    println("A new subject as been born.")
    return df, times

end
