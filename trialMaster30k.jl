# the trialMaster30k will use the 'df' and 'time' DFs from subjectMaker2000 to slice single trials out of the pupil labs stream
# it will then low pass the relevant signals
# it will then compute a bunch of candidate variables for subsequent analysis
# it returns a single DF called 'df_trial', which contains all the candidate variables

using Distances, DSP

function trialMaster30k(df::DataFrame, times::DataFrame, trial::Int)

    # the actual trial slicer
    function get_trial(i)

        df1 = filter(df) do x
            x.time >= times[i,2] &&
            x.time <= times[i,4] &&
            x.conf > 0
        end

        return(df1)

    end

    # choose trial
    df_trial = get_trial(trial)

    # define a filter
    f = digitalfilter(Lowpass(5, fs = 200), Butterworth(4))

    # use the filter
    filt!(df_trial.gaze_point_3d_x ,f, df_trial.gaze_point_3d_x)
    filt!(df_trial.gaze_point_3d_y ,f, df_trial.gaze_point_3d_y)

## iterate through the eye tracking data and calculate the distance between two successive gaze points

    dist_x, dist_y = [0.0], [0.0]

    for i in 2:nrow(df_trial)
        append!(dist_x, df_trial.gaze_point_3d_x[i] - df_trial.gaze_point_3d_x[i-1])
        append!(dist_y, df_trial.gaze_point_3d_y[i] - df_trial.gaze_point_3d_y[i-1])
    end

    df_trial.dist_x = dist_x
    df_trial.dist_y = dist_y

## calculate angular change as in (Aks, Zielinski & Sprott, 2002)
    angchange = []
    for i in 1:nrow(df_trial)
        append!(angchange, atan(df_trial.dist_y[i], df_trial.dist_x[i]))
    end

    df_trial.angchange = angchange

## calculate angular position as..well..nowhere

    angles = []

    for i in 1:nrow(df_trial)
        append!(angles, atan(df_trial.gaze_point_3d_x[i], df_trial.gaze_point_3d_y[i]))
    end

    df_trial.angles = angles
## calculate plain ol' ecuclidian distances for two consecutive gazepoints

    edistance = [0.0]

    for i in 2:(nrow(df_trial))
        append!(edistance, euclidean(df_trial[i-1, 4:5], df_trial[i, 4:5]))
    end
    df_trial.edistance = edistance

## debug info and output
    println("The trialmaster has succeeded. Again.")
    return df_trial

end
