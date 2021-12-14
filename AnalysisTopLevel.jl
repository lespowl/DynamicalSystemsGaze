using CSV, Statistics, MixedModels

include("IntergalacticRQAtor.jl")
include("MassiveFastFourierTransformer.jl")

## define an array of data set names do be included in the analysis
datasets = ["003", "004","005","006","007","008","009","010","011"]
datasets = ["012"]
datasets = ["013","014","015","016","017","018"]
datasets = ["008"]

## This is the full nelson. it will perform both, RQA and FFT on all the datasets in the array selected above

allout = DataFrame()
fixed, windowsize, indicator, ε = 0, 1200, "angchange", 10

for dataset in datasets
    out = DataFrame()

    # select dataset for subject
    df, times = subjectMaker2000(dataset)

    # loop through all the trials, return RQA measures and plots
    for trial in 1:21
        # slice trial from dataset
        df_trial = trialMaster30k(df, times, trial)

        # get columnnames ready (could also be done globally but ended up being here for histroric reasons)
        columnnames_RQA = [:Window, :Length, :DIV, :LAM, :Vmax, :VENTR, :Lmax, :MRT,:NMPRT,:RR,:RTE,:TT,:L,:ENTR,:DET,:TREND]
        columnnames_FIT = [:Window, :Slope, :Intercept]

        # store ouputs from RQA and FFT and join them into one DF
        RQAdata = DataFrame(RQAtor(df_trial; fixed = fixed, windowsize = windowsize, indicator = indicator), columnnames_RQA)
        FITdata =  DataFrame(MassiveFourierTranformer(df_trial; windowsize = windowsize), columnnames_FIT)
        outi = innerjoin(RQAdata, FITdata, on = :Window)

        outi.trial = zeros(nrow(outi))
        outi.trial .= trial

        append!(out, outi)

        # organize plotting if needed
        #fn = "output/" * dataset * string(trial)
        #plot(ploti)
        #savefig(fn)
        println("Trial number ", trial, " done.")
    end

    out.subject = zeros(nrow(out))
    out.subject .= tryparse(Int, dataset)


    # use this to prduce a CSV for every subject
    #fn = "output/rqa_out" * dataset * "fixed_1000.csv"
    #CSV.write(fn, out)

    # use that to put all the subjects into on CSV
    append!(allout, out)
    CSV.write("allet.csv", allout)

    println("Agent " * dataset * " successfully terminated.")
end

plot(allout.Slope)
CSV.write("../output/allofthem.csv", allout)

## Export for Tehran

 for dataset in datasets
     out = DataFrame()
     df, times = subjectMaker2000(dataset)

     for trial in 1:4
         df_trial = trialMaster30k(df, times, trial)

         outi = df_trial[!,[:conf, :gaze_point_3d_x, :gaze_point_3d_y, :angchange]]
         outi.trial = zeros(nrow(outi))
         outi.trial .= trial + 100
         append!(out, outi)
     end

     for trial in 5:25
         df_trial = trialMaster30k(df, times, trial)

         outi = df_trial[!,[:conf, :gaze_point_3d_x, :gaze_point_3d_y, :angchange]]
         outi.trial = zeros(nrow(outi))
         outi.trial .= trial - 4
         append!(out, outi)
     end
         fn = "raw/data_" * dataset * ".csv"
         CSV.write(fn, out)
end

## This section can be used to find out the length of a trial

datasets = ["002","003", "004","005","006","007","008","009","010","011", "012", "013","014","015","016","017","018"]

out = DataFrame(ID=String[], Length=Int[])

for dataset in datasets
    df, times = subjectMaker2000(dataset)

     v = 0

    for trial in 1:21
        df_trial = trialMaster30k(df, times, trial)
        v = v + size(df_trial, 1)
    end
    new = [dataset, v]
    push!(out, new)

end

triallength = mean(out.Length)/200/60

## This section can be used to perform FFT on a set of subjects/trials

include("MassiveFastFourierTransformer.jl")

datasets = ["002","003", "004","005","006","007","008","009","010","011", "012", "013","014","015","016","017","018"]

allout = DataFrame()
fixed, windowsize, indicator, ε = 0, 100000, "angchange", 10

for dataset in datasets
    out = DataFrame()
    # select dataset
    df, times = subjectMaker2000(dataset)


    # loop through all the trials, return RQA measures and plots
    for trial in 1:21
        df_trial = trialMaster30k(df, times, trial)
        columnnames_FIT = [:Window, :Slope, :Intercept, :p1, :p2, :p3]
        outi =  DataFrame(MassiveFourierTranformer(df_trial; windowsize = 50000), columnnames_FIT)
        #for key in keyset push!(columnnames,  key) end
        outi.trial = zeros(nrow(outi))
        outi.trial .= trial
        append!(out, outi)

        println("Trial number ", trial, " done.")
    end

    out.subject = zeros(nrow(out))
    out.subject .= tryparse(Int, dataset)

    append!(allout, out)

    #fn = "output/rqa_out" * dataset * "fixed_1000.csv"
    #CSV.write(fn, out)
    CSV.write("../output/allet.csv", allout)
    println("Agent " * dataset * " successfully terminated.")

end
