using XDF, Plots, StatsKit, DataFrames, DynamicalSystems, ElasticArrays, LinearAlgebra, Distances, CSV, FFTW, StatsPlots, DSP

include("subjectMaker2000.jl")
include("trialMaster30k.jl")

##
# this is the all in one RQA solution to get an RQA
# of a specific trial from a specific subject
function IntergalacticRQAtor(subject::String, trial::Int; maxlag::Int = 200)

        # get the data ready and
        df, times = subjectMaker2000(subject)
        df_trial = trialMaster30k(df, times, trial)
        maxlag = maxlag
        println("Data parsed.")
        RQAtor(df_trial)
end

## this is the regular function also used in the analysis
function RQAtor(
        df_trial::DataFrame;
        maxlag::Int = 100, fixed::Int = 0, windowsize::Int = 100000, indicator::String="angles", ε::Int = 3
        )

        # Use edistance or angles or diameter or distvector...defined by the keyword 'indicator'
        S = convert(Vector{Float64},df_trial[!,indicator])

        #remove the finite values
        S = filter(S) do x
                x != pi &&
                x != -pi
        end

        # chop suey!
        # cut the time series S into windows of site 'windowsize'

        # check how long the trial is and initialize everything for chopping
        triallength = size(S,1)

        windows =  []
        q = triallength

        # generate an array of windows

        while q > 100
                if q >= windowsize
                        prepend!(windows, windowsize)
                else
                        prepend!(windows, q)
                end
                q = q - windowsize
        end

        # initialize output df and iterator
        output = zeros(16)
        i = 0

        # RQA element wise throug the windows array
        for window in windows
                start = triallength-window >= 1 ? triallength-window : 1
                stop = triallength
                s = S[start:stop]

                # debug feature
                println(i)

                triallength -= window

                ncols = size(s,2);

                # parameters:
                maxlag = ((size(s,1) * .15) <= 500 ? (size(s,1) * .15) : 500) |> round |> Int;

                # plotting features. used for diagnosis but disabled in analysis for speed
                # plotMI = zeros(maxlag+1, ncols);
                # plotMI = selfmutualinfo(s, 1:maxlag) # |> plot
                # τ = estimate_delay(s, "mi_min", 1:maxlag, binwidth = 6)
                # D = delay_afnn(s, τ)

                # This embeds the time series
                #Comment: changes in τs and dmax produce chatotic changes in the output...
                embetty = optimal_traditional_de(s, "afnn", "mi_min"; τs = 1:maxlag, dmax = 20)
                println("Embedding done.")

                # This retrieves the parameters for RQA
                embed_TS = embetty[1]
                τ = embetty[2]
                D = size(embetty[1],2)

                # construct Recurrence Matrix
                # fixedrate leads to occassional crashes of julia for timeseries greater than 15k samples
                R = fixed == 1 ? RecurrenceMatrix(embed_TS, 0.05 ;fixedrate = true) : R = RecurrenceMatrix(embed_TS, ε)

                # recurrence quantification
                rqaOUT = rqa(R; lmin = 10, theiler = 2)

                #organize dict into df
                global keyset = keys(rqaOUT)

                # organize output for the current window
                v = []
                i += 1
                push!(v, i)
                push!(v, stop-start)
                for j in keyset
                        push!(v, rqaOUT[j])
                end

                output = [output v]

                println("Window done.")

        end

        output = transpose(output)
        output = output[2:end,:]

        println("Trial done.")

        return output
end


        # Step 7 ------- Construct output data frame
#
#        ami_plot = plot(plotMI[:,:]);
#        plot!(ami_plot, [τ], seriestype="vline", linestyle=:dash,label = "delay selected")
#
#        # false nearest neighbors plot
#       fnn_plot = plot(fnnPerc, label = "%FNN");
#        plot!(fnn_plot, [D], seriestype="vline", linestyle=:dash,label = "embed selected")
#
#        # recurrence plot
#        xs, ys = coordinates(R);
#        rec_plot = scatter(xs, ys, markersize = 0.1, markercolor = :black, legend = false);
        #plot!(rec_plot, size = (1000,1000))
#
#        # combine plots
#        l = @layout [a{0.5w} b]
#        diag_plots = plot(rec_plot, ami_plot, layout = l);
#        plot!(diag_plots, size=(800,400))
#        println("Plotting done.")

#        return output, diag_plots, length(s)
