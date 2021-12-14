# the MassiveFastFourierTransformer will

using XDF, Plots, StatsKit, DataFrames, DynamicalSystems, ElasticArrays, LinearAlgebra, Distances, CSV, FFTW, StatsPlots
using GLM, Random

function MassiveFourierTranformer(df_trial::DataFrame;  windowsize::Int)
        # same as in IntergalacticRQAtor
        # could actually be stored as a stand alone method
        # anyway....constructs an array of window sizes

        S = convert(Vector{Float64},df_trial[!,:angchange])
        S = filter(S) do x
                        x != pi &&
                        x != -pi
        end


        triallength = size(S,1)
        windows =  []
        q = triallength

        while q > 100
                if q >= windowsize
                        prepend!(windows, windowsize)
                else
                        prepend!(windows, q)
                end
                q = q - windowsize
        end

        output = zeros(6)
        i = 0

        # this function does an FFT on every window
        # and the fits a slope to the power spectrum
        # and then fits three different distributions to the power spectrum
        # and then compare the fits with an Kolomogorov Smirnoff Test. This last step is complete BS.

        for window in windows
                start = triallength-window >= 1 ? triallength-window : 1
                stop = triallength
                s = S[start:stop]
## all the FFT magic
                fs = 200
                data = s
                F = fft(data) |> fftshift
                freqs = FFTW.fftfreq(length(data), fs) |> fftshift


                start = Int(floor(length(F)/2))
                stop = length(F)

                freqs = freqs[start+15:stop]
                Fspec = abs.(F)[start+15:stop]

                d = [log.(freqs), log.(Fspec)]
                fftout = DataFrame(d, [:freq, :power])

## all the slope magic

                model = @formula(freq ~ power)
                modelParams = lm(model, fftout)

                v = []
                i += 1
                push!(v, i)
                push!(v, coef(modelParams)[2])
                push!(v, coef(modelParams)[1])

## all the distribution fitting stuff

                expfit = fit_mle(Pareto, Fspec)
                parfit = fit_mle(Exponential, Fspec)
                gamfit = fit_mle(Normal, Fspec)

                test1 = ExactOneSampleKSTest(Fspec, expfit)
                test2 = ExactOneSampleKSTest(Fspec, parfit)
                test3 = ExactOneSampleKSTest(Fspec, gamfit)

                push!(v, pvalue(test1))
                push!(v, pvalue(test2))
                push!(v, pvalue(test3))

## output and debug per window

                output = [output v]
                println("Time series got FFFT'ed!!")

        end

        ## ouput and debug per trial

        output = transpose(output)
        output = output[2:end,:]
        return output
        println("Trial got FFT'ed. Pow.")
end
