
using Revise
using trophicspinglass

using StatsBase
using RCall
using Graphs
using IterTools
using Combinatorics
using GraphPlot
using Colors
using Random
using Distributions 
using DataFrames
using LinearAlgebra 
using JLD2 
using LightGraphs
using Base.Threads

using CSV
using Plots
using UnicodePlots
using StatsPlots
using LsqFit

C = 0.02



topreps = 100

reps = 1000;

mavg_tlspi = Array{Float64}(undef,0)
mavg_asymm_ratio = Array{Float64}(undef,0)
mavg_nichei = Array{Float64}(undef,0)

@threads for rr = 1:topreps
    nichei = Array{Float64}(undef,0)
    tlspi = Array{Float64}(undef,0)
    asymm_ratio = Array{Float64}(undef,0)

    for r = 1:reps
        # Generate network
        A = nothing  # Placeholder for adjacency matrix
        niche = nothing
        tl = Float64[]                    # Placeholder for trophic levels
        mintl = 0.0                       # Initial minimum trophic level

        while mintl < 1.0
            A, niche = nichemodelweb(S, C)
            tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6)
            mintl = minimum(tl)
        end
        
        tlsp = sortperm(tl)
        tlsort = tl[tlsp]
        Asort = A[tlsp, tlsp]
        nichesort = niche[tlsp]

        sp_withpreds = findall(!iszero,vec(sum(Asort,dims=2)))

        tldiff_pred = Array{Vector{Float64}}(undef,length(sp_withpreds))
        tlspivec = Array{Vector{Float64}}(undef,length(sp_withpreds))
        
        
        for i in eachindex(sp_withpreds)

            spi = sp_withpreds[i]

            prey_i = findall(x->x==1,Asort[:,spi])
            preds_i = findall(x->x==1,Asort[spi,:])

            tldiff_pred[i] = tlsort[preds_i] .- tlsort[spi]

            tlspivec[i] = repeat([tlsort[spi]],length(tldiff_pred[i]))

            if length(prey_i) > 0
                push!(nichei, nichesort[spi])
                push!(tlspi, tlsort[spi])
                push!(asymm_ratio,length(preds_i)/length(prey_i))
            end
        end

        # tldiff_pred_all = reduce(vcat,tldiff_pred)
        # tlspi_all = reduce(vcat,tlspi)
        # UnicodePlots.scatterplot(tlspi_all,tldiff_pred_all)

        # sum(tldiff_pred_all .< 0)/length(tldiff_pred_all)
        # mean(tldiff_pred_all)

    end
    
    avg_tlspi, avg_asymm_ratio = moving_average(tlspi,log.(asymm_ratio),0.05)
    
    Plots.scatter(tlspi,log.(asymm_ratio),alpha=0.5)
    Plots.plot!(moving_average(tlspi,log.(asymm_ratio),0.05),color=:red)

    avg_nichei, avg_asymm_ratio = moving_average(nichei,log.(asymm_ratio),0.05)

    Plots.scatter(nichei,log.(asymm_ratio),alpha=0.5)
    Plots.plot!(moving_average(nichei,log.(asymm_ratio),0.05),color=:red)
    expniche = collect(0:0.01:1);
    Plots.plot!(expniche,log.((1 .- expniche)./expniche),color=:black)
    
    Plots.scatter(nichei,tlspi)
    
    Plots.plot!(avg_tlspi,avg_asymm_ratio,color=:red)

    append!(mavg_tlspi,avg_tlspi)
    append!(mavg_asymm_ratio,avg_asymm_ratio)
    append!(mavg_nichei,avg_nichei)

end

# UnicodePlots.scatterplot(tlspi,asymm_ratio)

Plots.scatter(tlspi,log.(asymm_ratio),alpha=0.1)


Plots.scatter(mavg_tlspi,mavg_asymm_ratio,xlims=(1.5,5))
Plots.plot!(moving_average(mavg_tlspi,mavg_asymm_ratio,0.1),color=:red)


Plots.scatter(mavg_nichei,mavg_asymm_ratio,xlims=(0,1))
Plots.plot!(moving_average(mavg_nichei,mavg_asymm_ratio,0.1),color=:red)

Plots.scatter(mavg_nichei,mavg_tlspi,xlims=(0,1))
Plots.plot!(moving_average(mavg_nichei,mavg_tlspi,0.1),color=:red)
