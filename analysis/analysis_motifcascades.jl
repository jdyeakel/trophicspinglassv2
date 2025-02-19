using Revise
using trophicspinglass

using RCall
using Graphs
using GraphPlot
using Colors
using Random
using Distributions 
using DataFrames
using LinearAlgebra 
using JLD2 
using LightGraphs
using Base.Threads

using Plots
using UnicodePlots
using StatsPlots
using StatsBase
using ProgressMeter
using IterTools

subfile = "_lowC"
subfile = ""
subfile = "_highC"
subfile = "_veryhighC"
subfile = "_veryhigh2xC"

# Load variables from saved file
filename_preamble = smartpath(string("data/gvec_kout",subfile,"/preamble.jld2"))
@load filename_preamble reps gvec koutvec S C kin tmax sigma

lgvec = length(gvec);
lkoutvec = length(koutvec);

motifs = motifs_generate();

mean_predenriched = Array{Float64}(undef, lkoutvec, lgvec, reps);

@showprogress 1 "Computing..." for i in 1:lkoutvec
    # @threads for j in 1:lgvec

    #Choose a value of g
    j = 5;
    println(gvec[j])

    filename = smartpath(string("data/gvec_kout",subfile,"/g$(j)_kout$(i).jld2"))
    @load filename Asort_list tlsort_list s_list gprim kout kin sigma tmax
    # println([i,j])

    #WHICH MOTIF?
    motiflabel = :m2

    maxtl_motif_reps = Array{Vector{Float64}}(undef,reps)
    meanratio_cascade_motif_reps = Array{Vector{Float64}}(undef,reps)

    @threads for r in 1:reps
        
        s = s_list[r][end-100:end, :]
        Asort = Asort_list[r]
        tlsort = tlsort_list[r];
        g = Graphs.DiGraph(Asort);

        #find the 2-chain interactions
        emb_motif = find_motifs(g, motifs[motiflabel])
        
        # find predator and associated tl for each pair
        pred_id_tl = Vector{Tuple{Int64,Float64}}(undef,length(emb_motif));  # Initialize predtl as an empty array
        # proppos = Array{Float64}(undef,length(predpreymotif),2);
        predenriched = Array{Float64}(undef,length(emb_motif));
        
        meanratio_cascade_motif = Array{Float64}(undef,length(emb_motif));
        maxtl_motif = Array{Float64}(undef,length(emb_motif));
        for k = 1:length(emb_motif)

            #Species in the motif
            spmotif = emb_motif[k]

            maxtl_motif[k] = maximum(tlsort[spmotif])
            
            #List each species' predators present in the current motif
            sp_preds = [intersect(spmotif,findall(x->x==1,Asort[spmotif[l],:])) for l=1:length(spmotif)]

            cascade_ratio_motif = Array{Float64}(undef,0)
            predpreypairs = 0;
            for l=1:length(spmotif)
                if length(sp_preds[l]) > 0
                    preds = sp_preds[l]
                    

                    offset = 1;
                    coeff = 2;
                    
                    
                    cascade_ratio_perpred = Array{Float64}(undef,length(preds))
                    for p=1:length(preds)

                        predpreypairs += 1
                        
                        # +P positions
                        posP = findall(x->x==1,s[1:(end-offset),preds[p]])
                        # -P positions
                        negP = setdiff(collect(1:1:length(posP)),posP)

                        avg_state_posP = mean(s[posP .+ offset,spmotif[l]])
                        avg_state_negP = mean(s[negP .+ offset,spmotif[l]])
                        cascade_ratio_perpred = log((coeff + avg_state_posP) / (coeff + avg_state_negP))
                        push!(cascade_ratio_motif,cascade_ratio_perpred)
                    end
                    # propneg_cascade_ratio_motif[l] = sum(filter(!isnan,cascade_ratio_perpred) .< 0) / length(spmotif)
                else
                    push!(cascade_ratio_motif,NaN)
                end
            end

            # propneg_cascade_motif[k] = sum(filter(!isnan,cascade_ratio_motif) .< 0) / predpreypairs
            meanratio_cascade_motif[k] = mean(filter(!isnan,cascade_ratio_motif))

        end 

        # UnicodePlots.scatterplot(maxtl_motif,propneg_cascade_motif)

        #This may not be thread safe!!!!
        maxtl_motif_reps[r] = maxtl_motif
        meanratio_cascade_motif_reps[r] = meanratio_cascade_motif

        # append!(maxtl_motif_reps,maxtl_motif)
        # append!(meanratio_cascade_motif_reps,meanratio_cascade_motif)

        # checkmin = findmin(maxtl_motif)[1]

        # if checkmin < 2.0
        #     println(r)
        # end

    end #end reps

    #collapse
    maxtl_motif_flattened = reduce(vcat, maxtl_motif_reps);
    meanratio_cascade_motif_flattened =  reduce(vcat, meanratio_cascade_motif_reps);

    UnicodePlots.scatterplot(maxtl_motif_flattened,meanratio_cascade_motif_flattened)

    maxtl_avg, cascadeindex_avg = moving_average(maxtl_motif_flattened .- 1,meanratio_cascade_motif_flattened,0.1);
    
    UnicodePlots.lineplot(maxtl_avg,cascadeindex_avg)

    # end
end



avgreps_predenriched = mapslices(x -> mean(filter(!isnan, x)), mean_predenriched; dims=3)[:,:];

UnicodePlots.heatmap(avgreps_predenriched)