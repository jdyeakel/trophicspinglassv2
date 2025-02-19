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
# Initialize 3D arrays
prop_changes = Array{Float64}(undef, lkoutvec, lgvec, reps);
prop_oscillating = Array{Float64}(undef, lkoutvec, lgvec, reps);

avg_time_pos1 = Array{Float64}(undef, lkoutvec, lgvec, reps);
avg_time_neg1 = Array{Float64}(undef, lkoutvec, lgvec, reps);

avg_dur_pos1 = Array{Float64}(undef, lkoutvec, lgvec, reps);
avg_dur_neg1 = Array{Float64}(undef, lkoutvec, lgvec, reps);

std_time_pos1 = Array{Float64}(undef, lkoutvec, lgvec, reps);
std_time_neg1 = Array{Float64}(undef, lkoutvec, lgvec, reps);

std_dur_pos1 = Array{Float64}(undef, lkoutvec, lgvec, reps);
std_dur_neg1 = Array{Float64}(undef, lkoutvec, lgvec, reps);

transition_pospos = Array{Float64}(undef, lkoutvec, lgvec, reps);
transition_posneg = Array{Float64}(undef, lkoutvec, lgvec, reps);
transition_negpos = Array{Float64}(undef, lkoutvec, lgvec, reps);
transition_negneg = Array{Float64}(undef, lkoutvec, lgvec, reps);

cascade_ratio_offset = Array{Float64}(undef, lkoutvec, lgvec, reps);
cascade_ratio = Array{Float64}(undef, lkoutvec, lgvec, reps);

propneg_effect_ratio = Array{Float64}(undef, lkoutvec, lgvec, reps);
propneg_effect_ratio_offset = Array{Float64}(undef, lkoutvec, lgvec, reps);

propneg_effect_interactions = Array{Float64}(undef, lkoutvec, lgvec, reps);
propneg_effect_interactions_offset = Array{Float64}(undef, lkoutvec, lgvec, reps);

@showprogress 1 "Computing..." for i in 1:lkoutvec
    @threads for j in 1:lgvec
        filename = smartpath(string("data/gvec_kout",subfile,"/g$(j)_kout$(i).jld2"))
        @load filename Asort_list tlsort_list s_list gprim kout kin sigma tmax
        # println([i,j])
        for r in 1:reps
            s = s_list[r][end-100:end, :]
            Asort = Asort_list[r]

            # State changes
            n_changes = compute_state_changes(s)
            prop_changes[i, j, r] = mean(n_changes / 101)
            prop_oscillating[i, j, r] = compute_proportion_oscillating(s)

            # Average time in each state
            avg_time_pos1[i, j, r] = mean(filter(!isnan,compute_total_time_in_states(s)[!,:TotalTime_pos1] ./ 101))
            avg_time_neg1[i, j, r] = mean(filter(!isnan,compute_total_time_in_states(s)[!,:TotalTime_neg1] ./ 101))

            std_time_pos1[i, j, r] = std(filter(!isnan,compute_total_time_in_states(s)[!,:TotalTime_pos1] ./ 101))
            std_time_neg1[i, j, r] = std(filter(!isnan,compute_total_time_in_states(s)[!,:TotalTime_neg1] ./ 101))

            # Average consecutive durations
            avg_durations = compute_average_durations(compute_state_durations(s))
            avg_dur_pos1[i, j, r] = mean(filter(!isnan,avg_durations[!,:pos1]))
            avg_dur_neg1[i, j, r] = mean(filter(!isnan,avg_durations[!,:neg1]))

            std_dur_pos1[i, j, r] = std(filter(!isnan,avg_durations[!,:pos1]))
            std_dur_neg1[i, j, r] = std(filter(!isnan,avg_durations[!,:neg1]))

            # Transition probabilities
            transition_probs = compute_state_transition_probabilities(s)
            transition_pospos[i, j, r] = mean(filter(!isnan,transition_probs[!,:pos1_pos1]))
            transition_posneg[i, j, r] = mean(filter(!isnan,transition_probs[!,:pos1_neg1]))
            transition_negpos[i, j, r] = mean(filter(!isnan,transition_probs[!,:neg1_pos1]))
            transition_negneg[i, j, r] = mean(filter(!isnan,transition_probs[!,:neg1_neg1]))

            # Cascade effects
            cascade_effect_ratio, cascade_effect_interactions = cascade_effect(s, Asort, 0)
            cascade_ratio[i, j, r] = mean(filter(!isnan, cascade_effect_ratio))
            propneg_effect_ratio[i, j, r] = sum(cascade_effect_ratio .< 0)/length(cascade_effect_ratio)
            propneg_effect_interactions[i, j, r] = sum(filter(!isnan,cascade_effect_interactions) .< 0)/length(filter(!isnan,cascade_effect_interactions))

            cascade_effect_ratio_offset, cascade_effect_interactions_offset = cascade_effect(s, Asort, 1)
            cascade_ratio_offset[i, j, r] = mean(filter(!isnan, cascade_effect_ratio_offset))
            propneg_effect_ratio_offset[i, j, r] = sum(cascade_effect_ratio_offset .< 0)/length(cascade_effect_ratio_offset)
            propneg_effect_interactions_offset[i, j, r] = sum(filter(!isnan,cascade_effect_interactions_offset) .< 0)/length(filter(!isnan,cascade_effect_interactions_offset))

        end
    end
end


# Compute mean across dimension 3 while ignoring NaNs
mean_prop_oscillating = mapslices(x -> mean(filter(!isnan, x)), prop_oscillating; dims=3)[:,:];
mean_cascade_ratio = mapslices(x -> mean(filter(!isnan, x)), cascade_ratio; dims=3)[:,:];
mean_cascade_ratio_offset = mapslices(x -> mean(filter(!isnan, x)), cascade_ratio_offset; dims=3)[:,:];

#Probability of a negative (-) effect:
#At Network Level
propneg_cascade_ratio = mapslices(x -> sum(filter(!isnan, x) .< 0) / reps, cascade_ratio; dims=3)[:,:];
propneg_cascade_ratio_offset = mapslices(x -> sum(filter(!isnan, x) .< 0) / reps, cascade_ratio_offset; dims=3)[:,:];
#At Species Level
propneg_cascade_ratio_species = mapslices(x -> mean(filter(!isnan, x)), propneg_effect_ratio; dims=3)[:,:];
propneg_cascade_ratio_species_offset = mapslices(x -> mean(filter(!isnan, x)), propneg_effect_ratio_offset; dims=3)[:,:];
#By interaction
propneg_cascade_ratio_int = mapslices(x -> mean(filter(!isnan, x)), propneg_effect_interactions; dims=3)[:,:];
propneg_cascade_ratio_int_offset = mapslices(x -> mean(filter(!isnan, x)), propneg_effect_interactions_offset; dims=3)[:,:];

#### Save output
filename = smartpath(string("data/gvec_kout",subfile,"/post_analysis.jld2"))
@save filename mean_prop_oscillating mean_cascade_ratio mean_cascade_ratio_offset propneg_cascade_ratio propneg_cascade_ratio_offset propneg_cascade_ratio_species propneg_cascade_ratio_species_offset propneg_cascade_ratio_int propneg_cascade_ratio_int_offset avg_dur_pos1 avg_dur_neg1 transition_pospos transition_posneg transition_negpos transition_negneg avg_time_pos1 avg_time_neg1 std_time_pos1 std_time_neg1 std_dur_pos1 std_dur_neg1






############# Some Plots ################




UnicodePlots.heatmap(mapslices(x -> mean(filter(!isnan, x)), avg_dur_neg1; dims=3)[:,:])

UnicodePlots.heatmap(mapslices(x -> mean(filter(!isnan, x)), transition_negneg; dims=3)[:,:])

UnicodePlots.heatmap(mapslices(x -> mean(filter(!isnan, x)), transition_posneg; dims=3)[:,:])

UnicodePlots.heatmap(mapslices(x -> mean(filter(!isnan, x)), avg_time_pos1; dims=3)[:,:])

StatsPlots.heatmap(mapslices(x -> mean(filter(!isnan, x)), avg_time_neg1; dims=3)[:,:] .> mapslices(x -> mean(filter(!isnan, x)), avg_time_pos1; dims=3)[:,:])

UnicodePlots.heatmap(mean_prop_oscillating)

StatsPlots.histogram(filter(!isnan,vec(cascade_ratio)))
StatsPlots.histogram(filter(!isnan,vec(cascade_ratio_offset)))

UnicodePlots.heatmap(mean_cascade_ratio)

UnicodePlots.heatmap(mean_cascade_ratio_offset)

UnicodePlots.scatterplot(vec(mean_cascade_ratio),vec(mean_cascade_ratio_offset))

UnicodePlots.heatmap(propneg_cascade_ratio)

UnicodePlots.heatmap(propneg_cascade_ratio_offset)

UnicodePlots.scatterplot(vec(propneg_cascade_ratio),vec(propneg_cascade_ratio_offset))

UnicodePlots.heatmap(propneg_cascade_ratio_species)
UnicodePlots.heatmap(propneg_cascade_ratio_species_offset)
UnicodePlots.heatmap(propneg_cascade_ratio_int)
UnicodePlots.heatmap(propneg_cascade_ratio_int_offset)




####

UnicodePlots.heatmap(prop_changes')
UnicodePlots.heatmap(prop_oscillating')

UnicodePlots.heatmap(avg_dur_pos1')
UnicodePlots.heatmap(avg_dur_neg1')

UnicodePlots.heatmap(transition_pospos')
UnicodePlots.heatmap(transition_posneg')
UnicodePlots.heatmap(transition_negpos')
UnicodePlots.heatmap(transition_negneg')

UnicodePlots.heatmap(mean(cascade_ratio, dims=3)[:,:])
UnicodePlots.heatmap(mean(cascade_ratio_offset, dims=3)[:,:])




UnicodePlots.heatmap(median(cascade_ratio, dims=3)[:,:])
UnicodePlots.heatmap(median(cascade_ratio_offset, dims=3)[:,:])




UnicodePlots.scatterplot(vec(cascade_ratio),vec(cascade_ratio_offset))
UnicodePlots.scatterplot(vec(mean(cascade_ratio, dims=3)[:,:]),vec(mean(cascade_ratio_offset, dims=3)[:,:]))
UnicodePlots.scatterplot(vec(median(cascade_ratio, dims=3)[:,:]),vec(median(cascade_ratio_offset, dims=3)[:,:]))


Plots.scatter(vec(mean(cascade_ratio, dims=3)[:,:]),vec(mean(cascade_ratio_offset, dims=3)[:,:]))
