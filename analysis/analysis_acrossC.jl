

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
using ColorSchemes

subfilelist = ["_lowC", "", "_highC", "_veryhighC", "_veryhigh2xC"]

# Create a dictionary to store results for each run
results = Dict()

for subfile in subfilelist
    # Load the file
    filename = smartpath(string("data/gvec_kout", subfile, "/post_analysis.jld2"))

    # Initialize a sub-dictionary for this specific run
    run_data = Dict()

    # Load variables into the sub-dictionary
    # @load filename mean_prop_oscillating mean_cascade_ratio mean_cascade_ratio_offset propneg_cascade_ratio propneg_cascade_ratio_offset propneg_cascade_ratio_species propneg_cascade_ratio_species_offset propneg_cascade_ratio_int propneg_cascade_ratio_int_offset
    # run_data[:mean_prop_oscillating] = mean_prop_oscillating
    # run_data[:mean_cascade_ratio] = mean_cascade_ratio
    # run_data[:mean_cascade_ratio_offset] = mean_cascade_ratio_offset
    # run_data[:propneg_cascade_ratio] = propneg_cascade_ratio
    # run_data[:propneg_cascade_ratio_offset] = propneg_cascade_ratio_offset
    # run_data[:propneg_cascade_ratio_species] = propneg_cascade_ratio_species
    # run_data[:propneg_cascade_ratio_species_offset] = propneg_cascade_ratio_species_offset
    # run_data[:propneg_cascade_ratio_int] = propneg_cascade_ratio_int
    # run_data[:propneg_cascade_ratio_int_offset] = propneg_cascade_ratio_int_offset

    
    # Load variables into the sub-dictionary
    @load filename mean_prop_oscillating mean_cascade_ratio mean_cascade_ratio_offset propneg_cascade_ratio propneg_cascade_ratio_offset propneg_cascade_ratio_species propneg_cascade_ratio_species_offset propneg_cascade_ratio_int propneg_cascade_ratio_int_offset avg_dur_pos1 avg_dur_neg1 transition_pospos transition_posneg transition_negpos transition_negneg avg_time_pos1 avg_time_neg1 std_time_pos1 std_time_neg1 std_dur_pos1 std_dur_neg1

    run_data[:mean_prop_oscillating] = mean_prop_oscillating
    run_data[:mean_cascade_ratio] = mean_cascade_ratio
    run_data[:mean_cascade_ratio_offset] = mean_cascade_ratio_offset
    run_data[:propneg_cascade_ratio] = propneg_cascade_ratio
    run_data[:propneg_cascade_ratio_offset] = propneg_cascade_ratio_offset
    run_data[:propneg_cascade_ratio_species] = propneg_cascade_ratio_species
    run_data[:propneg_cascade_ratio_species_offset] = propneg_cascade_ratio_species_offset
    run_data[:propneg_cascade_ratio_int] = propneg_cascade_ratio_int
    run_data[:propneg_cascade_ratio_int_offset] = propneg_cascade_ratio_int_offset
    run_data[:avg_dur_neg1] = avg_dur_neg1
    run_data[:avg_dur_pos1] = avg_dur_pos1
    run_data[:transition_pospos] = transition_pospos
    run_data[:transition_posneg] = transition_posneg
    run_data[:transition_negpos] = transition_negpos
    run_data[:transition_negneg] = transition_negneg
    run_data[:avg_time_pos1] = avg_time_pos1
    run_data[:avg_time_neg1] = avg_time_neg1
    run_data[:std_time_pos1] = std_time_pos1
    run_data[:std_time_neg1] = std_time_neg1
    run_data[:std_dur_pos1] = std_dur_pos1
    run_data[:std_dur_neg1] = std_dur_neg1


    # Use the subfile (or a cleaned version of it) as the key
    results[subfile] = run_data
end

#Systemic analysis - early in the paper
#Proportion cycling
UnicodePlots.heatmap(log.(mean(results[""][:avg_dur_neg1],dims=3)[:,:,1]))

UnicodePlots.heatmap(log.(mean(results[""][:avg_dur_pos1],dims=3)[:,:,1]))

UnicodePlots.scatter(log.(vec(results[""][:avg_dur_pos1])),log.(vec(results[""][:avg_dur_neg1])))

UnicodePlots.scatterplot(vec(results[""][:avg_time_pos1]),log.(vec(results[""][:avg_dur_neg1])))

Plots.scatter(vec(results[""][:avg_time_pos1]),(vec(results[""][:avg_dur_neg1]))./101,
    # xlims=(0,0.3),
    # ylims=(0,0.3),
    xlabel="1 - Probability S-",
    ylabel="Duration S-",
    alpha=0.3,
    markerstrokewidth=0,
    label=false,
    frame=:box)

expected_state = mean.(1*vec(results[""][:avg_time_pos1]) + (-1)*vec(results[""][:avg_time_neg1]));

q = vec(results[""][:avg_time_neg1]);
std_state = 2 .* sqrt.((1 .- q).*q);

Plots.scatter(expected_state,(vec(results[""][:std_time_neg1]))./101,
    # xlims=(0,0.3),
    # ylims=(0,0.3),
    xlabel="Expected State",
    ylabel="Std Time S-",
    alpha=0.3,
    markerstrokewidth=0,
    label=false,
    frame=:box)


Plots.scatter(expected_state,(vec(results[""][:avg_dur_neg1]))./101,
    # xlims=(0,0.3),
    # ylims=(0,0.3),
    xlabel="Expected State",
    ylabel="Duration S-",
    alpha=0.3,
    markerstrokewidth=0,
    label=false,
    frame=:box)

Plots.scatter(expected_state,(vec(results[""][:std_dur_neg1]))./101,
    # xlims=(0,0.3),
    # ylims=(0,0.3),
    xlabel="Expected State",
    ylabel="Std Duration S-",
    alpha=0.3,
    markerstrokewidth=0,
    label=false,
    frame=:box)

UnicodePlots.heatmap(mean(results[""][:std_dur_pos1],dims=3)[:,:,1])

Plots.scatter(std_state,(vec(results[""][:std_dur_neg1]))./101,
    # xlims=(0,0.3),
    # ylims=(0,0.3),
    xlabel="Std State",
    ylabel="Std Duration S-",
    alpha=0.3,
    markerstrokewidth=0,
    label=false,
    frame=:box)


using LaTeXStrings
using PGFPlotsX
pgfplotsx()
plotoscillating = Plots.heatmap(results[""][:mean_prop_oscillating],
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Proportion~cycling}",
    colorbar_titlefontsize = 14,
    color=:thermal);
filename = smartpath("figures/fig_proposcillating.pdf");
Plots.savefig(plotoscillating,filename)


UnicodePlots.heatmap(mean(results[""][:avg_time_pos1],dims=3)[:,:,1])


UnicodePlots.scatterplot(vec(results[""][:avg_dur_neg1]),vec(results[""][:avg_time_pos1]))

using Measures
using LaTeXStrings
using PGFPlotsX
pgfplotsx()
dur_pos = Plots.heatmap(mean(results[""][:avg_time_pos1],dims=3)[:,:,1],
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Duration,}~S^+",
    colorbar_titlefontsize = 14,
    color=:thermal,
    clim=(0, 1),
    left_margin=15mm,  # Add more space to the left
    bottom_margin=5mm);
dur_neg = Plots.heatmap(mean(results[""][:avg_time_neg1],dims=3)[:,:,1],
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Duration,}~S^-",
    colorbar_titlefontsize = 14,
    color=:thermal,
    clim=(0,1),
    left_margin=15mm,  # Add more space to the left
    bottom_margin=5mm);    

durplot = Plots.plot(dur_pos, dur_neg, layout=(1, 2), size=(1000, 400))
filename = smartpath("figures/fig_duration.pdf");
Plots.savefig(durplot,filename)


#Structural Analysis - later in the paper
#Comparing cascade ratios
p1 = UnicodePlots.scatterplot(vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_veryhigh2xC"][:propneg_cascade_ratio_species_offset]))
UnicodePlots.scatterplot!(p1,vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_veryhighC"][:propneg_cascade_ratio_species_offset]))
UnicodePlots.scatterplot!(p1,vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_highC"][:propneg_cascade_ratio_species_offset]))
UnicodePlots.lineplot!(p1,collect(0:0.1:1),collect(0:0.1:1))

p1fig = Plots.scatter(vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_veryhigh2xC"][:propneg_cascade_ratio_species_offset]),
    xlims=(0,0.3),
    ylims=(0,0.3),
    xlabel="Prob(Cascade); C=0.01",
    ylabel="Prob(Cascade); C=Calt",
    alpha=0.3,
    markerstrokewidth=0,
    label="Calt=0.04",
    frame=:box);
Plots.scatter!(p1fig,vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_veryhighC"][:propneg_cascade_ratio_species_offset]),
    alpha=0.3,
    markerstrokewidth=0,
    label="Calt=0.03");
Plots.scatter!(p1fig,vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_highC"][:propneg_cascade_ratio_species_offset]),
    alpha=0.3,
    markerstrokewidth=0,
    label="Calt=0.025");
Plots.plot!(p1fig,collect(0:0.1:1),collect(0:0.1:1),label=false,color=:black);
p1fig

filename = smartpath("figures/fig_probcascade_conn.pdf")
Plots.savefig(p1fig,filename)



p2 = UnicodePlots.scatterplot(vec(results["_lowC"][:propneg_cascade_ratio_int_offset]),vec(results["_veryhigh2xC"][:propneg_cascade_ratio_int_offset]))
UnicodePlots.scatterplot!(p2,vec(results["_lowC"][:propneg_cascade_ratio_int_offset]),vec(results["_veryhighC"][:propneg_cascade_ratio_int_offset]))
UnicodePlots.scatterplot!(p2,vec(results["_lowC"][:propneg_cascade_ratio_int_offset]),vec(results["_highC"][:propneg_cascade_ratio_int_offset]))
UnicodePlots.lineplot!(p2,collect(0:0.1:1),collect(0:0.1:1))


#what is the effect of connectance across kout and g?
ratio = results["_veryhigh2xC"][:propneg_cascade_ratio_int_offset] ./ results["_lowC"][:propneg_cascade_ratio_int_offset];

Plots.heatmap(results["_veryhigh2xC"][:propneg_cascade_ratio_int_offset])

Plots.heatmap(results["_lowC"][:propneg_cascade_ratio_int_offset])

Plots.heatmap((ratio),clim=(0,2))

Plots.heatmap((results["_veryhigh2xC"][:propneg_cascade_ratio_int_offset] .- results["_lowC"][:propneg_cascade_ratio_int_offset]).^2)

Plots.scatter(vec(results["_lowC"][:propneg_cascade_ratio_species_offset]),vec(results["_veryhigh2xC"][:propneg_cascade_ratio_species_offset]))


mean_cascade_ratio_offset

Plots.scatter(vec(results["_lowC"][:mean_cascade_ratio_offset]),vec(results["_veryhigh2xC"][:mean_cascade_ratio_offset]))
Plots.heatmap(results["_lowC"][:mean_cascade_ratio_offset])
Plots.heatmap(results["_veryhigh2xC"][:mean_cascade_ratio_offset])