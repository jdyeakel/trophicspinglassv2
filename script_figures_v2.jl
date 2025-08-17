#=

THESE ARE REVISED FIGURES/CODE IN RESPONSE TO MANUSCRIPT REVISIONS

A note that this script reproduces the figures in the manuscript, but is not written incredibly efficiently, as these scripts were compiled from a larger assortment of scripts to investigate the model at different times... so there is a lot of duplication of computational effort.

Note also that Figures 5 and 6 are constructed using /math/analytical_solutions.nb as well as simulation export in the below script.

I have not included scripts for the supplementary figures.

I use the function src/smartpath() to direct figure and data files to a central directory. It's currently set for the path on my computer, but you will want to replace with the path to your file repository.


Coarse-Graining Cascades Within Food Webs

Quantifying population dynamics is a fundamental challenge in ecology and evolutionary biology, particularly for species that are cryptic, microscopic, or extinct. Traditional approaches rely on continuous representations of population size, but in many cases, the precise number of individuals is unknowable. Here, we present a coarse-grained population model that simplifies population dynamics to binary states -- high or low -- determined by the balance of bottom-up resource availability and top-down predation pressure. This Boolean framework provides a minimal yet analytically tractable alternative to traditional Lotka-Volterra-based models, enabling direct insights into the role of food web structure in shaping community stability. Using this approach, we investigate how trophic interactions influence population persistence, cyclic dynamics, and extinction risk across model food webs. We find that top-down effects are a primary driver of cycling, aligning with theoretical expectations from traditional population models, and that trophic position strongly influences extinction risk, with higher-trophic species more prone to persistent low-population states. Additionally, we explore the role of trophic short-circuits -- direct interactions between apex predators and low-trophic prey -- and find that they can buffer cascades and alter extinction patterns in ways that are often overlooked in classical models. By simplifying population dynamics to a two-state system, this framework provides a powerful tool for disentangling the structural drivers of community stability. These results highlight the potential of coarse-grained approaches to complement existing models, offering new insights into trophic interactions, extinction risks, and the susceptibility of species to trophic cascades.
=#



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
using ProgressMeter
using CSV
using Printf



#########################################
######### FIGURES 2, 3 ##################
#########################################

#NOTE: for Figures 2/3, fix the indexing and 'oldcode' issue


subfilelist = ["_lowC", "", "_highC", "_veryhighC", "_veryhigh2xC"]
Cnorm = 0.02;
Cvec = [Cnorm*(1 - 0.25), Cnorm, Cnorm*(1 + 0.25), Cnorm*(1 + 0.75), Cnorm*(1 + 1.)]

reps = 200
gvec = collect(0:1:25)
lgvec = length(gvec)
koutvec = collect(0:0.5:15)
lkoutvec = length(koutvec)

S = 100; 
kin = 5.0
tmax = 300
sigma = 0.01

for c in 1:length(Cvec)

    C = Cvec[c];
    subfile = subfilelist[c];

    # Optional: Save preamble or settings if needed
    filename_preamble = smartpath(string("data/gvec_kout",subfile,"/preamble.jld2"))

    @save filename_preamble reps gvec koutvec S C kin tmax sigma

    # Loop over parameter combinations
    # @showprogress 1 "Computing..." 
    for i in 1:lgvec
        @threads for j in 1:lkoutvec
            gprim = gvec[i]
            kout = koutvec[j]

            # Prepare to collect outputs from all repetitions
            Asort_list = Vector{Matrix{Int64}}(undef, reps)
            tlsort_list = Vector{Vector{Float64}}(undef, reps)
            s_list = Vector{Matrix{Int64}}(undef, reps)

            # Use threading over repetitions
            for r in 1:reps
                
                # Generate network
                A = nothing  # Placeholder for adjacency matrix
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

                # Run simulation
                s = cascade(Asort, kout, kin, gprim, sigma, tmax)

                # Collect outputs
                Asort_list[r] = Asort
                tlsort_list[r] = tlsort
                s_list[r] = s
            end

            # After all reps, save data to a single file

            filename = smartpath(string("data/gvec_kout",subfile,"/g$(i)_kout$(j).jld2"))
            @save filename Asort_list tlsort_list s_list gprim kout kin sigma tmax


        end

        println(i/lgvec)

    end
end


# CALCULATE STATUSICS FROM SIMULATIONS

for c=1:length(Cvec)
    subfile = subfilelist[c];
    # Load variables from saved file
    filename_preamble = smartpath(string("data/gvec_kout",subfile,"/preamble.jld2"))
    @load filename_preamble reps gvec koutvec S C kin tmax sigma

    lgvec = length(gvec);
    lkoutvec = length(koutvec);
    # Initialize 3D arrays
    prop_changes = Array{Float64}(undef, lgvec, lkoutvec, reps);
    prop_oscillating = Array{Float64}(undef, lgvec, lkoutvec, reps);

    avg_time_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    avg_time_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    avg_dur_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    avg_dur_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    std_time_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    std_time_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    std_dur_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    std_dur_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    transition_pospos = Array{Float64}(undef, lgvec, lkoutvec, reps);
    transition_posneg = Array{Float64}(undef, lgvec, lkoutvec, reps);
    transition_negpos = Array{Float64}(undef, lgvec, lkoutvec, reps);
    transition_negneg = Array{Float64}(undef, lgvec, lkoutvec, reps);

    cascade_ratio_offset = Array{Float64}(undef, lgvec, lkoutvec, reps);
    cascade_ratio = Array{Float64}(undef, lgvec, lkoutvec, reps);

    propneg_effect_ratio = Array{Float64}(undef, lgvec, lkoutvec, reps);
    propneg_effect_ratio_offset = Array{Float64}(undef, lgvec, lkoutvec, reps);

    propneg_effect_interactions = Array{Float64}(undef, lgvec, lkoutvec, reps);
    propneg_effect_interactions_offset = Array{Float64}(undef, lgvec, lkoutvec, reps);

    @showprogress 1 "Computing..." for i in 1:lgvec
        @threads for j in 1:lkoutvec
            filename = smartpath(string("data/gvec_kout",subfile,"/g$(i)_kout$(j).jld2"))
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

            end
        end
    end

    # Compute mean across dimension 3 while ignoring NaNs
    mean_prop_oscillating = mapslices(x -> mean(filter(!isnan, x)), prop_oscillating; dims=3)[:,:];
    mean_cascade_ratio = mapslices(x -> mean(filter(!isnan, x)), cascade_ratio; dims=3)[:,:];
    mean_cascade_ratio_offset = mapslices(x -> mean(filter(!isnan, x)), cascade_ratio_offset; dims=3)[:,:];

    #### Save output
    # newdir = replace(smartpath(string("data/gvec_kout",subfile,".jl")), ".jl" => "/");
    # mkpath(newdir);
    filename = smartpath(string("data/gvec_kout",subfile,"/post_analysis.jld2"))

    @save filename mean_prop_oscillating avg_dur_pos1 avg_dur_neg1 transition_pospos transition_posneg transition_negpos transition_negneg avg_time_pos1 avg_time_neg1 std_time_pos1 std_time_neg1 std_dur_pos1 std_dur_neg1

end


#OR LOAD DATAFILES DIRECTLY
subfilelist = ["_lowC", "", "_highC", "_veryhighC", "_veryhigh2xC"]
# filename_preamble = smartpath(string("data/gvec_kout",subfile,"/preamble.jld2"))
# @load filename_preamble reps gvec koutvec S C kin tmax sigma

# Create a dictionary to store results for each run
results = Dict()
for subfile in subfilelist
    # Load the file
    filename = smartpath(string("data/gvec_kout", subfile, "/post_analysis.jld2"))

    # Initialize a sub-dictionary for this specific run
    run_data = Dict()

    @load filename mean_prop_oscillating avg_dur_pos1 avg_dur_neg1 transition_pospos transition_posneg transition_negpos transition_negneg avg_time_pos1 avg_time_neg1 std_time_pos1 std_time_neg1 std_dur_pos1 std_dur_neg1

    run_data[:mean_prop_oscillating] = mean_prop_oscillating
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


##################################
######### FIGURE 2 ###############
##################################

using LaTeXStrings
using PGFPlotsX
pgfplotsx()
plotoscillating = Plots.heatmap(results[""][:mean_prop_oscillating]',
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Proportion~cycling}",
    colorbar_titlefontsize = 14,
    color=:thermal);
filename = smartpath("figures_prefinal/fig_proposcillating.pdf");
Plots.savefig(plotoscillating,filename)


##################################
######### FIGURE 3 ###############
##################################

using Measures
using LaTeXStrings
using PGFPlotsX
pgfplotsx()
dur_pos = Plots.heatmap(mean(results[""][:avg_time_pos1],dims=3)[:,:,1]',
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Proportional~residence~time,}~S^+",
    colorbar_titlefontsize = 14,
    color=:thermal,
    clim=(0, 1),
    left_margin=15mm,  # Add more space to the left
    bottom_margin=5mm);
dur_neg = Plots.heatmap(mean(results[""][:avg_time_neg1],dims=3)[:,:,1]',
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Proportional~residence~time,}~S^-",
    colorbar_titlefontsize = 14,
    color=:thermal,
    clim=(0,1),
    left_margin=15mm,  # Add more space to the left
    bottom_margin=5mm);    

durplot = Plots.plot(dur_pos, dur_neg, layout=(1, 2), size=(1000, 400));

#Save Panel A only
filename = smartpath("figures_prefinal/fig_avgtime_rev.pdf");
Plots.savefig(dur_pos,filename)



####################################
######### FIGURE S?? ###############
####################################

plotcycleperiod = Plots.heatmap((mean(results[""][:avg_dur_pos1],dims=3)[:,:]./101)',
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Cycle~period}",
    colorbar_titlefontsize = 14,
    color=:thermal)
filename = smartpath("figures_prefinal/fig_cycleperiod.pdf");
Plots.savefig(plotcycleperiod,filename)




##################################
######### FIGURE 4 ###############
##################################

#Run for a specific implementation and for C = 0.02 and C= 0.1
kin = 5.;
kout = 10.;
gprim = 10.; #try 20

S=100;
C_values = [0.02,0.04,0.1] # or a range like 0.01:0.01:0.2
tmax = 300
sigma = 0.01
reps = 5000

n_all = Float64[]
C_all = Float64[]
propneg_all = Float64[]
propneg_std_all = Float64[]
dur_pos1_all = Float64[]
dur_neg1_all = Float64[]
trophic_all = Float64[]

for C in C_values

    # Prepare to collect outputs from all repetitions
    Asort_r = Vector{Matrix{Int64}}(undef, reps);
    tlsort_r = Vector{Vector{Float64}}(undef, reps);
    nichesort_r = Vector{Vector{Float64}}(undef, reps);
    s_r = Vector{Matrix{Int64}}(undef, reps);

    @threads for r = 1:reps

        # Generate network
        A = nothing  # Placeholder for adjacency matrix
        tl = Float64[]                    # Placeholder for trophic levels
        mintl = 0.0                       # Initial minimum trophic level
        niche = Float64[]

        while mintl < 1.0
            A, niche = nichemodelweb(S, C)
            tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6)
            mintl = minimum(tl)
        end
        # g = Graphs.DiGraph(A);

        tl = round.(TrophInd(Array{Float64}(A))[!, :TL]; sigdigits=6)

        tlsp = sortperm(tl)
        tlsort = tl[tlsp]
        Asort = A[tlsp, tlsp]
        nichesort = niche[tlsp]

        # Run simulation
        s = cascade(Asort, kout, kin, gprim, sigma, tmax)

        # Collect outputs
        Asort_r[r] = Asort
        tlsort_r[r] = tlsort
        nichesort_r[r] = nichesort
        s_r[r] = s

    end

    propneg_r = Array{Vector{Float64}}(undef,reps);
    niche_r = Array{Vector{Float64}}(undef,reps);
    trophic_r = Array{Vector{Float64}}(undef,reps);
    avg_dur_pos1_r = Array{Vector{Float64}}(undef,reps);
    avg_dur_neg1_r = Array{Vector{Float64}}(undef,reps);

    #What is the state for each species?
    for r = 1:reps
        Asort = Asort_r[r]
        tlsort = tlsort_r[r]
        nichesort = nichesort_r[r]
        s = s_r[r] #NOTE: Not doing anything with this

        propneg_sp = Vector{Float64}(undef,size(Asort)[1]);
        for i = 1:size(Asort)[1]
            propneg_sp[i] = sum(s[end-100:end,i] .== -1) / 101
        end

        avg_durations = compute_average_durations(compute_state_durations(s[end-100:end,:]))
        avg_dur_pos1_r[r] = avg_durations[!,:pos1]
        avg_dur_neg1_r[r] = avg_durations[!,:neg1]

        propneg_r[r] = propneg_sp
        niche_r[r] = nichesort
        trophic_r[r] = tlsort

    end

    propneg_flattened = reduce(vcat, propneg_r);
    niche_flattened = reduce(vcat, niche_r);
    trophic_flattened = reduce(vcat, trophic_r);

    avg_dur_pos1_flattened = reduce(vcat, avg_dur_pos1_r);
    avg_dur_neg1_flattened = reduce(vcat, avg_dur_neg1_r);

    # Sort by niche value
    sorted_idx = sortperm(niche_flattened)
    niche_sorted = niche_flattened[sorted_idx]
    propneg_sorted = propneg_flattened[sorted_idx]
    trophic_sorted = trophic_flattened[sorted_idx]
    avg_dur_pos1_sorted = avg_dur_pos1_flattened[sorted_idx]
    avg_dur_neg1_sorted = avg_dur_neg1_flattened[sorted_idx]

    # Define a window width (e.g., 0.05 means ±0.025 around each point)
    window_width = 0.05

    # Prepare arrays for the sliding mean results
    mean_propneg = similar(niche_sorted, Float64)
    mean_trophic = similar(niche_sorted, Float64)
    mean_dur_pos1 = similar(niche_sorted, Float64)
    mean_dur_neg1 = similar(niche_sorted, Float64)
    std_propneg = similar(niche_sorted, Float64)

    # Calculate sliding averages using a similar approach
    @threads for i in eachindex(niche_sorted)
        center = niche_sorted[i]
        low_bound = center - window_width / 2
        high_bound = center + window_width / 2

        # Find the indices of points in the window using binary search
        low_idx = searchsortedfirst(niche_sorted, low_bound)
        high_idx = searchsortedlast(niche_sorted, high_bound)

        # Compute mean values within this window
        mean_propneg[i] = mean(propneg_sorted[low_idx:high_idx])
        std_propneg[i] = std(propneg_sorted[low_idx:high_idx])

        mean_trophic[i] = mean(trophic_sorted[low_idx:high_idx])
        mean_dur_pos1[i] = mean(filter(!isnan,avg_dur_pos1_sorted[low_idx:high_idx]))
        mean_dur_neg1[i] = mean(filter(!isnan,avg_dur_neg1_sorted[low_idx:high_idx]))
    end

    append!(trophic_all, mean_trophic)
    append!(n_all, niche_sorted)
    append!(C_all, fill(C, length(niche_sorted)))
    append!(propneg_all, mean_propneg)
    append!(propneg_std_all, std_propneg)
    append!(dur_pos1_all, mean_dur_pos1)
    append!(dur_neg1_all, mean_dur_neg1)

end

# filename = smartpath(string("data/duration_trophiclevel.jld2"))
# @save filename kin kout gprim S C_values tmax sigma reps trophic_all n_all C_all propneg_all propneg_std_all dur_pos1_all dur_neg1_all

filename = smartpath(string("data/duration_trophiclevel.jld2"))
@load filename kin kout gprim S C_values tmax sigma reps trophic_all n_all C_all propneg_all propneg_std_all dur_pos1_all dur_neg1_all


using LaTeXStrings
using PGFPlotsX
pgfplotsx()
cvalue = findall(x->x==0.02,C_all);
meanpos = moving_average(trophic_all[cvalue],dur_pos1_all[cvalue] ./ tmax,0.04);
meanneg = moving_average(trophic_all[cvalue],dur_neg1_all[cvalue] ./ tmax,0.04);
durplotS = Plots.plot(
    # n_all[cvalue],
    # trophic_clow,
    # trophic_all[cvalue],
    # dur_pos1_all[cvalue] ./ tmax,
    meanpos[1],
    meanpos[2],
    # title = "Positive State Durations",
    xlabel = L"\mathrm{Trophic~level},~\tau",
    ylabel = L"\mathrm{Within-cycle~sequential~duration}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    yscale = :log10,  # Use a logarithmic scale for the y-axis,
    label = L"S^+",
    legendfontsize=14,
    framestyle = :box,
    foreground_color_legend = nothing,
    linewidth=2
);
Plots.plot!(durplotS,
    # trophic_clow,
    # trophic_all[cvalue],
    # dur_neg1_all[cvalue] ./ tmax,
    meanneg[1],
    meanneg[2],
    label = L"S^-",
    legendfontsize=14,
    linewidth=2
);
durplotS

filename = smartpath("figures_prefinal/fig_duration_trophic2.pdf")
Plots.savefig(durplotS,filename)



########################################################
######### Generate sim data for Figure 5 ###############
########################################################

#Run for a specific implementation and for C = 0.02 and C= 0.1
kin = 5.;
kout = 10.;
gprim = 10.; #try 20

S=100;
C_values = [0.02,0.04,0.1] # or a range like 0.01:0.01:0.2
tmax = 300
sigma = 0.01
reps = 5000

n_all = Float64[]
C_all = Float64[]
propneg_all = Float64[]
propneg_std_all = Float64[]
dur_pos1_all = Float64[]
dur_neg1_all = Float64[]
trophic_all = Float64[]

for C in C_values

    # Prepare to collect outputs from all repetitions
    Asort_r = Vector{Matrix{Int64}}(undef, reps);
    tlsort_r = Vector{Vector{Float64}}(undef, reps);
    nichesort_r = Vector{Vector{Float64}}(undef, reps);
    s_r = Vector{Matrix{Int64}}(undef, reps);

    @threads for r = 1:reps

        # Generate network
        A = nothing  # Placeholder for adjacency matrix
        tl = Float64[]                    # Placeholder for trophic levels
        mintl = 0.0                       # Initial minimum trophic level
        niche = Float64[]

        while mintl < 1.0
            A, niche = nichemodelweb(S, C)
            tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6)
            mintl = minimum(tl)
        end
        # g = Graphs.DiGraph(A);

        tl = round.(TrophInd(Array{Float64}(A))[!, :TL]; sigdigits=6)

        tlsp = sortperm(tl)
        tlsort = tl[tlsp]
        Asort = A[tlsp, tlsp]
        nichesort = niche[tlsp]

        # Run simulation
        s = cascade(Asort, kout, kin, gprim, sigma, tmax)

        # Collect outputs
        Asort_r[r] = Asort
        tlsort_r[r] = tlsort
        nichesort_r[r] = nichesort
        s_r[r] = s

    end

    propneg_r = Array{Vector{Float64}}(undef,reps);
    niche_r = Array{Vector{Float64}}(undef,reps);
    trophic_r = Array{Vector{Float64}}(undef,reps);
    avg_dur_pos1_r = Array{Vector{Float64}}(undef,reps);
    avg_dur_neg1_r = Array{Vector{Float64}}(undef,reps);

    #What is the state for each species?
    for r = 1:reps
        Asort = Asort_r[r]
        tlsort = tlsort_r[r]
        nichesort = nichesort_r[r]
        s = s_r[r] #NOTE: Not doing anything with this

        propneg_sp = Vector{Float64}(undef,size(Asort)[1]);
        for i = 1:size(Asort)[1]
            propneg_sp[i] = sum(s[end-100:end,i] .== -1) / 101
        end

        avg_durations = compute_average_durations(compute_state_durations(s[end-100:end,:]))
        avg_dur_pos1_r[r] = avg_durations[!,:pos1]
        avg_dur_neg1_r[r] = avg_durations[!,:neg1]

        propneg_r[r] = propneg_sp
        niche_r[r] = nichesort
        trophic_r[r] = tlsort

    end

    propneg_flattened = reduce(vcat, propneg_r);
    niche_flattened = reduce(vcat, niche_r);
    trophic_flattened = reduce(vcat, trophic_r);

    avg_dur_pos1_flattened = reduce(vcat, avg_dur_pos1_r);
    avg_dur_neg1_flattened = reduce(vcat, avg_dur_neg1_r);

    #test 02/20/25
    # Moving average with non-overlapping window
    # avgtrophic,avgpropneg = moving_average(trophic_flattened,propneg_flattened,0.05);

    # Moving average with overlapping window
    # Sort by niche value
    sorted_idx = sortperm(niche_flattened)
    niche_sorted = niche_flattened[sorted_idx]
    propneg_sorted = propneg_flattened[sorted_idx]
    trophic_sorted = trophic_flattened[sorted_idx]
    avg_dur_pos1_sorted = avg_dur_pos1_flattened[sorted_idx]
    avg_dur_neg1_sorted = avg_dur_neg1_flattened[sorted_idx]

    # Define a window width (e.g., 0.05 means ±0.025 around each point)
    window_width = 0.05

    # Prepare arrays for the sliding mean results
    mean_propneg = similar(niche_sorted, Float64)
    mean_trophic = similar(niche_sorted, Float64)
    mean_dur_pos1 = similar(niche_sorted, Float64)
    mean_dur_neg1 = similar(niche_sorted, Float64)
    std_propneg = similar(niche_sorted, Float64)

    # Calculate sliding averages using a similar approach
    @threads for i in eachindex(niche_sorted)
        center = niche_sorted[i]
        low_bound = center - window_width / 2
        high_bound = center + window_width / 2

        # Find the indices of points in the window using binary search
        low_idx = searchsortedfirst(niche_sorted, low_bound)
        high_idx = searchsortedlast(niche_sorted, high_bound)

        # Compute mean values within this window
        mean_propneg[i] = mean(propneg_sorted[low_idx:high_idx])
        std_propneg[i] = std(propneg_sorted[low_idx:high_idx])

        mean_trophic[i] = mean(trophic_sorted[low_idx:high_idx])
        mean_dur_pos1[i] = mean(filter(!isnan,avg_dur_pos1_sorted[low_idx:high_idx]))
        mean_dur_neg1[i] = mean(filter(!isnan,avg_dur_neg1_sorted[low_idx:high_idx]))
    end

    append!(trophic_all, mean_trophic)
    append!(n_all, niche_sorted)
    append!(C_all, fill(C, length(niche_sorted)))
    append!(propneg_all, mean_propneg)
    append!(propneg_std_all, std_propneg)
    append!(dur_pos1_all, mean_dur_pos1)
    append!(dur_neg1_all, mean_dur_neg1)

end


#Export to fit with mathematica
df = DataFrame(n = n_all, C = C_all, propneg = propneg_all, durpos1 = dur_pos1_all, durneg1 = dur_neg1_all,trophic = trophic_all,propneg_sd = propneg_std_all);

gfilevalue = Int64(round(gprim));
koutfilevalue = Int64(round(kout));

file = string("data/n_c_propneg_data_","g",gfilevalue,"_kout",koutfilevalue,".csv")
filename = smartpath(file)
CSV.write(filename, df)




########################################################
######### GENERATE SIM DATA FOR FIGURE 6 ###############
########################################################

reps = 500

Cvec = [0.02,0.03,0.04];
lCvec = length(Cvec);

S = 100; 
gprim = 10
kout = 10
kin = 5.0
tmax = 300
sigma = 0.01

motifs = motifs_generate();

cascade_ratio = Array{Array{Float64}}(undef,lCvec,reps);

cascade_tl = Array{Array{Float64}}(undef,lCvec,reps);

tlreps = Array{Array{Float64}}(undef,lCvec,reps);

# Loop over parameter combinations
@showprogress 1 "Computing..." for c in 1:lCvec

    C = Cvec[c];

    # Prepare to collect outputs from all repetitions
    # Asort_list = Vector{Matrix{Int64}}(undef, reps)
    # tlsort_list = Vector{Vector{Float64}}(undef, reps)
    # s_list = Vector{Matrix{Int64}}(undef, reps)

    cascade_ratio_motif = Array{Float64}(undef,0)
    cascade_effect_tl = Array{Float64}(undef,0)

    for r in 1:reps
        
        # Generate network
        A = nothing  # Placeholder for adjacency matrix
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
        g = Graphs.DiGraph(Asort);

        # Run simulation
        s = cascade(Asort, kout, kin, gprim, sigma, tmax)

        #Calculate cascade effect across pairwise interactions
        
        sprime = s[tmax-100:tmax,:];

        #find the 2-chain interactions
        motiflabel = :m1; #pred-prey interactions
        emb_motif = find_motifs(g, motifs[motiflabel])
        lmotif = length(emb_motif)
        
        for k = 1:lmotif

            #Species in the motif
            spmotif = emb_motif[k]

            #NOTE: FIX!
            mintl_motif, sp_prey_pos = findmin(tlsort[spmotif])
            maxtl_motif, sp_pred_pos = findmax(tlsort[spmotif])

            sp_prey = spmotif[sp_prey_pos]
            sp_pred = spmotif[sp_pred_pos]
            
            offset = 1;
            coeff = 2;
                            
            # +P positions
            posP = findall(x->x==1,sprime[1:(end-offset),sp_pred])
            
            # -P positions
            negP = setdiff(collect(1:1:length(posP)),posP)
            
            if length(posP) > 0 && length(negP) > 0
                
                avg_state_posP = mean(sprime[posP .+ offset,sp_prey])
            
                avg_state_negP = mean(sprime[negP .+ offset,sp_prey])
                
                cascade_ratio_value = log((coeff + avg_state_posP) / (coeff + avg_state_negP))

                push!(cascade_ratio_motif,cascade_ratio_value)

                push!(cascade_effect_tl,tlsort[sp_prey])


            else

                push!(cascade_ratio_motif,NaN)
                push!(cascade_effect_tl,NaN)

            end


            # propneg_cascade_motif[k] = sum(filter(!isnan,cascade_ratio_motif) .< 0) / predpreypairs
            # meanratio_cascade_motif[k] = mean(filter(!isnan,cascade_ratio_motif))

        end #end motifs

        cascade_ratio[c,r] = filter(!isnan,cascade_ratio_motif);
        cascade_tl[c,r] = filter(!isnan,cascade_effect_tl);

        tlreps[c,r] = tlsort;

    end #end reps

    # StatsPlots.histogram(filter(!isnan,cascade_ratio_motif))

end

tlcascademean = moving_average(vcat(cascade_tl[1,:]...),vcat(cascade_ratio[1,:]...),0.01);
tlcascadesim = DataFrame(tl=tlcascademean[1],cascade=tlcascademean[2])
filename = smartpath("data/tlcascadeC02.csv")
CSV.write(filename, tlcascadesim)

tlcascademean = moving_average(vcat(cascade_tl[3,:]...),vcat(cascade_ratio[3,:]...),0.01);
tlcascadesim = DataFrame(tl=tlcascademean[1],cascade=tlcascademean[2])
filename = smartpath("data/tlcascadeC04.csv")
CSV.write(filename, tlcascadesim)




##################################
######### FIGURE 7 ###############
##################################


reps = 200
gvec = collect(0:1:25)
lgvec = length(gvec)
koutvec = collect(0:0.5:15)
lkoutvec = length(koutvec)

Cvec = [0.02,0.03,0.04];
lCvec = length(Cvec);

S = 100; 
kin = 5.0
tmax = 300
sigma = 0.01

motifs = motifs_generate();

mean_cascade_ratio = Array{Float64}(undef,lCvec,lgvec,lkoutvec);

# Loop over parameter combinations
@showprogress 1 "Computing..." for c in 1:lCvec

    C = Cvec[c];

    for i in 1:lgvec

        @threads for j in 1:lkoutvec

            gprim = gvec[i]
            kout = koutvec[j]

            # Prepare to collect outputs from all repetitions
            # Asort_list = Vector{Matrix{Int64}}(undef, reps)
            # tlsort_list = Vector{Vector{Float64}}(undef, reps)
            # s_list = Vector{Matrix{Int64}}(undef, reps)

            cascade_ratio_motif = Array{Float64}(undef,0)

            for r in 1:reps
                
                # Generate network
                A = nothing  # Placeholder for adjacency matrix
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
                g = Graphs.DiGraph(Asort);

                # Run simulation
                s = cascade(Asort, kout, kin, gprim, sigma, tmax)

                # Collect outputs
                # Asort_list[r] = Asort
                # tlsort_list[r] = tlsort
                # s_list[r] = s

                #Calculate cascade effect across pairwise interactions
                
                sprime = s[tmax-100:tmax,:];

                #find the 2-chain interactions
                motiflabel = :m1; #pred-prey interactions
                emb_motif = find_motifs(g, motifs[motiflabel])
                lmotif = length(emb_motif)
                
                # mintl_motif = Array{Float64}(undef,lmotif);
                # maxtl_motif = Array{Float64}(undef,lmotif);

                for k = 1:lmotif
                    #Species in the motif
                    spmotif = emb_motif[k]

                    #NOTE: FIX!
                    mintl_motif, sp_prey_pos = findmin(tlsort[spmotif])
                    maxtl_motif, sp_pred_pos = findmax(tlsort[spmotif])

                    sp_prey = spmotif[sp_prey_pos]
                    sp_pred = spmotif[sp_pred_pos]
                    
                    # #List each species' predators present in the current motif
                    # sp_preds = [intersect(spmotif,findall(x->x==1,Asort[spmotif[l],:])) for l=1:length(spmotif)]
                    # cascade_ratio_motif = Array{Float64}(undef,0)
                    
                    offset = 1;
                    coeff = 2;

                    # cascade_ratio_perpred = Array{Float64}(undef,length(preds))

                    # +P positions
                    posP = findall(x->x==1,sprime[1:(end-offset),sp_pred])
                    # -P positions
                    negP = setdiff(collect(1:1:length(posP)),posP)
                    
                    if length(posP) > 0 && length(negP) > 0
                        avg_state_posP = mean(sprime[posP .+ offset,sp_prey])
                        avg_state_negP = mean(sprime[negP .+ offset,sp_prey])
                        cascade_ratio_value = log((coeff + avg_state_posP) / (coeff + avg_state_negP))
                        push!(cascade_ratio_motif,cascade_ratio_value)
                    else
                        push!(cascade_ratio_motif,NaN)
                    end


                    # propneg_cascade_motif[k] = sum(filter(!isnan,cascade_ratio_motif) .< 0) / predpreypairs
                    # meanratio_cascade_motif[k] = mean(filter(!isnan,cascade_ratio_motif))

                end #end motifs
            end #end reps

            # StatsPlots.histogram(filter(!isnan,cascade_ratio_motif))
            mean_cascade_ratio[c,i,j] = mean(filter(!isnan,cascade_ratio_motif))

        end #end koutvec
        # println(i/lgvec)
    end #end gvec
end

# After all reps, save data to a single file
filename = smartpath("data/cascade_connectance.jld2")
@save filename reps gvec lgvec koutvec lkoutvec Cvec lCvec S kin tmax sigma mean_cascade_ratio


# After all reps, load data from a single file
filename = smartpath("data/cascade_connectance.jld2")
@load filename reps gvec lgvec koutvec lkoutvec Cvec lCvec S kin tmax sigma mean_cascade_ratio


p1 = Plots.scatter(vec(mean_cascade_ratio[1,:,:]),vec(mean_cascade_ratio[2,:,:]),
    xlims=(-0.9,0.0),
    ylims=(-0.9,0.0),
    xlabel="Cascade effect; C=0.02",
    ylabel="Cascade effect; C=Calt",
    alpha=0.3,
    markerstrokewidth=0,
    label="Calt=0.03",
    frame=:box);
Plots.scatter!(p1,vec(mean_cascade_ratio[1,:,:]),vec(mean_cascade_ratio[3,:,:]),
    alpha=0.3,
    markerstrokewidth=0,
    label="Calt=0.04");
Plots.plot!(p1,collect(-1:0.1:1),collect(-1:0.1:1),label=false,color=:black);
# Plots.scatter!(p1,
#     [mean_cascade_ratio[1,11,12]], [mean_cascade_ratio[2,11,12]],
#     color=:gray, label=false);
Plots.scatter!(p1,
    [mean_cascade_ratio[1,11,12]], [mean_cascade_ratio[3,11,12]],
    color=:black, label=false);

filename = smartpath("figures_prefinal/fig_cascade_effect_conn.pdf")
Plots.savefig(p1,filename)


ratio = 1 ./ (mean_cascade_ratio[3,:,:] ./ mean_cascade_ratio[1,:,:]);

ratioplot = Plots.heatmap(gvec,koutvec,ratio',
    clims=(0.5,1.5),
    xlabel=L"\mathrm{Productivity,~}g",
    ylabel=L"\mathrm{Top-down~magnitude,~}k^{\rm out}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_titlefontsize = 14,
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Cascade~effect~ratio,~}\bar{\psi}_{\rm high}/\bar{\psi}_{\rm low}");

filename = smartpath("figures_prefinal/fig_cascade_effect_ratio.pdf");
Plots.savefig(ratioplot,filename)

using Measures
combinedplot = plot(p1, ratioplot, 
    layout = (1, 2), 
    size = (1000, 400),
    left_margin=10mm,  # Add more space to the left
    bottom_margin=5mm);
filename = smartpath("figures_prefinal/fig_cascade_effect_connratio.pdf");
Plots.savefig(combinedplot,filename)














# NEW FIGURES/CODE - 05/27/2025

########################################################
######### Appendix 1: Lotka-Volterra Map ###############
########################################################
using Statistics, Printf, Plots, LaTeXStrings

# Parameterise a *stable* LV oscillator (daily units)

# August 15 2025
r       = [ 0.3,  0.0 ]                 # r_1>0 prey, r_2 used as mortality d = -r_2 -0.01
a_pos   = [ 0.0   0.0 ;                 # a_21 > 0 : prey → predator conversion
            0.09  0.0 ]
b_neg   = [ 0.02 0.19 ;                # b_11 self-lim prey  ; b_12 predation
            0.0   0.02 ]               # b_22 self-lim predator

Δt      = 1                         # integration step (days)
steps   = 1000

function lv_euler!(X, r, a, b, Δt, steps)
    n = length(X)
    out = zeros(n, steps)
    out[:,1] .= X
    for t in 1:steps-1
        Xi = out[:,t]
        gain =  r .* Xi .+  sum(a .* (Xi .* Xi'), dims=2)
        loss = sum(b .* (Xi .* Xi'), dims=2)
        out[:,t+1] .= Xi .+ Δt .* (gain .- loss)   # Euler step
    end
    return out
end

X0    = [1.0, 1.0]
Xtraj = lv_euler!(X0, r, a_pos, b_neg, Δt, steps)

plot(Xtraj[1,:])
plot!(Xtraj[2,:])

############################################################
# Compute one reference biomass per species
############################################################
X̂_low  = [quantile(Xtraj[1,:],0.10), quantile(Xtraj[2,:],0.10)]
X̂_high = [quantile(Xtraj[1,:],0.90), quantile(Xtraj[2,:],0.90)]

X̂ = 0.5 .* (X̂_low .+ X̂_high)          # <— single Boolean threshold

############################################################
# Derive C, σ, k and Boolean update
############################################################
c  =  a_pos .- b_neg
ciidi = diag(c)
Q  =  r .+ c*X̂
# σ  =  1 .+ Q #The approximation without self-limitation (works just as well)
σ  = 1 .+ Δt .* (Q .+ ciidi .* X̂)
k  =  c .* X̂'          # k_{ij}=c_{ij}·X̂_j
g  =  Q ./ σ

# time-step factors: multiply the non-self terms by Δt
A  = 1 .+ Δt .* Q              # diagonal coeffs  (1 + Δt C_i)
B  =  Δt .* Q                  # constant term
K  =  Δt .* k                  # interaction matrix
K_norm = K ./ σ                 # element-wise row division: k_{ij}/σ_i

############################################################
# Simulate the linearised map  s_lin(t) - note NOT σ-normalised
############################################################
s_lv = (Xtraj .- X̂) ./ X̂      # “ground-truth’’ deviations from LV run

s_lin       = zeros(2, steps)
s_lin[:,1] .= s_lv[:,1]        # start from same initial deviation
for t in 1:steps-1
    s_lin[1,t+1] = A[1]*s_lin[1,t] + B[1] + K[1,2]*s_lin[2,t]
    s_lin[2,t+1] = A[2]*s_lin[2,t] + B[2] + K[2,1]*s_lin[1,t]
end

# Simulate the σ-normalised linear map
#     s_i(t+1) = s_i(t) + g_i + Σ_j (k_{ij}/σ_i) s_j(t)
############################################################
s_norm        = zeros(2, steps)
s_norm[:,1]  .= s_lv[:,1] ./ σ      # scaled deviation of LV
for t in 1:steps-1
    # explicit form 
    s_norm[1,t+1] = s_norm[1,t] + g[1] + K_norm[1,2]*s_norm[2,t]
    s_norm[2,t+1] = s_norm[2,t] + g[2] + K_norm[2,1]*s_norm[1,t]
end

############################################################
# Boolean map
############################################################
S_lv = sign.(s_lv)
S_lin  = sign.(s_lin)    
S_norm = sign.(s_norm)

tmax = 250; #copy(steps)

#THIS IS THE PLOT
pltA = plot(s_lv[1,1:tmax],
    label="Simulated",     
    lw=2,  
    color=:blue,
    xlabel=L"Time step, $t$", 
    ylabel=L"Relative dynamic, $s(t)$",
    frame=:box,
    foreground_color_legend = nothing
    )
# plot!(pltA,s_lin[1,1:tmax],
#     label="Linearized map",     
#     lw=2,  
#     color=:red)
plot!(pltA,s_norm[1,1:tmax],
    label="Linearized map",     
    lw=2,  
    color=:red)
# Boolean sequence implied by s_lin

pltB = plot(S_lv[1,1:tmax],
    ylim=[-1.5,1.5],
    label="Simulated",     
    lw=2,  
    color=:blue,
    xlabel=L"Time step, $t$", 
    ylabel=L"Boolean state dynamic, $S(t)$",
    frame=:box,
    foreground_color_legend = nothing
    )
# plot!(pltB,S_lin[1,1:tmax],
#     label="Linearized map",
#     color=:red,
#     lw=2)
plot!(pltB,S_norm[1,1:tmax],
    label="Linearized map",
    color=:red,
    lw=2)

using Measures
combLVplot = plot(pltA, pltB, layout=(1,2), 
    size = (1000, 400),
    left_margin=10mm,  # Add more space to the left
    bottom_margin=5mm);
display(combLVplot)
filename = smartpath("figures_prefinal/suppfig_LV.pdf");
Plots.savefig(combLVplot,filename)




###########################################################
######### Sensitivity to species-specific k ###############
###########################################################
ksdvec = [0.0,0.5,1.0,1.5,2.0,5.0,10.0]
ksdstrs = [ @sprintf("%02d", Int(round(10*x))) for x in ksdvec ]
subfile = "";

for k in 1:length(ksdvec)

    ksd = ksdvec[k];
    ksdname = ksdstrs[k];

    Cnorm = 0.02;
    reps = 200
    gvec = collect(0:1:25)
    lgvec = length(gvec)
    koutvec = collect(0:0.5:15)
    lkoutvec = length(koutvec)

    S = 100; 
    kinmean = 5.0
    tmax = 300
    sigma = 0.01

    C = Cnorm;

    # Optional: Save preamble or settings if needed
    filename_preamble = smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout",subfile,"/preamble.jld2"))

    @save filename_preamble reps gvec koutvec S C kinmean tmax sigma

    # Loop over parameter combinations
    # @showprogress 1 "Computing..." 
    for i in 1:lgvec
        @threads for j in 1:lkoutvec
            gprim = gvec[i]
            koutmean = koutvec[j]

            # Prepare to collect outputs from all repetitions
            Asort_list = Vector{Matrix{Int64}}(undef, reps)
            tlsort_list = Vector{Vector{Float64}}(undef, reps)
            s_list = Vector{Matrix{Int64}}(undef, reps)

            # Use threading over repetitions
            for r in 1:reps
                
                # Generate network
                A = nothing  # Placeholder for adjacency matrix
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

                # Run simulation
                s = cascade_kdiverse(Asort, koutmean, kinmean, ksd, gprim, sigma, tmax)

                # Collect outputs
                Asort_list[r] = Asort
                tlsort_list[r] = tlsort
                s_list[r] = s
            end

            # After all reps, save data to a single file

            filename = smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout",subfile,"/g$(i)_kout$(j).jld2"))
            @save filename Asort_list tlsort_list s_list gprim koutmean kinmean ksd sigma tmax
        end
        println(i/lgvec)
    end
end

# CALCULATE STATUSICS FROM SIMULATIONS

subfile = "";
for k=1:length(ksdvec)
    ksd = ksdvec[k];
    ksdname = ksdstrs[k];
    # Load variables from saved file
    filename_preamble = smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout",subfile,"/preamble.jld2"))
    @load filename_preamble reps gvec koutvec S C kinmean tmax sigma

    lgvec = length(gvec);
    lkoutvec = length(koutvec);
    # Initialize 3D arrays
    prop_changes = Array{Float64}(undef, lgvec, lkoutvec, reps);
    prop_oscillating = Array{Float64}(undef, lgvec, lkoutvec, reps);

    avg_time_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    avg_time_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    avg_dur_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    avg_dur_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    std_time_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    std_time_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    std_dur_pos1 = Array{Float64}(undef, lgvec, lkoutvec, reps);
    std_dur_neg1 = Array{Float64}(undef, lgvec, lkoutvec, reps);

    transition_pospos = Array{Float64}(undef, lgvec, lkoutvec, reps);
    transition_posneg = Array{Float64}(undef, lgvec, lkoutvec, reps);
    transition_negpos = Array{Float64}(undef, lgvec, lkoutvec, reps);
    transition_negneg = Array{Float64}(undef, lgvec, lkoutvec, reps);

    cascade_ratio_offset = Array{Float64}(undef, lgvec, lkoutvec, reps);
    cascade_ratio = Array{Float64}(undef, lgvec, lkoutvec, reps);

    propneg_effect_ratio = Array{Float64}(undef, lgvec, lkoutvec, reps);
    propneg_effect_ratio_offset = Array{Float64}(undef, lgvec, lkoutvec, reps);

    propneg_effect_interactions = Array{Float64}(undef, lgvec, lkoutvec, reps);
    propneg_effect_interactions_offset = Array{Float64}(undef, lgvec, lkoutvec, reps);

    @showprogress 1 "Computing..." for i in 1:lgvec
        @threads for j in 1:lkoutvec
            filename = smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout",subfile,"/g$(i)_kout$(j).jld2"))
            @load filename Asort_list tlsort_list s_list gprim koutmean kinmean ksd sigma tmax
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

            end
        end
    end

    # Compute mean across dimension 3 while ignoring NaNs
    mean_prop_oscillating = mapslices(x -> mean(filter(!isnan, x)), prop_oscillating; dims=3)[:,:];
    mean_cascade_ratio = mapslices(x -> mean(filter(!isnan, x)), cascade_ratio; dims=3)[:,:];
    mean_cascade_ratio_offset = mapslices(x -> mean(filter(!isnan, x)), cascade_ratio_offset; dims=3)[:,:];

    #### Save output
    # newdir = replace(smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout",subfile,".jl")), ".jl" => "/");
    # mkpath(newdir);
    filename = smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout",subfile,"/post_analysis.jld2"))

    @save filename mean_prop_oscillating avg_dur_pos1 avg_dur_neg1 transition_pospos transition_posneg transition_negpos transition_negneg avg_time_pos1 avg_time_neg1 std_time_pos1 std_time_neg1 std_dur_pos1 std_dur_neg1

end


# Create a dictionary to store results for each run
results = Dict()
for k in 1:length(ksdvec)
    ksdname = ksdstrs[k];
    # Load the file
    filename = smartpath(string("data/kdiverse_sd$(ksdname)_gvec_kout", subfile, "/post_analysis.jld2"))
    # Initialize a sub-dictionary for this specific run
    run_data = Dict()
    @load filename mean_prop_oscillating avg_dur_pos1 avg_dur_neg1 transition_pospos transition_posneg transition_negpos transition_negneg avg_time_pos1 avg_time_neg1 std_time_pos1 std_time_neg1 std_dur_pos1 std_dur_neg1
    run_data[:mean_prop_oscillating] = mean_prop_oscillating
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
    results[k] = run_data
end


using LaTeXStrings
using PGFPlotsX
pgfplotsx()
plotoscillating = Plots.heatmap(results[1][:mean_prop_oscillating]',
    xlabel=L"\mathrm{Productivity,}~g",
    ylabel=L"\mathrm{Top-down~magnitude,}~k^{\mathrm{out}}",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    colorbar_tickfontsize = 12,
    colorbar_title = L"\mathrm{Proportion~cycling}",
    colorbar_titlefontsize = 14,
    color=:thermal);
display(plotoscillating)


baseline = results[1][:mean_prop_oscillating]';
summederror = Array{Float64}(undef,length(ksdvec));
for k=1:length(ksdvec)
    sd_results = results[k][:mean_prop_oscillating]';
    summederror[k] = sum(sqrt.((baseline .- sd_results).^2)) / length(baseline);
end

pltsd = scatter(ksdvec,summederror,
    xlab = L"k standard deviation, $\sigma_k$",
    ylab = "Mean absolute error",
    frame = :box,
    label = false);
display(pltsd)
filename = smartpath("figures_prefinal/fig_sderror.pdf");
Plots.savefig(pltsd,filename)



#Trophic level vs. Connectance
using Graphs          # works with Graphs v2 or LightGraphs ≤v1.3
# ------------------------------------------------------------------
#  Breadth-first search that returns a vector of shortest‐path lengths
# ------------------------------------------------------------------
function bfs_distances(g::Graphs.SimpleDiGraph, src::Int)
    n     = Graphs.nv(g)
    dist  = fill(-1, n)         # -1  ≡  “unreached”
    dist[src] = 0
    queue = [src]

    while !isempty(queue)
        v = popfirst!(queue)
        for w in Graphs.outneighbors(g, v)   # prey → predator direction
            dist[w] == -1 || continue  # already discovered
            dist[w] = dist[v] + 1
            push!(queue, w)
        end
    end
    return dist                 # shortest paths from src (in # links)
end

# ------------------------------------------------------------------
#  Minimal food-chain height  (= longest *shortest* basal→top path)
# ------------------------------------------------------------------
function chain_height(A::AbstractMatrix{Bool})
    g = Graphs.SimpleDiGraph(A)                   # edge: prey  →  predator

    # Basals are columns with no in-links (no prey of their own)
    basals = findall(i -> sum(@view A[:, i]) == 0, 1:Graphs.nv(g))
    maxlen = 0

    for b in basals
        dists = bfs_distances(g, b)        # shortest paths from basal b
        for d in dists
            d == -1 && continue            # unreachable nodes
            maxlen = max(maxlen, d)
        end
    end
    return maxlen                          # # links; add +1 for levels
end

S = 100;
reps = 500;
Cvec = collect(0.02:0.005:0.1)
maxtl = Array{Float64}(undef,length(Cvec),reps);
avgtl = Array{Float64}(undef,length(Cvec),reps);
chainlen = Array{Float64}(undef,length(Cvec),reps);

for c = 1:length(Cvec)
    C = Cvec[c];

    @threads for r=1:reps
        # Generate network
        A = nothing  # Placeholder for adjacency matrix
        tl = Float64[]                    # Placeholder for trophic levels
        mintl = 0.0                       # Initial minimum trophic level
        niche = Float64[]

        while mintl < 1.0
            A, niche = nichemodelweb(S, C)
            tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6)
            mintl = minimum(tl)
        end
        # g = Graphs.DiGraph(A);

        tl = round.(TrophInd(Array{Float64}(A))[!, :TL]; sigdigits=6)

        tlsp = sortperm(tl)
        tlsort = tl[tlsp]
        Asort = A[tlsp, tlsp]
        nichesort = niche[tlsp]

        maxtl[c,r] = maximum(tlsort)
        avgtl[c,r] = mean(tlsort)
        chainlen[c, r] = chain_height(A)
    end
end

exp_maxtl = mean(maxtl,dims=2)[:,1]
exp_avgtl = mean(avgtl,dims=2)[:,1]
exp_chainlen = mean(chainlen,dims=2)[:,1]


nicheTLvCplot = plot(Cvec,exp_maxtl,
    xlab = "Food web connectance, C",
    ylab = "Food web trophic level",
    label = "Maximum TL",
    width = 2,
    frame = :box
    );
plot!(nicheTLvCplot,Cvec,exp_avgtl,
    label = "Average TL",
    width = 2);
plot!(nicheTLvCplot,Cvec,exp_chainlen,
    label = "Chain length",
    width = 2);
filename = smartpath("figures_prefinal/fig_TLvC.pdf");
Plots.savefig(nicheTLvCplot,filename)
