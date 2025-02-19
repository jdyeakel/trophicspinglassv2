
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

using Plots
using UnicodePlots
using StatsPlots


#Build 
# 
# S=100;
# C=0.01;
# @time A,n = nichemodelweb(S,C);
#R"""
#image($(A),col=grey(c(0,1)))
#"""

# A food chain example
S = 10;
A = motif_chain(S);

# The Haberman et al. example
A = motif_haberman();

# A food web model exzmple
S=100; C=0.02;
A,niche = nichemodelweb(S,C);
g = Graphs.DiGraph(A);
tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6);

# Plot the adjacency matrix
UnicodePlots.heatmap(A)

plotweb(A,tl)

# # Node names
# node_names = ["Primary", "Consumer", "Predator"]

# # Dead nodes (optional)
# Dead = ["Decomposer"]

@time result = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6);

# @time result_old = round.(trophic_old(A); sigdigits=6);





S=size(A)[1];
tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6); 
tlsp = sortperm(tl); 
tlsort = tl[tlsp];
Asort = A[tlsp,tlsp];
tmax=200;
#Predation coupling
kout=2.;
#Consumption coupling
kin=5.;
#Global influence of primary producers
gprim = 5.;
#Noise
sigma=0.01;
@time s = cascade(Asort,kout,kin,gprim,sigma,tmax);
# R"""
# image($(s),col=c('black','white'))
# """

UnicodePlots.heatmap(s')

Plots.heatmap(s')


plotweb(A,tl,s[tmax,:])


# Detect motifs in Haberman network
net_adj, _ = motif_haberman();
g = Graphs.DiGraph(net_adj);
plotweb(net_adj,tl)


motif_adj, _ = smallwebs(6)
motif = Graphs.DiGraph(motif_adj);

find_motifs(g, motif)

motifs = motifs_generate()

find_motifs(g, motifs[:m1])

# m1_g, m2_g, m3_g, m4_g, m5_g, m6_g, m7_g, m8_g, m9_g = motifs_generate()


plotweb(Asort,tlsort)

#Exploring dynamics

# Given s is the state matrix

#State change metric
sprime = s[end-100:end,:]
prop_spin = compute_proportion_oscillating(sprime)
n_changes = compute_state_changes(sprime)
avg_durations = compute_average_durations(compute_state_durations(sprime))

#Time metrics
change_times = compute_state_change_times(sprime)
total_time_in_states = compute_total_time_in_states(sprime)[!,:TotalTime_neg1]

#Correlation metrics
cospin_matrix = compute_cospin_matrix(sprime,0)
transition_probs = compute_state_transition_probabilities(sprime)

cascade_effect_ratio = cascade_effect(sprime,Asort,1);
mean(filter(!isnan,cascade_effect_ratio))


#Plot example sim
A = motif_haberman();
S = size(A)[1]
tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6); 
tlsp = sortperm(tl); 
tlsort = tl[tlsp];
Asort = A[tlsp,tlsp];
tmax=50;
#Predation coupling
kout=2.;
#Consumption coupling
kin=5.;
#Global influence of primary producers
gprim = 5.;
#Noise
sigma=0.01;
@time s = cascade(Asort,kout,kin,gprim,sigma,tmax);

Asort_bool = Matrix{Bool}(Bool.(Asort))
plotweb(Asort_bool,tlsort,s[tmax,:])

using PGFPlotsX
pgfplotsx()

simfig = Plots.heatmap(s',
    xlabel="Time",
    ylabel="Species sorted by trophic level",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    color=cgrad([:white,:gray]),
    frame=:box,
    legend=false)

filename = smartpath("figures/fig_sim.pdf")
Plots.savefig(simfig,filename)