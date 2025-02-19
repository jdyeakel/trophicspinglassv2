module trophicspinglass

using StatsBase
using RCall
using Graphs
using IterTools #Can remove
using Combinatorics
using GraphPlot
using Colors
using Random
using Distributions 
using DataFrames
using LinearAlgebra 
using JLD2 
# using LightGraphs
using Base.Threads

include("cascade.jl")
include("detectcascade.jl")
include("nichemodelweb.jl")
include("plotweb.jl")
include("quantitativeweb.jl")
include("motif_chain.jl")
include("trophic.jl")
include("smartpath.jl")
include("motif_haberman.jl")
include("find_motifs.jl")
include("is_isomorphic.jl")
include("motifs_generate.jl")
include("dynamics.jl")
include("moving_average.jl")
include("moving_sd.jl")
include("analytical_cascade_expressions.jl")

export 
cascade,
detectcascade,
nichemodelweb,
plotweb,
quantitativeweb,
motif_chain,
InternalNetwork,
Diet,
TrophInd,
# trophic
smartpath,
motif_haberman,
find_motifs,
is_isomorphic,
motifs_generate,
moving_average,
moving_sd,

compute_proportion_oscillating,
compute_state_changes,
compute_state_durations,
compute_average_durations,
compute_state_change_times,
compute_cospin_matrix,
compute_total_time_in_states,
compute_state_transition_probabilities,
cascade_effect,

tl_to_niche,
pr_state_neg,
prob_state,
exp_cascaderatio_conditional,

pr_state_negfit,
exp_cascaderatio_conditionalfit


end # module trophicspinglass
