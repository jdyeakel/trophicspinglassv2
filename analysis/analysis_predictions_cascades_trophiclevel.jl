
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


# function exp_cascaderatio_simple(tl,C,tldiff)

#     #Calculate expected cascade ratio
#     num = 2 + (pr_state_neg(tl,C)*(-1) + (1-pr_state_neg(tl,C)*(+1)))*(1-pr_state_neg(tl+tldiff.tdpred,C))

#     denom = 2 + (pr_state_neg(tl,C)*(-1) + (1-pr_state_neg(tl,C)*(+1)))*(pr_state_neg(tl+tldiff.tdpred,C))

#     expcascade = log(num/denom)

#     return expcascade

# end

# function exp_cascaderatio_structure(tl,C,tldiff)

#     nichevalue = tl_to_niche(tl,C)

#     pr_Nstate_neg_num =  kout*(1-nichevalue)*(1-2*pr_state_neg(tl+tldiff.tdprey,C))
    
#     pr_Nstate_neg_denom = kin*nichevalue*(1-2*pr_state_neg(tl-tldiff.tdprey,C)) + pr_Nstate_neg_num
    
#     pr_Nstate_neg = pr_Nstate_neg_num / pr_Nstate_neg_denom

#     # exp_state = (1 - 2*pr_state_neg(tl,C)) + 2*(S-1)*C * (kin*nichevalue*(1 - 2*pr_state_neg(tl-tldiff,C)) - kout*(1-nichevalue)*(1 - 2*pr_state_neg(tl+tldiff,C)))

#     nichevalue_pred = tl_to_niche(tl+tldiff.tdpred,C)

#     pr_predstate_neg_num =  kout*(1-nichevalue_pred)*(1-2*pr_state_neg(tl+2*tldiff.tdpred,C))
    
#     pr_predstate_neg_denom = kin*nichevalue_pred*(1-2*pr_state_neg(tl,C)) + pr_predstate_neg_num
    
#     pr_predstate_neg = pr_predstate_neg_num / pr_predstate_neg_denom

#     expcascade_num = 2 + ((-1)*pr_Nstate_neg + (1)*(1 - pr_Nstate_neg))*(1 - pr_predstate_neg)

#     expcascade_denom = 2 + ((-1)*pr_Nstate_neg + (1)*(1 - pr_Nstate_neg))*pr_predstate_neg

#     expcascade = log(expcascade_num / expcascade_denom)

#     return expcascade

# end

    

# This is the cascade effect expectation if we IGNORE neighbors and just use the probability of state alone

kout = 4.
kin = 10.
g = 10
#What is the trophic level difference of the predator compared to the prey?
tldiff_pred = 1
tldiff_prey = 1
tldiff = (tdprey=tldiff_prey,tdpred=tldiff_pred)

S = 100;
Cvalues = [0.02,0.1]
trophiclevels = collect(1:0.05:5);
# cascaderatio_simple = Array{Float64}(undef,length(Cvalues),length(trophiclevels));
# cascaderatio_structure = Array{Float64}(undef,length(Cvalues),length(trophiclevels));
cascaderatio_conditional = Array{Float64}(undef,length(Cvalues),length(trophiclevels));

for c in eachindex(Cvalues)
    
    C = Cvalues[c]
    
    for i in eachindex(trophiclevels)

        tl = trophiclevels[i]

        # cascaderatio_simple[c,i] = exp_cascaderatio_simple(tl,C,tldiff)
        # cascaderatio_structure[c,i] = exp_cascaderatio_structure(tl,C,tldiff)
        # println([C,i])
        cascaderatio_conditional[c,i] = exp_cascaderatio_conditional(tl,S,C,tldiff,kin,kout,g)

    end
end

ccplot = Plots.plot(trophiclevels,cascaderatio_conditional[1,:],
    ylims=(-1.5,0.01),
    ylabel="Cascade Index",
    xlabel="Trophic Level",
    linewidth=2,
    label="C=0.02",
    framestyle=:box); # Set the y-axis range)
Plots.plot!(ccplot,trophiclevels,cascaderatio_conditional[2,:],
    linewidth=2,
    label="C=0.10");
ccplot

# Plots.plot!(ccplot,maxtl_avg,cascadeindex_avg)
# Plots.scatter!(ccplot,maxtl_motif_flattened .- 1,meanratio_cascade_motif_flattened)

# csplot = Plots.plot(trophiclevels,cascaderatio_structure[1,:],ylims=(minimum(cascaderatio_structure)*(1+0.1), maximum(cascaderatio_structure)*(1+0.1))); # Set the y-axis range)
# Plots.plot!(csplot,trophiclevels,cascaderatio_structure[2,:]);
# csplot

# cplot = Plots.plot(trophiclevels,cascaderatio_simple[1,:],ylims=(minimum(cascaderatio_simple), maximum(cascaderatio_simple))); # Set the y-axis range)
# Plots.plot!(cplot,trophiclevels,cascaderatio_simple[2,:]);
# cplot

# See https://chatgpt.com/share/675b345b-f624-8009-b936-51943ecf54e8

