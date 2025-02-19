
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
using GLM
using CSV


#Set parameterizations to equal that of analytical predictions

kout = 10.;
kin = 5.;
gprim = 10.;

S = 100;
Cvalues = [0.02]; #,0.1

tmax = 300;
sigma = 0.0;
reps = 20000;

cascade_ratio_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);
trophiclevel_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);
tldiffpred_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);
tldiffprey_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);

asymmetry_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);
probneg_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);
niche_sim = Array{Vector{Float64}}(undef,length(Cvalues),reps);


for c in eachindex(Cvalues)

    C = Cvalues[c]

    @threads for r=1:reps

        #Sample freely, as was done to create the probneg relationship
        # kout = rand(collect(0:0.1:15));
        # kin = rand(collect(0:0.1:15));
        # gprim = rand(collect(0:0.1:25));

        # Generate network
        A = nothing  # Placeholder for adjacency matrix
        niche = Float64[]
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

        # Run simulation
        s = cascade(Asort, kout, kin, gprim, sigma, tmax)

        sprime = s[end-100:end,:];

        #species with predators
        sp_all = collect(1:size(Asort)[1])
        # sp_withpreds = findall(!iszero,vec(sum(Asort,dims=2)))
        
        offset = 1;
        coeff = 2;
        
        cascade_ratio_A = Vector{Float64}(undef,0)
        trophiclevel_A = Vector{Float64}(undef,0)
        tldiffpred_A = Vector{Float64}(undef,0)
        tldiffprey_A = Vector{Float64}(undef,0)

        asymmetry_A = Vector{Float64}(undef,0)
        probneg_A = Vector{Float64}(undef,0)
        niche_A = Vector{Float64}(undef,0)
        
        #calculate the cascade index
        for i in eachindex(sp_all)
            
            # preyi = sp_withpreds[i]
            preyi = sp_all[i]
            trophiclevelsp = tlsort[preyi]

            i_res = findall(!iszero,Asort[:,preyi])
            i_preds = findall(!iszero,Asort[preyi,:])

            #save asymmetry relationship
            #Primary producers will be Inf - can skip those later
            asymmetry_i = length(i_preds) / length(i_res)
            #save probability negative relationship
            probneg_i = sum(sprime[:,preyi] .== -1) / 101;
            niche_i = nichesort[preyi]

            #Don't save Inf! (primary producers)
            # if asymmetry_i != Inf
            push!(asymmetry_A,asymmetry_i)
            push!(probneg_A,probneg_i)
            push!(niche_A,niche_i)
            push!(trophiclevel_A,trophiclevelsp)
            # end

            #For prey i, what is the average tldiff of its resources?
            
            # tlpreds = mean(tlsort[i_preds])
            
            if length(i_preds) > 0
                for j in eachindex(i_preds)
                    
                    predj = i_preds[j]

                    tldiff_pred = tlsort[predj] - tlsort[preyi]
                    tldiff_prey = tlsort[preyi] - mean(tlsort[i_res])
                    
                    # +P positions
                    posP = findall(x->x==1,sprime[1:(end-offset),predj])
                    
                    # -P positions
                    negP = setdiff(collect(1:1:length(posP)),posP)

                    cascade_ratio_perpred = 0.0

                    if length(posP) > 0 && length(negP) > 0
                        avg_state_posP = mean(sprime[posP .+ offset,preyi])
                        
                        avg_state_negP = mean(sprime[negP .+ offset,preyi])
                        
                        cascade_ratio_perpred = log((coeff + avg_state_posP) / (coeff + avg_state_negP))
                    else
                        cascade_ratio_perpred = NaN
                    end
                    
                    # push!(cascade_ratio_sim,cascade_ratio_perpred)
                    # push!(trophiclevel_sim,trophiclevelsp)

                    # #Don't save Inf! (primary producers)
                    # if asymmetry_i != Inf
                    #     push!(asymmetry_A,asymmetry_i)
                    #     push!(probneg_A,probneg_i)
                    #     push!(niche_A,niche_i)
                    # end

                    #Don't save NaN!
                    if cascade_ratio_perpred != NaN
                        push!(cascade_ratio_A,cascade_ratio_perpred)
                        push!(tldiffpred_A,tldiff_pred)
                        push!(tldiffprey_A,tldiff_prey)
                    end
                end #end j
            end
        end #end i
        
        # plotweb(Asort,tlsort)
        # UnicodePlots.scatterplot(trophiclevel_A,tldiffpred_A)
        # mean(tldiffpred_A)
        niche_sim[c,r] = niche_A
        probneg_sim[c,r] = probneg_A
        asymmetry_sim[c,r] = asymmetry_A
        cascade_ratio_sim[c,r] = cascade_ratio_A
        trophiclevel_sim[c,r] = trophiclevel_A
        tldiffpred_sim[c,r] = tldiffpred_A
        tldiffprey_sim[c,r] = tldiffprey_A
    end #end r
end #end c

UnicodePlots.scatterplot(log.(reduce(vcat, asymmetry_sim[1,:])),reduce(vcat, probneg_sim[1,:]))
Plots.scatter(log.(reduce(vcat, asymmetry_sim[1,:])),reduce(vcat, probneg_sim[1,:]))

Plots.scatter(moving_average(log.(reduce(vcat, asymmetry_sim[1,:])),reduce(vcat, probneg_sim[1,:]),0.1))

nichevalues = reduce(vcat, niche_sim[1,:]);
highniche_pos = findall(x->x>0.8,nichevalues)
asymmvalues = (reduce(vcat, asymmetry_sim[1,:])); #log
probnegvalues = reduce(vcat, probneg_sim[1,:]);
bheavy_pos = findall(x->x<0,asymmvalues);
theavy_pos = findall(x->x>0,asymmvalues);


#NOTE: This relationship is important
# asymmvalues = replace(asymmvalues, Inf => 0, -Inf => 0)

positions = findall(x -> x != Inf && x != -Inf, asymmvalues)
# asymmvalues[positions]
# highniche_pos = findall(x->x>0.0,nichevalues)
asymmavg = moving_average(log.(asymmvalues[positions]),probnegvalues[positions],0.1)
psym = Plots.scatter(asymmavg,
    ylims=(0.48,0.51),
    frame=:box,
    xlabel="Asymmetry index",
    ylabel="Probability S-",
    legend=false)
filename=smartpath("figures/fig_asymm_probneg.pdf")
Plots.savefig(psym,filename)


positions = findall(x -> x != Inf && x != -Inf, asymmvalues)
# asymmvalues[positions]
# highniche_pos = findall(x->x>0.0,nichevalues)
nicheavg = moving_average(nichevalues[positions],probnegvalues[positions],0.05)
psym = Plots.scatter(nicheavg,
    ylims=(0.48,0.51),
    frame=:box,
    xlabel="Niche value",
    ylabel="Probability S-",
    legend=false)



Plots.scatter(moving_average(asymmvalues[positions],nichevalues[positions],0.1))


# Separate probnegvalues into bottom-heavy and top-heavy groups
bheavy_probneg = probnegvalues[intersect(highniche_pos,bheavy_pos)]
theavy_probneg = probnegvalues[intersect(highniche_pos,theavy_pos)]
# Create labels for the groups
group_labels = vcat(fill("Bottom Heavy", length(bheavy_probneg)),
                    fill("Top Heavy", length(theavy_probneg)))
# Combine values into a single array for plotting
combined_values = vcat(bheavy_probneg, theavy_probneg)
# Plot the boxplot
Plots.boxplot(group_labels, combined_values,
    xlabel="Asymmetry Category",
    ylabel="Probneg Values",
    title="Comparison of Probneg Values by Asymmetry Category",
    legend=false)



UnicodePlots.scatterplot(reduce(vcat, trophiclevel_sim[1,:]),reduce(vcat, tldiffpred_sim[1,:]))

tl_avg, tldiffpred_avg = moving_average(reduce(vcat, trophiclevel_sim[1,:]),reduce(vcat, tldiffpred_sim[1,:]),0.2);

UnicodePlots.lineplot(tl_avg,tldiffpred_avg)
Plots.plot(tl_avg,tldiffpred_avg)

cascade_ratio_avg = Array{Vector{Float64}}(undef,length(Cvalues));
cascade_ratio_sd = Array{Vector{Float64}}(undef,length(Cvalues));
trophiclevel_avg = Array{Vector{Float64}}(undef,length(Cvalues));

samplesize = 1000;
sample_cascade_ratio_all = Array{Float64}(undef,length(Cvalues),samplesize);
sample_trophiclevel_all = Array{Float64}(undef,length(Cvalues),samplesize);

for c in eachindex(Cvalues)
        
    cascade_ratio_all= reduce(vcat, cascade_ratio_sim[c,:]);
    trophiclevel_all = reduce(vcat, trophiclevel_sim[c,:]);

    # tldiff_sim_all = reduce(vcat, tldiff_sim[c,:]);

    trophiclevel_avg[c], cascade_ratio_avg[c] = moving_average(trophiclevel_all,cascade_ratio_all,0.1);

    _ , cascade_ratio_sd[c] = moving_sd(trophiclevel_all,cascade_ratio_all,0.1);

    #sample points:
    indices = round.(Int, range(1, length(cascade_ratio_all[1,:]), length=samplesize))
    sample_cascade_ratio_all[c,:] = cascade_ratio_all[c,indices];
    sample_trophiclevel_all[c,:] = trophiclevel_all[c,indices];
end

Cvalues = [0.02,0.1];
tlpred_avg = Array{Vector{Float64}}(undef,2);
tldiffpred_avg = Array{Vector{Float64}}(undef,2);
tldiffpred_intercept = Array{Float64}(undef,2)
tldiffpred_slope = Array{Float64}(undef,2)
for c in eachindex(Cvalues)
    tlpred_avg[c], tldiffpred_avg[c] = moving_average(reduce(vcat,trophiclevel_sim[c,:]),reduce(vcat,tldiffpred_sim[c,:]),0.01)

    tldiffpreddata = DataFrame(tlpred_avg=tlpred_avg[c], tldiffpred_avg=tldiffpred_avg[c])
    lmpred_model = lm(@formula(tldiffpred_avg ~ tlpred_avg), tldiffpreddata)
    tldiffpred_intercept[c], tldiffpred_slope[c] = coef(lmpred_model) 

end

Plots.scatter(tlpred_avg[1], tldiffpred_avg[1])
plot!(tlpred_avg[1], tldiffpred_intercept[1] .+ tldiffpred_slope[1] .* tlpred_avg[1], label="Fitted Line", linewidth=2, color=:red)


tlprey_avg, tldiffprey_avg = moving_average(reduce(vcat,trophiclevel_sim[2,:]),reduce(vcat,tldiffprey_sim[2,:]),0.01);

UnicodePlots.scatterplot(tlprey_avg, tldiffprey_avg)

StatsPlots.histogram(tldiffpred_avg)
StatsPlots.histogram(tldiffprey_avg)


#Read in empirical fit
filename = smartpath("data/n_c_propneg_data_g10_kout10.csv")
fitdata = CSV.read(filename, DataFrame)
fit_trophic = Array{Vector{Float64}}(undef,2)
fit_propneg = Array{Vector{Float64}}(undef,2)
for c in eachindex(Cvalues)
    posC = findall(x->x==Cvalues[c],fitdata[!,:C])
    fit_trophic[c], fit_propneg[c] = moving_average(fitdata[!,:trophic][posC],fitdata[!,:propneg][posC],0.01);
end
Plots.plot(fit_trophic[1],fit_propneg[1],xlims=(1,5),ylims=(0,0.6))
Plots.plot!(fit_trophic[2],fit_propneg[2])
pr_state_neg(4,0.02,fit_trophic,fit_propneg)

#Compare to expectation

kout = 4.
kin = 10.
g = 10
#What is the trophic level difference of the predator compared to the prey?


S = 100;
Cvalues = [0.02,0.1]
trophiclevels = collect(1:0.05:5);

reps = 5000;
cascaderatio_conditional = Array{Float64}(undef,length(Cvalues),reps,length(trophiclevels));
cascaderatio_conditionalfit = Array{Float64}(undef,length(Cvalues),reps,length(trophiclevels));
valuevec = collect(-4:2);
probvec = Array{Vector{Float64}}(undef,2)
for c in eachindex(Cvalues)
    probvec_unnormalized = [sum((tldiffpred_avg[1] .>= i) .&& (tldiffpred_avg[1] .< i+1))/length(tldiffpred_avg[1]) for i=valuevec]
    probvec[c] = probvec_unnormalized ./ sum(probvec_unnormalized)
end

@threads for r=1:reps
    
    for c in eachindex(Cvalues)
        
        tldiff = (tdprey=1,tdpred=sample(valuevec, Weights(probvec[c]), 1)[1])

        C = Cvalues[c]
        
        # tldiff_sim_all = reduce(vcat, tldiff_sim[c,:]);

        #Sample from sim tldiffs 
        # diffvec = Int64.(round.(tldiff_sim_all[rand(collect(1:1:length(tldiff_sim_all)),2)]))
        # values = valuevec
        # prob_pred = []
        # probdiff_pred = 1
        # probdiff_prey = [0.,1/3,2/3]

        # diffvec = [1,sample(valuevec, Weights(probvec), 1)[1]];

        
        
        
        for i in eachindex(trophiclevels)

            tl = trophiclevels[i]

            


            # exp_tldiff = tldiffpred_intercept[c] + tldiffpred_slope[c] * tl
            # tldiff = (tdprey=1,tdpred=exp_tldiff)
            # tldiff = (tdprey=1,tdpred=1)

            cascaderatio_conditional[c,r,i] = exp_cascaderatio_conditional(tl,S,C,tldiff,kin,kout,g)
            cascaderatio_conditionalfit[c,r,i] = exp_cascaderatio_conditionalfit(tl,S,C,tldiff,kin,kout,fit_trophic,fit_propneg,g)

        end
    end
end

# cascaderatio_conditional_avg = mean(cascaderatio_conditional,dims=2)[:,1,:]

cascaderatio_conditional_avg = mean(cascaderatio_conditionalfit,dims=2)[:,1,:]


# stopplot = findfirst(x->x==2.5,trophiclevels)
stopplot = length(trophiclevels)
ccplot = Plots.plot(trophiclevels[1:stopplot],cascaderatio_conditional_avg[1,1:stopplot],
    xlims=(1,4),
    ylims=(-1.2,1.2),
    ylabel="Cascade Index",
    xlabel="Trophic Level",
    linewidth=2,
    label="Exp C=0.02",
    framestyle=:box); # Set the y-axis range)
Plots.plot!(ccplot,trophiclevels,cascaderatio_conditional_avg[2,:],
    linewidth=2,
    label="Exp C=0.10");

Plots.scatter!(ccplot,trophiclevel_avg[1],cascade_ratio_avg[1],
    linewidth=2,
    label="Sim C=0.02",
    color="#1f77b4");
Plots.scatter!(ccplot,trophiclevel_avg[2] .+ 0.02,cascade_ratio_avg[2],
    linewidth=2,
    label="Sim C=0.10",
    color="#ff7f0e");

for i=1:length(trophiclevel_avg[1])
    plot!(ccplot,repeat([trophiclevel_avg[1][i]],2),[cascade_ratio_avg[1][i]-cascade_ratio_sd[1][i],cascade_ratio_avg[1][i]+cascade_ratio_sd[1][i]],
        color="#1f77b4",
        linewidth=2,
        alpha=0.5,
        label=false);
end

for i=1:length(trophiclevel_avg[2])
    plot!(ccplot,repeat([trophiclevel_avg[2][i]],2) .+ 0.02,[cascade_ratio_avg[2][i]-cascade_ratio_sd[2][i],cascade_ratio_avg[2][i]+cascade_ratio_sd[2][i]],
        color="#ff7f0e",
        linewidth=2,
        alpha=0.5,
        label=false);
end

ccplot

filename = smartpath("figures/fig_analyticalcascade_v3.pdf")
Plots.savefig(ccplot,filename)


Plots.scatter(
    sample_trophiclevel_all, 
    sample_cascade_ratio_all, 
    alpha=0.1,
    markerstrokewidth=0,
    label="",
    color=:lightblue
)