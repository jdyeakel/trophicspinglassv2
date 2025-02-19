
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

Plots.scatter(vcat(cascade_tl[3,:]...),vcat(cascade_ratio[3,:]...))

tlcascademean = moving_average(vcat(cascade_tl[3,:]...),vcat(cascade_ratio[3,:]...),0.01);
Plots.scatter(tlcascademean)

tlcascadesim = DataFrame(tl=tlcascademean[1],cascade=tlcascademean[2])
filename = smartpath("data/tlcascadeC04.csv")
CSV.write(filename, tlcascadesim)



[mean(vcat(cascade_ratio[i,:]...)) for i=1:length(Cvec)]

dataC1 = vcat(tlreps[1, :]...);
dataC2 = vcat(tlreps[2, :]...);
dataC3 = vcat(tlreps[3, :]...);

StatsPlots.histogram(dataC1)
StatsPlots.histogram(dataC2)
StatsPlots.histogram(dataC3)

using StatsBase
histC1 = fit(Histogram,dataC1)
histC2 = fit(Histogram,dataC2)
histC3 = fit(Histogram,dataC3)

histC1.edges

#Interpolate for TL 1.5
histC1.weights[2] = Int64(round(mean([histC1.weights[1],histC1.weights[3]])))
histC2.weights[2] = Int64(round(mean([histC2.weights[1],histC2.weights[3]])))
histC3.weights[2] = Int64(round(mean([histC3.weights[1],histC3.weights[3]])))

wC1 = histC1.weights ./ sum(histC1.weights)
wC2 = histC2.weights ./ sum(histC2.weights)
wC3 = histC3.weights ./ sum(histC3.weights)

weightplot = Plots.scatter(collect(histC1.edges[1])[1:end-1],wC1)
Plots.plot!(weightplot,collect(histC1.edges[1])[1:end-1],wC1)
Plots.scatter!(weightplot,collect(histC2.edges[1])[1:end-1],wC2)
Plots.plot!(weightplot,collect(histC2.edges[1])[1:end-1],wC2)
Plots.scatter!(weightplot,collect(histC3.edges[1])[1:end-1],wC3)
Plots.plot!(weightplot,collect(histC3.edges[1])[1:end-1],wC3)

TL = collect(histC3.edges[1])[1:end-1]
df = DataFrame(TL=TL, wC1=[wC1;0], wC2=wC2, wC3=wC3)

filename = smartpath("data/tl_distributions.csv")
CSV.write(filename, df)




Plots.heatmap(gvec,koutvec,mean_cascade_ratio[3,:,:]')



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

filename = smartpath("figures/fig_cascade_effect_conn.pdf")
Plots.savefig(p1,filename)

# Flip ratio bc both negative numbers!
# Now a value > 1 will mean an increase in cascade effect with increased connectance
ratio = 1 ./ (mean_cascade_ratio[3,:,:] ./ mean_cascade_ratio[1,:,:]);

ratio[11,12]

ratioplot = Plots.heatmap(gvec,koutvec,ratio',
clims=(0,1.5),
xlabel="Primary productivity g",
ylabel="Top-down effect kout",
colorbar_title = "Cascade effect ratio Chigh/Clow");

filename = smartpath("figures/fig_cascade_effect_ratio.pdf");
Plots.savefig(ratioplot,filename)




size(mean_cascade_ratio)

koutvec[21]
gvec[11]
#(11,12) is the coordinate for g = 10, kout = 10
mean_cascade_ratio[1,11,12]
mean_cascade_ratio[2,11,12]
mean_cascade_ratio[3,11,12]