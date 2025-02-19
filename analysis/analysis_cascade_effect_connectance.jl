
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

@load filename reps gvec lgvec koutvec lkoutvec Cvec lCvec S kin tmax sigma mean_cascade_ratio



Plots.heatmap(gvec,koutvec,mean_cascade_ratio[3,:,:]')


using Measures
using LaTeXStrings
using PGFPlotsX
pgfplotsx()
p1 = Plots.scatter(vec(mean_cascade_ratio[1,:,:]),vec(mean_cascade_ratio[2,:,:]),
    xlims=(-0.9,0.0),
    ylims=(-0.9,0.0),
    xlabel=L"\mathrm{Cascade~effect,~}\bar{\psi}~(C=0.02)",
    ylabel=L"\mathrm{Cascade~effect,~}\bar{\psi}~(C=C_{\rm alt})",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    alpha=0.3,
    markerstrokewidth=0,
    label=L"C_{\rm alt}=0.03",
    frame=:box,
    legendfontsize=14,
    framestyle = :box,
    foreground_color_legend = nothing,
    legend = (0.6, 0.2));
Plots.scatter!(p1,vec(mean_cascade_ratio[1,:,:]),vec(mean_cascade_ratio[3,:,:]),
    alpha=0.3,
    markerstrokewidth=0,
    label=L"C_{\rm alt}=0.04",
    legendfontsize=14);
Plots.plot!(p1,collect(-1:0.1:1),collect(-1:0.1:1),label=false,color=:black);
# Plots.scatter!(p1,
#     [mean_cascade_ratio[1,11,12]], [mean_cascade_ratio[2,11,12]],
#     color=:gray, label=false);
# Plots.scatter!(p1,
#     [mean_cascade_ratio[1,11,12]], [mean_cascade_ratio[3,11,12]],
#     color=:black, label=false);

filename = smartpath("figures/fig_cascade_effect_conn.pdf")
Plots.savefig(p1,filename)

# Flip ratio bc both negative numbers!
# Now a value > 1 will mean an increase in cascade effect with increased connectance
ratio = 1 ./ (mean_cascade_ratio[3,:,:] ./ mean_cascade_ratio[1,:,:]);

#What percentage is increasing?
sum(vec(filter(!isnan,ratio)) .>= 1)/length(vec(filter(!isnan,ratio)))

ratio[11,12]

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

filename = smartpath("figures/fig_cascade_effect_ratio.pdf");
Plots.savefig(ratioplot,filename)

using Measures
combinedplot = plot(p1, ratioplot, 
    layout = (1, 2), 
    size = (1000, 400),
    left_margin=10mm,  # Add more space to the left
    bottom_margin=5mm);
filename = smartpath("figures/fig_cascade_effect_connratio.pdf");
Plots.savefig(combinedplot,filename)


size(mean_cascade_ratio)

koutvec[21]
gvec[11]
#(11,12) is the coordinate for g = 10, kout = 10
mean_cascade_ratio[1,11,12]
mean_cascade_ratio[2,11,12]
mean_cascade_ratio[3,11,12]


