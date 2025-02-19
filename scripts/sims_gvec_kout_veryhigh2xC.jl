
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

# using Plots
# using UnicodePlots
# using StatsPlots

using ProgressMeter
# using Logging

# global_logger(SimpleLogger(stderr))

subfile = "_veryhigh2xC"

reps = 200
gvec = collect(0:1:25)
lgvec = length(gvec)
koutvec = collect(0:0.5:15)
lkoutvec = length(koutvec)

S = 100; C = 0.02*(1 + 1.)
kin = 5.0
tmax = 300
sigma = 0.01


# Load variables from saved file
filename_preamble = smartpath(string("data/gvec_kout",subfile,"/preamble.jld2"))
@save filename_preamble reps gvec koutvec S C kin tmax sigma

# Before the loops
# @info "Starting simulation"

# Loop over parameter combinations
# @showprogress 1 "Computing..." 
for i in 1:lgvec
    # @info "Starting gvec index $i"
    @threads for j in 1:lkoutvec
        # @info "Starting koutvec index $j"
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
            # @info "Thread $(Threads.threadid()) completed repetition $r"
        end

        # After all reps, save data to a single file
        filename = smartpath(string("data/gvec_kout",subfile,"/g$(i)_kout$(j).jld2"))
        @save filename Asort_list tlsort_list s_list gprim kout kin sigma tmax

        # @info "Completed koutvec index $j"


    end
    # @info "Completed gvec index $i"

    # GC.gc()  # Force garbage collection after each rep
    println(i/lgvec)

end

# @info "Simulation completed"




# namespace = string("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/propcycle.jld");
# @save namespace propcycle;

# # @load namespace propcycle;

# mpropcycle = mean(propcycle,dims=3)[:,:,1];

# namespace = string("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/propcycle.pdf");
# R"""
# library(RColorBrewer)
# library(fields)
# pal=colorRampPalette(brewer.pal(9,"Spectral"))(100);
# pdf($namespace,width=8,height=7)
# image.plot(x=$gvec,y=$koutvec/$kin,z=$mpropcycle,col=pal,xlab='Primary productivity',ylab='Top-down influence (k_out/k_in)',legend.args=list(text='  Prop. cycle',side=3,line=1))
# dev.off()
# """



# namespace = string("$(homedir())/Dropbox/PostDoc/2018_trophicspinglass/spin.pdf");
# R"""
# library(RColorBrewer)
# library(igraph)
# pal = brewer.pal(9,"Blues")
# #pal = rev(brewer.pal($tlout-1,'Spectral'))
# pdf($namespace,width=8,height=10)
# par(mfrow=c(2,1))
# image(x=$(collect(1:tmax)),y=$(collect(1:S)),$(s),col=c(pal[2],pal[7]),xlab='Time',ylab='Species ID')
# #plot($(Delta[1,:]),type='l',ylim=c(min($Delta),max($Delta)),col=pal[1],lwd=2)
# g = graph_from_adjacency_matrix($Asort)
# plot(g,vertex.size=2,edge.arrow.size=0.2,vertex.label=NA,vertex.color='lightblue')
# dev.off()
# """

#for i=tlout-1:-1:2
#    R"lines($(Delta[i,:]),col=pal[$i],lwd=2)"
#end
# R"dev.off()"
#Pattern of traveling peaks shows the cascade direction/strength






#t=50;
#plotweb(Asort,s[t,:],tlsort)






#webs
#tend = 100;
#trange = collect(tend-4+1:tend);
#R"par(mfrow=c(2,2))"
#for t=tend-4+1:tend
#    plotweb(Asort,s[t,:],tlsort)
#end
