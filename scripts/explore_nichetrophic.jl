
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

using CSV
using Plots
using UnicodePlots
using StatsPlots
using LsqFit


# A food web model exzmple
S=100; C=0.02;
A,niche = nichemodelweb(S,C);
g = Graphs.DiGraph(A);
tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6);
plotweb(A,tl)

# Compute row sums and column sums
A = A_ex
tl = round.(TrophInd(Array{Float64}(A))[!,:TL]; sigdigits=6);
row_sums = sum(A_ex, dims=2)[:, 1] # Sum along rows
col_sums = sum(A_ex, dims=1)[1, :] # Sum along columns
# Find indices where both row sums and column sums are zero
indices = findall(row_sums .== 0 .&& col_sums .== 0)




S=100; C=0.02;
reps = 5000
niche_r = Array{Vector{Float64}}(undef,reps);
tl_r = Array{Vector{Float64}}(undef,reps);
numpreds_r = Array{Vector{Float64}}(undef,reps);
numprey_r = Array{Vector{Float64}}(undef,reps);

# Container to store problematic cases
problem_cases = []

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

    if length(tl) !== length(niche)
        println("Uh Oh at iteration $r")
        
        # Save the problematic case
        push!(problem_cases, (A = A, niche = niche, tl = tl, iteration = r))
    end

    # How many predators for species i?
    numpreds = vec(sum(A,dims=2))
    # How many prey for species i?
    numprey = vec(sum(A,dims=1))
    

    numpreds_r[r] = numpreds;
    numprey_r[r] = numprey;
    niche_r[r] = niche;
    tl_r[r] = tl;

end

#collapse
niche_all = reduce(vcat, niche_r);
tl_all =  reduce(vcat, tl_r);
prey_all = reduce(vcat,numprey_r);
preds_all = reduce(vcat,numpreds_r);


# Sort by niche value
sorted_idx = sortperm(niche_all)
niche_sorted = niche_all[sorted_idx]
prey_sorted = prey_all[sorted_idx]
preds_sorted = preds_all[sorted_idx]
trophic_sorted = tl_all[sorted_idx]

mean_niche , mean_prey = moving_average(niche_sorted,prey_sorted,0.05)
_, mean_preds = moving_average(niche_sorted,preds_sorted,0.05)
_, mean_trophic = moving_average(niche_sorted,trophic_sorted,0.05)

# # Define a window width (e.g., 0.05 means ±0.025 around each point)
# window_width = 0.05

# # We'll iterate over a set of points along the niche axis to compute a sliding mean.
# # For simplicity, let's just compute it at each data point:
# mean_prey = similar(niche_sorted, Float64);
# mean_preds = similar(niche_sorted, Float64);
# mean_trophic = similar(niche_sorted, Float64);


# # For efficient searching, we can use binary search.
# # using Base: searchsortedfirst, searchsortedlast

# @threads for i in eachindex(niche_sorted)
#     center = niche_sorted[i]
#     low_bound = center - window_width/2
#     high_bound = center + window_width/2

#     # Find the indices of points in the window using binary search
#     low_idx = searchsortedfirst(niche_sorted, low_bound)
#     high_idx = searchsortedlast(niche_sorted, high_bound)

#     # Compute mean prey count within this window
#     mean_prey[i] = mean(prey_sorted[low_idx:high_idx])
#     mean_preds[i] = mean(preds_sorted[low_idx:high_idx])

#     mean_trophic[i] = mean(trophic_sorted[low_idx:high_idx])
# end


# Now mean_prey[i] holds the mean prey count in the neighborhood of niche_sorted[i]
# If you want a smoother curve, you might plot niche_sorted vs mean_prey:
# using Plots

using LaTeXStrings
using PGFPlotsX
pgfplotsx()
nichepreyplot = Plots.plot(mean_niche, mean_prey, 
    linewidth=2,
    xlabel=L"\mathrm{Niche~value,~} \eta", 
    ylabel="Mean prey count",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    legendfontsize=14,
    label="Simulated",
    framestyle = :box,
    foreground_color_legend = nothing,
    legend = (0.1, 0.9)
    );
analytical_niche = collect(0:0.1:1)
analytical_prey = (S - 1) .* 2 .* C .* analytical_niche;
Plots.plot!(nichepreyplot,analytical_niche,analytical_prey,
    linewidth=2,
    legendfontsize=14,
    label="Expected");

nichepredplot = Plots.plot(mean_niche, mean_preds, 
    linewidth=2,
    xlabel=L"\mathrm{Niche~value,~} \eta", 
    ylabel="Mean predator count",
    xlabelfontsize = 14,  # Set x-axis label font size
    ylabelfontsize = 14,   # Set y-axis label font size
    tickfontsize = 12,    # Font size for axis tick numbers
    framestyle = :box,
    label=false
    # legendfontsize=14,
    # label="Simulated",
    # framestyle = :box,
    # foreground_color_legend = nothing
    );
analytical_niche = collect(0:0.1:1)
analytical_pred = (S - 1) .* 2 .* C .* (1 .- analytical_niche);
Plots.plot!(nichepredplot,analytical_niche,analytical_pred,
    linewidth=2,
    label=false);


predpreyplot = plot(nichepreyplot, nichepredplot,  layout = (1, 2), size=(900, 400));

filename = smartpath("figures/fig_predpreyniche.pdf")
Plots.savefig(predpreyplot,filename)






avg_niche_all, avg_prey_all = moving_average(niche_all,prey_all,0.02);
nichepreyplot = Plots.plot(avg_niche_all, avg_prey_all, xlabel="Niche value", ylabel="Mean prey count");
Plots.plot!(nichepreyplot,analytical_niche,analytical_prey);
Plots.plot!(nichepreyplot,niche_sorted,mean_prey);
nichepreyplot

nichepredsplot = Plots.plot(niche_sorted, mean_preds, xlabel="Niche value", ylabel="Mean pred count")
analytical_preds = (S - 1) .* 2 .* C .* (1 .- analytical_niche);
Plots.plot!(analytical_niche,analytical_preds)

nichetrophicplot = Plots.plot(niche_sorted, mean_trophic, xlabel="Niche value", ylabel="Mean trophic level")





# Data points
x = niche_sorted
y = mean_trophic

# Define a power-law model: y = a * x^b + c
power_law_model(x, p) = 1 .+ p[1] * x.^p[2] #.+ p[3]

# Initial guesses for parameters [a, b, c]
p0 = [1.0, 0.5] #, 1.0]

# Fit the model to the data
power_law_fit = curve_fit(power_law_model, x, y, p0)

# Extract fitted parameters
a, b = power_law_fit.param

println("Fitted power-law equation: y = 1 + $(a) * x^$(b)")

# Generate fitted curve
x_fit = range(minimum(x), stop=maximum(x), length=100)  # Smooth x range
y_fit = power_law_model(x_fit, power_law_fit.param)  # Predicted y values

# Plot the data and fitted curve
plot(x, y, label="Data Points", xlabel="x", ylabel="y", title="Power-Law Fit")
plot!(x_fit, y_fit, label="Fitted Power-Law", lw=2)



# Fit trophic level across C and n

using Statistics
using LsqFit # or LsqFit or another fitting package

# Suppose you have functions:
# run_simulation(C) -> returns arrays n_values, TL_values for given C
# For demonstration, let's assume you have these arrays from somewhere.
S=100;
C_values = collect(0.01:0.01:0.1) # or a range like 0.01:0.01:0.2
n_all = Float64[]
C_all = Float64[]
TL_all = Float64[]

for C in C_values

    reps = 5000
    niche_r = Array{Vector{Float64}}(undef,reps);
    tl_r = Array{Vector{Float64}}(undef,reps);
    numpreds_r = Array{Vector{Float64}}(undef,reps);
    numprey_r = Array{Vector{Float64}}(undef,reps);

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

        # How many predators for species i?
        numpreds = vec(sum(A,dims=2))
        # How many prey for species i?
        numprey = vec(sum(A,dims=1))

        numpreds_r[r] = numpreds;
        numprey_r[r] = numprey;
        niche_r[r] = niche;
        tl_r[r] = tl;

    end

    #collapse
    niche_all = reduce(vcat, niche_r);
    tl_all =  reduce(vcat, tl_r);
    prey_all = reduce(vcat,numprey_r);
    preds_all = reduce(vcat,numpreds_r);

    # Sort by niche value
    sorted_idx = sortperm(niche_all)
    niche_sorted = niche_all[sorted_idx]
    prey_sorted = prey_all[sorted_idx]
    preds_sorted = preds_all[sorted_idx]
    trophic_sorted = tl_all[sorted_idx]

    # Define a window width (e.g., 0.05 means ±0.025 around each point)
    window_width = 0.2

    # We'll iterate over a set of points along the niche axis to compute a sliding mean.
    # For simplicity, let's just compute it at each data point:
    mean_prey = similar(niche_sorted, Float64)
    mean_preds = similar(niche_sorted, Float64)
    mean_trophic = similar(niche_sorted, Float64)


    # For efficient searching, we can use binary search.
    # using Base: searchsortedfirst, searchsortedlast

    @threads for i in eachindex(niche_sorted)
        center = niche_sorted[i]
        low_bound = center - window_width/2
        high_bound = center + window_width/2

        # Find the indices of points in the window using binary search
        low_idx = searchsortedfirst(niche_sorted, low_bound)
        high_idx = searchsortedlast(niche_sorted, high_bound)

        # Compute mean prey count within this window
        mean_prey[i] = mean(prey_sorted[low_idx:high_idx])
        mean_preds[i] = mean(preds_sorted[low_idx:high_idx])

        mean_trophic[i] = mean(trophic_sorted[low_idx:high_idx])
    end




    # Run simulation at this C
    # n_values, TL_values = run_simulation(C)
    # Append results
    append!(n_all, niche_sorted)
    append!(C_all, fill(C, length(niche_sorted)))
    append!(TL_all, mean_trophic)
end

df = DataFrame(n = n_all, C = C_all, TL = TL_all)
filename = smartpath("data/n_c_tl_data.csv")
CSV.write(filename, df)




#####################################
# Probability of STATE given niche value/trophic level


S=100; C=0.2;
reps = 10000;

# Prepare to collect outputs from all repetitions
Asort_r = Vector{Matrix{Int64}}(undef, reps);
tlsort_r = Vector{Vector{Float64}}(undef, reps);
nichesort_r = Vector{Vector{Float64}}(undef, reps);
s_r = Vector{Matrix{Int64}}(undef, reps);

@threads for r = 1:reps


    gprim = rand(collect(0.0:0.1:25))
    kout = rand(collect(0.0:0.1:15.0))
    kin = rand(collect(0.0:0.1:15.0))
    tmax = 300
    sigma = 0.01


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

#What is the state for each species?
for r = 1:reps
    Asort = Asort_r[r]
    tlsort = tlsort_r[r]
    nichesort = nichesort_r[r]
    s = s_r[r]

    propneg_sp = Vector{Float64}(undef,size(Asort)[1]);
    for i = 1:size(Asort)[1]
        propneg_sp[i] = sum(s[end-100:end,i] .== -1) / 101
    end

    propneg_r[r] = propneg_sp
    niche_r[r] = nichesort
    trophic_r[r] = tlsort

end

propneg_all = reduce(vcat, propneg_r);
niche_all = reduce(vcat, niche_r);
trophic_all = reduce(vcat, trophic_r);

# Sort by niche value
sorted_idx = sortperm(niche_all)
niche_sorted = niche_all[sorted_idx]
propneg_sorted = propneg_all[sorted_idx]
trophic_sorted = trophic_all[sorted_idx]

# Define a window width (e.g., 0.05 means ±0.025 around each point)
window_width = 0.05

# Prepare arrays for the sliding mean results
mean_propneg = similar(niche_sorted, Float64)
mean_trophic = similar(niche_sorted, Float64)

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
    mean_trophic[i] = mean(trophic_sorted[low_idx:high_idx])
end

propnegplot = Plots.plot(niche_sorted, mean_propneg, xlabel="Niche value", ylabel="Mean proportion negative")
Plots.plot!(propnegplot,niche_sorted, 1 .- mean_propneg, xlabel="Niche value", ylabel="Mean proportion negative")

analytical_trophic = 1 .+ 11.64 .* (niche_sorted .* C) .^(0.49);
propnegplot_trophic = Plots.plot(analytical_trophic, mean_propneg, xlabel="Expected Trophic level", ylabel="Mean proportion negative")

trophicplot = Plots.plot(niche_sorted, mean_trophic, xlabel="Niche value", ylabel="Mean trophic")

propnegplot_trophic = Plots.plot(trophic_sorted, mean_propneg, xlabel="Trophic level", ylabel="Mean proportion negative")


#Run this over many values of C

S=100;
kin = 5.;
C_values = collect(0.01:0.01:0.1) # or a range like 0.01:0.01:0.2
gprim_lowerbound = [0.,12.];
gprim_upperbound = [12.,25.];
kout_lowerbound = [0.,7.]
kout_upperbound = [7.,15.];

index_name = ["low","high"]

#Across quadrants
for gindex in eachindex(gprim_upperbound)
    for koutindex in eachindex(kout_upperbound)


        n_all = Float64[]
        C_all = Float64[]
        propneg_all = Float64[]

        for C in C_values

            reps = 5000
            S=100; 

            # Prepare to collect outputs from all repetitions
            Asort_r = Vector{Matrix{Int64}}(undef, reps);
            tlsort_r = Vector{Vector{Float64}}(undef, reps);
            nichesort_r = Vector{Vector{Float64}}(undef, reps);
            s_r = Vector{Matrix{Int64}}(undef, reps);

            @threads for r = 1:reps

                gprim = rand(collect(gprim_lowerbound[gindex]:0.1:gprim_upperbound[gindex]))
                kout = rand(collect(kout_lowerbound[koutindex]:0.1:kout_lowerbound[koutindex]))
                tmax = 300
                sigma = 0.01


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

            #What is the state for each species?
            for r = 1:reps
                Asort = Asort_r[r]
                tlsort = tlsort_r[r]
                nichesort = nichesort_r[r]
                s = s_r[r]

                propneg_sp = Vector{Float64}(undef,size(Asort)[1]);
                for i = 1:size(Asort)[1]
                    propneg_sp[i] = sum(s[end-100:end,i] .== -1) / 101
                end

                propneg_r[r] = propneg_sp
                niche_r[r] = nichesort
                trophic_r[r] = tlsort

            end

            propneg_flattened = reduce(vcat, propneg_r);
            niche_flattened = reduce(vcat, niche_r);
            trophic_flattened = reduce(vcat, trophic_r);

            # Sort by niche value
            sorted_idx = sortperm(niche_flattened)
            niche_sorted = niche_flattened[sorted_idx]
            propneg_sorted = propneg_flattened[sorted_idx]
            trophic_sorted = trophic_flattened[sorted_idx]

            # Define a window width (e.g., 0.05 means ±0.025 around each point)
            window_width = 0.05

            # Prepare arrays for the sliding mean results
            mean_propneg = similar(niche_sorted, Float64)
            mean_trophic = similar(niche_sorted, Float64)

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
                mean_trophic[i] = mean(trophic_sorted[low_idx:high_idx])
            end

            append!(n_all, niche_sorted)
            append!(C_all, fill(C, length(niche_sorted)))
            append!(propneg_all, mean_propneg)

        end


        #Export to fit with mathematica
        df = DataFrame(n = n_all, C = C_all, propneg = propneg_all)

        file = string("data/n_c_propneg_data_",index_name[gindex],"g_",index_name[koutindex],"kout.csv")
        filename = smartpath(file)
        CSV.write(filename, df)
    end
end



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


#Export to fit with mathematica
df = DataFrame(n = n_all, C = C_all, propneg = propneg_all, durpos1 = dur_pos1_all, durneg1 = dur_neg1_all,trophic = trophic_all,propneg_sd = propneg_std_all);

gfilevalue = Int64(round(gprim));
koutfilevalue = Int64(round(kout));

file = string("data/n_c_propneg_data_","g",gfilevalue,"_kout",koutfilevalue,".csv")
filename = smartpath(file)
CSV.write(filename, df)

cvalue = findall(x->x==0.1,C_all);
UnicodePlots.scatterplot(n_all[cvalue],propneg_all[cvalue])

UnicodePlots.scatterplot(n_all[cvalue],dur_neg1_all[cvalue])
UnicodePlots.scatterplot(n_all[cvalue],dur_pos1_all[cvalue])

cvalue = findall(x->x==0.02,C_all);
est_trophic_clow = 1 .+ 11.6 .* (0.02 .* n_all[cvalue]) .^ 0.49;
cvalue = findall(x->x==0.1,C_all);
est_trophic_chigh = 1 .+ 11.6 .* (0.1 .* n_all[cvalue]) .^ 0.49;

cvalue = findall(x->x==0.1,C_all);
UnicodePlots.scatterplot(n_all[cvalue],trophic_all[cvalue])


# Create the first scatter plot
cvalue = findall(x->x==0.02,C_all);
durplotSminus = Plots.plot(
    # trophic_clow,
    trophic_all[cvalue],
    dur_neg1_all[cvalue],
    xlabel = "Trophic level",
    ylabel = "Duration in S-"
);
cvalue = findall(x->x==0.1,C_all);
Plots.plot!(durplotSminus,
    # n_all[cvalue],
    # trophic_chigh,
    trophic_all[cvalue],
    dur_neg1_all[cvalue]
);

# Create the second scatter plot
cvalue = findall(x->x==0.02,C_all);
durplotSplus = Plots.plot(
    # n_all[cvalue],
    # trophic_clow,
    trophic_all[cvalue],
    dur_pos1_all[cvalue],
    # title = "Positive State Durations",
    xlabel = "Trophic level",
    ylabel = "Duration in S+"
);
cvalue = findall(x->x==0.1,C_all);
Plots.plot!(durplotSplus,
    # n_all[cvalue],
    # trophic_chigh,
    trophic_all[cvalue],
    dur_pos1_all[cvalue]
);

# Combine both plots side by side
durplot = plot(durplotSplus, durplotSminus,  layout = (1, 2), size=(800, 400));

filename = smartpath("figures/fig_duration_trophic2.pdf")
Plots.savefig(durplot,filename)


#Just one panel with C = 0.02

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
    ylabel = L"\mathrm{Duration}",
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

filename = smartpath("figures/fig_duration_trophic2.pdf")
Plots.savefig(durplotS,filename)
