# Calculate the number of times its state differs from the previous time step.
function compute_proportion_oscillating(s::Array{Int64,2})
    n_species = size(s, 2)
    n_changes = zeros(Int64, n_species)
    for j in 1:n_species
        # Compute state changes after burn-in for each species
        species_states = s[:, j]
        n_changes[j] = sum(abs.(diff(species_states)) .> 0)
    end
    # Calculate the proportion of species that oscillate
    n_oscillating = count(n_changes .> 0)
    prop_oscillating = n_oscillating / n_species
    return prop_oscillating
end

function compute_state_changes(s::Array{Int64,2})
    n_species = size(s, 2)
    n_changes = zeros(Int64, n_species)
    for j in 1:n_species
        # Number of times the species changes state
        n_changes[j] = sum(abs.(diff(s[:, j])) .> 0)
    end
    return n_changes
end

# Track consecutive runs of the same state
function compute_state_durations(s::Array{Int64,2})
    n_species = size(s, 2)
    durations = Vector{Vector{Tuple{Int64, Int64}}}(undef, n_species)
    for j in 1:n_species
        species_states = s[:, j]
        durations[j] = []
        current_state = species_states[1]
        current_duration = 1
        for t in 2:length(species_states)
            if species_states[t] == current_state
                current_duration += 1
            else
                # Append the current state and duration
                push!(durations[j], (current_state, current_duration))
                # Reset for the next state
                current_state = species_states[t]
                current_duration = 1
            end
        end
        # Append the last run
        push!(durations[j], (current_state, current_duration))
    end
    return durations
end

# Aggregate durations for each state (+1 and -1) and compute the average duration in each state for each species.
function compute_average_durations(durations_with_states)
    n_species = length(durations_with_states)
    avg_duration_neg1 = Vector{Float64}(undef, n_species)
    avg_duration_pos1 = Vector{Float64}(undef, n_species)
    for j in 1:n_species
        durations = durations_with_states[j]
        state_durations = Dict{Int64, Vector{Int64}}()
        for (state, duration) in durations
            if !haskey(state_durations, state)
                state_durations[state] = Int[]
            end
            push!(state_durations[state], duration)
        end
        # Compute average durations for each state
        avg_duration_neg1[j] = haskey(state_durations, -1) ? mean(state_durations[-1]) : NaN
        avg_duration_pos1[j] = haskey(state_durations, 1) ? mean(state_durations[1]) : NaN
    end
    # Create DataFrame
    df = DataFrame(
        Species = 1:n_species,
        neg1 = avg_duration_neg1,
        pos1 = avg_duration_pos1
    )
    return df
end


# Record the time steps at which it changes state.
function compute_state_change_times(s::Array{Int64,2})
    n_species = size(s, 2)
    change_times = Vector{Vector{Int64}}(undef, n_species)
    for j in 1:n_species
        species_states = s[:, j]
        change_times[j] = []
        for t in 2:length(species_states)
            if species_states[t] != species_states[t - 1]
                push!(change_times[j], t)
            end
        end
    end
    return change_times
end

#Compute Co-Spin Matrix (Exact Synchrony Between Species)
# function compute_cospin_matrix(s::Array{Int64,2})
#     tmax, n_species = size(s)
#     # Compute state changes for all species
#     state_changes = zeros(Bool, tmax, n_species)
#     for t in 2:tmax
#         state_changes[t, :] .= s[t, :] .!= s[t - 1, :]
#     end
#     # Compute co-spin matrix
#     cospin_matrix = zeros(Int64, n_species, n_species)
#     for i in 1:n_species
#         for j in i + 1:n_species
#             # Number of times both species i and j change state together
#             cospin_matrix[i, j] = sum(state_changes[:, i] .& state_changes[:, j])
#             cospin_matrix[j, i] = cospin_matrix[i, j] # Symmetric matrix
#         end
#     end
#     return cospin_matrix
# end

#Compute Co-Spin Matrix (Delayed(?) Synchrony Between Species)
function compute_cospin_matrix(s::Array{Int64,2}, delay::Int=0)
    tmax, n_species = size(s)
    # Compute state changes for all species
    state_changes = zeros(Bool, tmax, n_species)
    for t in 2:tmax
        state_changes[t, :] .= s[t, :] .!= s[t - 1, :]
    end
    # Compute co-spin matrix with delay
    cospin_matrix = zeros(Int64, n_species, n_species)
    for i in 1:n_species
        for j in i + 1:n_species
            # Calculate indices for valid time steps considering the delay
            if delay >= 0
                # For non-negative delays
                valid_t = 2:(tmax - delay)
                counts = sum(state_changes[valid_t, i] .& state_changes[valid_t .+ delay, j])
            else
                # For negative delays
                valid_t = (2 - delay):tmax
                counts = sum(state_changes[valid_t .+ delay, i] .& state_changes[valid_t, j])
            end
            cospin_matrix[i, j] = counts
            cospin_matrix[j, i] = counts # Symmetric matrix
        end
    end
    #NOTE: Probably want to ./ tmax
    return cospin_matrix ./ tmax
end


# Compute Total Time in Each State
function compute_total_time_in_states(s::Array{Int64,2})
    n_species = size(s, 2)
    time_neg1 = zeros(Int64, n_species)
    time_pos1 = zeros(Int64, n_species)
    for j in 1:n_species
        species_states = s[:, j]
        counts = StatsBase.countmap(species_states)
        # Get counts for -1 and 1, defaulting to 0 if the key is not present
        time_neg1[j] = get(counts, -1, 0)
        time_pos1[j] = get(counts, 1, 0)
    end
    # Create DataFrame
    df = DataFrame(
        Species = 1:n_species,
        TotalTime_neg1 = time_neg1,
        TotalTime_pos1 = time_pos1
    )
    return df
end

# Compute State Transition Probabilities
using DataFrames

function compute_state_transition_probabilities(s::Array{Int64,2})
    n_species = size(s, 2)
    n_time_steps = size(s, 1)
    # Initialize arrays to store probabilities
    pos1_neg1_probs = zeros(Float64, n_species)
    pos1_pos1_probs = zeros(Float64, n_species)
    neg1_neg1_probs = zeros(Float64, n_species)
    neg1_pos1_probs = zeros(Float64, n_species)
    for j in 1:n_species
        species_states = s[:, j]
        transitions = Dict{Tuple{Int64, Int64}, Int64}()
        total_transitions = 0
        for t in 2:n_time_steps
            prev_state = species_states[t - 1]
            curr_state = species_states[t]
            key = (prev_state, curr_state)
            transitions[key] = get(transitions, key, 0) + 1
            total_transitions += 1
        end
        # Retrieve counts for each possible transition
        pos1_neg1 = get(transitions, (1, -1), 0)
        pos1_pos1 = get(transitions, (1, 1), 0)
        neg1_neg1 = get(transitions, (-1, -1), 0)
        neg1_pos1 = get(transitions, (-1, 1), 0)
        # Compute probabilities
        if total_transitions > 0
            pos1_neg1_probs[j] = pos1_neg1 / total_transitions
            pos1_pos1_probs[j] = pos1_pos1 / total_transitions
            neg1_neg1_probs[j] = neg1_neg1 / total_transitions
            neg1_pos1_probs[j] = neg1_pos1 / total_transitions
        else
            # Handle cases with no transitions
            pos1_neg1_probs[j] = NaN
            pos1_pos1_probs[j] = NaN
            neg1_neg1_probs[j] = NaN
            neg1_pos1_probs[j] = NaN
        end
    end
    # Create DataFrame
    df = DataFrame(
        Species = 1:n_species,
        pos1_pos1 = pos1_pos1_probs,
        pos1_neg1 = pos1_neg1_probs,
        neg1_pos1 = neg1_pos1_probs,
        neg1_neg1 = neg1_neg1_probs
    )
    return df
end

function cascade_effect(s::Array{Int64,2},Asort,offset)
    tmax = size(s,1)
    n_species = size(s, 2)
    # cascade_effect = Array{Float64}(undef,size(s))

    coeff = 2;
    # maxratio = log((coeff - 1)/(coeff + 1))
    # minratio = log((coeff + 1)/(coeff - 1))

    cascade_ratio = Array{Float64}(undef,n_species)
    cascade_ratio_interactions = Vector{Float64}(undef,0)

    for i = 1:n_species

        #how many predators does species i have?
        preds = findall(x->x==1,Asort[i,:])
        n_preds = length(preds)
        cascade_ratio_perpred = Array{Float64}(undef,n_preds)

        if n_preds > 0
            for j = 1:n_preds
                # +P positions
                posP = findall(x->x==1,s[1:(tmax-offset),j])
                # -P positions
                negP = setdiff(collect(1:1:(tmax-offset)),posP)

                #NOTE: consider including +2 offset to posP and negP
                avg_state_posP = mean(s[posP .+ offset,i])
                avg_state_negP = mean(s[negP .+ offset,i])

                cascade_ratio_perpred[j] = log((coeff + avg_state_posP) / (coeff + avg_state_negP))
            end

            #If predators never flip, there will be ratios = NaN.
            #Filter these out
            filtered_cascade_ratio_perpred = filter(!isnan,cascade_ratio_perpred)
            
            if length(filtered_cascade_ratio_perpred) > 0
                append!(cascade_ratio_interactions,vec(filtered_cascade_ratio_perpred))
            end

            cascade_ratio_mean = mean(filtered_cascade_ratio_perpred)
            cascade_ratio[i] = cascade_ratio_mean
            #Scaled ratio between 0 and 1
            # cascade_ratio[i] = (maxratio - cascade_ratio_mean) / (maxratio - minratio)

        else
            cascade_ratio[i] = NaN #Or nothing?
        end
    end

    return cascade_ratio, cascade_ratio_interactions

end