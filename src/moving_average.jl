function moving_average(xcoords, ycoords, window)
    
    # Ensure inputs are valid
    if length(xcoords) != length(ycoords)
        throw(ArgumentError("xcoords and ycoords must have the same length."))
    end
    if window <= 0
        throw(ArgumentError("Window size must be greater than 0."))
    end

    # Filter out pairs where either x or y is NaN
    valid_indices = .!isnan.(xcoords) .& .!isnan.(ycoords)
    xcoords = xcoords[valid_indices]
    ycoords = ycoords[valid_indices]

    # Sort the data based on xcoords to ensure sliding window works correctly
    sorted_indices = sortperm(xcoords)
    x_sorted = xcoords[sorted_indices]
    y_sorted = ycoords[sorted_indices]

    avg_xcoords = Float64[]
    avg_ycoords = Float64[]

    # Apply the moving average using the specified window size
    i = 1
    while i <= length(x_sorted)
        # Define the window range
        current_window_x = x_sorted[i]
        window_indices = findall(x -> abs(x - current_window_x) <= window / 2, x_sorted)
        
        # Calculate the mean of the x and y values in the window
        avg_x = mean(x_sorted[window_indices])
        avg_y = mean(y_sorted[window_indices])

        # Append the averages to the result arrays
        push!(avg_xcoords, avg_x)
        push!(avg_ycoords, avg_y)

        # Move the starting point of the sliding window
        i += length(window_indices)
    end

    return avg_xcoords, avg_ycoords
end