function plotweb(A, s, tl)
    # Ensure reproducibility of random x-coordinates
    Random.seed!(1)
    
    # Transpose A to match the direction of edges
    g = Graphs.DiGraph(A')

    n = Graphs.nv(g)
    x_coords = rand(n)  # Random x-coordinates between 0 and 1
    y_coords = tl       # y-coordinates based on trophic levels

    coords = hcat(x_coords, y_coords)

    # Map s from [-1, 1] to [0, 1] for greyscale mapping
    s_mapped = (s .+ 1) ./ 2

    # Create grey colors for vertices based on s
    vertex_colors = [Gray(s_i) for s_i in s_mapped]

    # Edge color with transparency (matching '#6495ED30' in R)
    edge_col = RGBA(100/255, 149/255, 237/255, 48/255)

    # Plot the graph using GraphPlot.jl
    # Plot the graph using GraphPlot.jl
    
    locs_x = x_coords;
    locs_y = y_coords;
    GraphPlot.gplot(
        g,locs_x, locs_y
    )
end

# nodefillc = vertex_colors,
# nodelabel = false,
# nodestrokec = :black,        # Replace nodeoutlinecolor
# nodesize = 0.1,              # Adjust nodesize for better scaling
# edgestrokec = edge_col,      # Replace edgecolor
# arrowlengthfrac = 0.05,      # Adjust arrow size
# linetype = :solid,           # Edge style
# outangle = 0                 # Arrow direction