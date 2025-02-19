function plotweb(A::Matrix{Bool}, tl::Vector{Float64}, s::Vector{Int64} = ones(Int64, size(A, 1)))
    R"""
    set.seed(1);
    library(igraph)
    library(RColorBrewer)
    pal <- brewer.pal(3,"Set1")
    # fw_g <- graph.adjacency($(A));
    fw_g <- graph_from_adjacency_matrix($(A));
    trophic <- as.numeric($tl);

    unique_trophic <- unique(trophic)
    # Get indices for nodes with this trophic level
    x_coords <- numeric(length(trophic))
    for (i in unique_trophic) {
        # Get the indices of nodes at this trophic level
        idx <- which(trophic == i)
        if (length(idx) == 1) {
            # If only one species, place it at 0.5
            x_coords[idx] <- jitter(0.5,3)
        } else {
            # Generate evenly spaced x-coordinates for this trophic level
            x_coords[idx] <- jitter(seq(0, 1, length.out = length(idx)),2)
        }
    }

    coords <- cbind(x_coords,trophic);
    plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color="black",vertex.label=NA,vertex.frame.color="black", vertex.color=grey(1-($s+1)/2))
    """
end

