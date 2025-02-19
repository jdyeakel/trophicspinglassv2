function find_motifs(g::DiGraph, motif::DiGraph)
    # Enumerate all subgraphs of size n (same as motif)
    n = Graphs.nv(motif)
    nodes = collect(Graphs.vertices(g))
    subgraphs = collect(Combinatorics.combinations(nodes, n))
    
    found = []
    # found_g = Array{Graphs.DiGraph}(undef,0)
    for sg in subgraphs
        subg, _ = Graphs.induced_subgraph(g, sg)
        if is_isomorphic(subg, motif)
            push!(found, sg)
            # push!(found_g, subg)
        end
    end
    return found #, found_g
end