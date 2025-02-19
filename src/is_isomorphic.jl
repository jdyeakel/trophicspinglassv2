function is_isomorphic(g1::DiGraph, g2::DiGraph)
    # Check if the number of vertices and edges are the same
    if nv(g1) != nv(g2) || ne(g1) != ne(g2)
        return false
    end
    
    # Get adjacency matrices
    adj1 = adjacency_matrix(g1)
    adj2 = adjacency_matrix(g2)
    
    # Check all permutations of vertices
    for perm in permutations(1:nv(g1))
        permuted_adj = adj1[perm, perm]
        if permuted_adj == adj2
            return true
        end
    end
    return false
end