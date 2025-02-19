function detectcascade(Asort,s,tlsort)
    
    G = DiGraph(Asort');
    
    #1) find all chains in A
    
    #Start from set of basal nodes
    basalnodes = find(iszero,round.(tlsort,2));
    
    sp = dijkstra_shortest_paths(G,basalnodes);
    
    # ep = enumerate_paths(dijkstra_shortest_paths(G,basalnodes));
    # ep[find(!iszero,length.(ep))];
    
    x=bfs_tree(G, basalnodes[30];dir=:out);
    ax = adjacency_matrix(x);
    # aax = Array(ax);
    nodesmotif = find(!iszero,ax);
    verts = union([ind2sub(ax,nodesmotif)[1];ind2sub(ax,nodesmotif)[2]]);
    verts = verts[find(x->x!=basalnodes[30],verts)];
    ep = enumerate_paths(dijkstra_shortest_paths(G,basalnodes))
    
    
    g = DiGraph(A');
    enumerate_paths(dijkstra_shortest_paths(g,302))
    
    
    rowsums = sum(aax,2);
    colsums = sum(aax,1);
    rowsumsin = find(!iszero,rowsums);
    colsumsin = find(!iszero,colsums);
    keepnodes = union(rowsumsin,colsumsin);
    motif = aax[keepnodes,keepnodes];
    # tl = trophic(motif');
    tl = TrophInd(motif')[!,:TL];
    basal = findmin(tl)[2]
    enumerate_paths(dijkstra_shortest_paths(motif,1))
    
    
    plotweb(motif',zeros(Int64,size(Asort)[1]),tl)

end
