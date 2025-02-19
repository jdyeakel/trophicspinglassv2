function motifs_generate()

    m1_a = motif_chain(2);
    m1_g = Graphs.DiGraph(m1_a);

    m2_a = [0 1 1 ; 0 0 0 ; 0 0 0];
    m2_g = Graphs.DiGraph(m2_a);

    m3_a = [0 0 1; 0 0 1; 0 0 0];
    m3_g = Graphs.DiGraph(m3_a);

    m4_a = motif_chain(3);
    m4_g = Graphs.DiGraph(m4_a);

    m5_a = [0 1 1; 0 0 1; 0 0 0];
    m5_g = Graphs.DiGraph(m5_a);

    m6_a = [0 1 1 0; 0 0 0 1; 0 0 0 1; 0 0 0 0];
    m6_g = Graphs.DiGraph(m6_a);

    m7_a = motif_chain(4);
    m7_g = Graphs.DiGraph(m7_a);

    m8_a = [0 0 1 1; 0 0 1 1; 0 0 0 0; 0 0 0 0];
    m8_g = Graphs.DiGraph(m8_a);

    m9_a = [0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 1; 0 0 0 0 0];
    m9_g = Graphs.DiGraph(m9_a);

    motifs = Dict(
        :m1 => m1_g,
        :m2 => m2_g,
        :m3 => m3_g,
        :m4 => m4_g,
        :m5 => m5_g,
        :m6 => m6_g,
        :m7 => m7_g,
        :m8 => m8_g,
        :m9 => m9_g
    )

    return motifs

    # return m1_g, m2_g, m3_g, m4_g, m5_g, m6_g, m7_g, m8_g, m9_g

end