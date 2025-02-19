function motif_chain(S)
    
    #Trophic chain
    Atc = zeros(Bool,S,S);
    Atc[diagind(Atc,-1)] .= true;
    Atc_t = Atc'
    niche = collect(1:S);
    
    return Atc_t
end
