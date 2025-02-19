function trophic_old(A)
    R"""
    library(MASS)  
    library(NetIndices)
    rtl<-TrophInd($(A'))
    """
    @rget rtl;
    tl = rtl[!,:TL]; #Array(undef,rtl[1]) - 1;
    return tl
end

    