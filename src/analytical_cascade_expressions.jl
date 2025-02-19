
# function tl_to_niche(tl,C)
#     #Our tl_to_niche function can't handle trophic levels below 1
#     #Treat them the same as nichevalue = 0
#     if tl > 1
#         nichevalue = (1/C)*((tl - 1)/11.6)^(1/0.49)
#     else
#         nichevalue = 0.0
#     end
#     nichevalue_clamped = clamp(nichevalue,0.0,1.0)
#     return nichevalue_clamped
# end

# Define analytical expressions:
# function pr_state_neg(tl,C)

#     nichevalue = tl_to_niche(tl,C)
    
#     # SPECIFICALLY for kout=10, kin=5, gprim=10
#     if C == 0.02
#         m0 = 0.355
#         m1 = 0.092
#         m2 = 0.552
#         m3 = 18.8
#     elseif C == 0.1
#         m0 = 0.4
#         m1 = 0.097
#         m2 = 0.218
#         m3 = 18.8
#     else
#         println("No parameterization is found for C=",C)
#     end
    

#     probneg = m0 + m1*(1 - exp(m2 - m3*nichevalue))

#     return probneg

# end

function tl_to_niche(tl,C)
    
    # SPECIFICALLY for kout=10, kin=5, gprim=10

    #Our tl_to_niche function can't handle trophic levels below 1
    #Treat them the same as nichevalue = 0
    if C == 0.02
        A = 1.00
        B = 1.956
        alpha = 0.74
    elseif C == 0.1
        A = 0.76
        B = 3.921
        alpha = 0.406
    end

    if tl > 1
        nichevalue = ((tl - A)/B)^(1/alpha)
    else
        nichevalue = 0.0
    end

    nichevalue_clamped = clamp(nichevalue,0.0,1.0)

    return nichevalue_clamped
end

function pr_state_neg(tl,C)

    # SPECIFICALLY for kout=10, kin=5, gprim=10
    if C == 0.02
        alpha = -0.8097
        beta = 2.372
        mu = 0.972
        nuparam = -18.209
        A = exp(alpha)
        B = exp(beta)
        nu = exp(nuparam)
    elseif C == 0.1
        alpha = -0.697
        beta = 1.434
        mu = 0.966
        nuparam = -0.053
        A = exp(alpha)
        B = exp(beta)
        nu = exp(nuparam)
    else
        println("No parameterization is found for C=",C)
    end

    probneg = A * (1 + nu * exp.(-B * (tl - mu)))^(-1/nu)

    return probneg

end

function pr_state_negfit(tl,C,fit_trophic,fit_propneg)

    if C == 0.02
        if tl < minimum(fit_trophic[1])
            closest_propneg = 0.
        else 
            # Find the closest position and value
            diffs = abs.(fit_trophic[1] .- tl)  # Broadcasted element-wise difference
            min_diff, closest_pos = findmin(diffs)
            closest_propneg = fit_propneg[1][closest_pos]
        end

    elseif C == 0.1
        if tl < minimum(fit_trophic[1])
            closest_propneg = 0.
        else
            # Find the closest position and value
            diffs = abs.(fit_trophic[2] .- tl)  # Broadcasted element-wise difference
            min_diff, closest_pos = findmin(diffs)
            closest_propneg = fit_propneg[2][closest_pos]
        end

    end

    return closest_propneg

end

function prob_state(DeltaSB,DeltaST)
    
    probSneg = 1 - (1 + exp(-(DeltaSB - DeltaST)))^(-1)

    return probSneg
end


function exp_cascaderatio_conditional(tl,S,C,tldiff,kin,kout,g=0)

    if tl > 1
        gprim = 0.;
    else
        gprim = g;
    end

    nichevalue = tl_to_niche(tl,C)

    DeltaSB = kin*2*(S-1)*C*nichevalue*(1 - 2*pr_state_neg(tl-tldiff.tdprey,C)) + gprim
    
    DeltaST_PredPos = kout*(1 + (2*(S-1)*C*(1-nichevalue)-1)*(1 - 2*pr_state_neg(tl+tldiff.tdpred,C)))
    
    DeltaST_PredNeg = kout*(-1 + (2*(S-1)*C*(1-nichevalue)-1)*(1 - 2*pr_state_neg(tl+tldiff.tdpred,C)))

    ProbSNeg_PredPos = prob_state(DeltaSB,DeltaST_PredPos)
    
    ProbSNeg_PredNeg = prob_state(DeltaSB,DeltaST_PredNeg)

    index_num = 2 + (1)*(1 - ProbSNeg_PredPos) + (-1)*ProbSNeg_PredPos

    index_denom = 2 + (1)*(1 - ProbSNeg_PredNeg) + (-1)*ProbSNeg_PredNeg

    index = log(index_num / index_denom)

    return index

end


function exp_cascaderatio_conditionalfit(tl,S,C,tldiff,kin,kout,fit_trophic,fit_propneg,g=0)

    if tl > 1
        gprim = 0.;
    else
        gprim = g;
    end

    nichevalue = tl_to_niche(tl,C)

    DeltaSB = kin*2*(S-1)*C*nichevalue*(1 - 2*pr_state_negfit(tl-tldiff.tdprey,C,fit_trophic,fit_propneg)) + gprim
    
    DeltaST_PredPos = kout*(1 + (2*(S-1)*C*(1-nichevalue)-1)*(1 - 2*pr_state_negfit(tl+tldiff.tdpred,C,fit_trophic,fit_propneg)))
    
    DeltaST_PredNeg = kout*(-1 + (2*(S-1)*C*(1-nichevalue)-1)*(1 - 2*pr_state_negfit(tl+tldiff.tdpred,C,fit_trophic,fit_propneg)))

    ProbSNeg_PredPos = prob_state(DeltaSB,DeltaST_PredPos)
    
    ProbSNeg_PredNeg = prob_state(DeltaSB,DeltaST_PredNeg)

    index_num = 2 + (1)*(1 - ProbSNeg_PredPos) + (-1)*ProbSNeg_PredPos

    index_denom = 2 + (1)*(1 - ProbSNeg_PredNeg) + (-1)*ProbSNeg_PredNeg

    index = log(index_num / index_denom)

    return index

end
