function cascade(Asort,kout,kin,gprim,sigma,tmax)
    
    # Asort,niche = nichemodelweb(S,C);
    n = size(Asort)[1];
    
    s = Array{Int64}(undef,tmax,n);
    #Define initial state
    s0 = rand([-1,1],n);
    s[1,:] = s0;

    G = Graphs.DiGraph(Asort);
    
    epsilondist = Normal(0,1);
    
    g = 0;

    #Build list of species' resources and predators
    nn_out = [Int[] for _ in 1:n]
    nn_in = [Int[] for _ in 1:n]
    for i=1:n
        #PREDATION
        #Predator influence is negative
        #nn_out = What predators eat species i
        nn_out[i] = Graphs.outneighbors(G,i); #NOTE: should run these outside of time loop and store them!
        
        #Resource consumption
        #Resource influence is positive
        #nn_in = What resources species i eats
        nn_in[i] = Graphs.inneighbors(G,i); #NOTE: should run these outside of time loop and store them!
    end
    
    for t=2:tmax
        
        past_s = copy(s[t-1,:]);
        future_s = zeros(Int64,n);
        
        #Determine the futures state of node i
        for i=1:n
            
            if length(nn_in[i]) == 0
                g = gprim;
            else
                g = 0;
            end
            
            #Spin function
            newstate = sign(past_s[i] - kout*sum(past_s[nn_out[i]]) + kin*sum(past_s[nn_in[i]]) + sigma*rand(epsilondist) + g);
            
            # newstate = sign(-kin*sum(past_s[nn_in]) + sigma*rand(epsilondist));
            
            future_s[i] = newstate;
        end
        
        # nn_out = outneighbors.(G,collect(1:n));
        # nn_in = inneighbors.(G,collect(1:n));
        # past_out = Array{Array}(n);
        # past_in = Array{Array}(n);
        # eprand = rand(epsilondist,n);
        # for i=1:n
        #     past_out[i] = past_s[nn_out[i]];
        #     past_in[i] = past_s[nn_in[i]];
        # end
        # #Spin function
        # newstate = sign.(-kout.*sum.(past_out) + kin.*sum.(past_in) + sigma.*eprand);
        # future_s = newstate;
        
        s[t,:] = future_s;
    end
    
    return(
    s
    )
    
end
