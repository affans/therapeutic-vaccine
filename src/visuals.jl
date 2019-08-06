function get_partnerships()
    ## returns two vectors of the same size, where each element corresponds to the ID of the partner. 
    ## can be used in a scatter plot. 
    totalpartners = findall(x -> x.partner > 0, humans)
    mh = Vector{Int64}(undef, 0)
    fh = Vector{Int64}(undef, 0)
    
    for p in totalpartners
        if humans[p].sex == MALE 
            push!(mh, humans[p].partner)                                   
        else 
            push!(fh, humans[p].partner)                                   
        end             
    end
    return mh, fh
end




function ds()
    ## this function is run once per time step.
    ## get all the paired partners first.

    allpartners = findall(x -> x.partner > 0, humans)
    while length(allpartners) > 0 
        p1 = humans[first(pairs)]
        p2 = humans[p1.partner]
        p1idx = 1  
        p2idx = findfirst(x -> x == p1.partner, pairs)
    end


 
    unique_tmp = unique(x -> Set(x), tmp) ## get the unique pairs.
    sick_with_partners = map(unique_tmp) do x
        if humans[x[1]] == INF 
            return x[1]
    end

    while length(unique_tmp) > 0 
        x = popfirst!(unique_tmp)
        if         

    end

    ctl = true
    while length(pairs) > 0 && ctl == true
        ## get the two humans
        p1 = humans[first(pairs)]
        p2 = humans[p1.partner]
        ## get their IDs IN THE PAIRS ARRAY. NOT IN THE HUMANS ARRAY!!
        p1idx = 1  
        p2idx = findfirst(x -> x == p1.partner, pairs)
        @debug "p1 has partner: $(p1.partner), p2 has partner: $(p2.partner)"
        @debug "p1idx: $p1idx, p2idx: $p2idx"
        @debug length(pairs)
      
        ## if one of them is sick...
        if xor(p1.health == INF, p2.health == INF)
            whois = p1.health == INF ? p1 : p2
            @debug "Person who is sick is: $whois"

        end

       
        #deleteat!(pairs, p1idx)
        #deleteat!(pairs, p2idx)
        ctl = false
    end


end


