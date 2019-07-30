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
