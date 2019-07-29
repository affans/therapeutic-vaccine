module thvaccine
    using Distributions
    using Parameters
    using Random
    include("./parameters.jl")
    include("./functions.jl")
    export P, humans, main, modelinfo
  
    const P = ModelParameters()
    const humans = Array{Human}(undef, P.num_of_humans)
    #const runsummary = Dict{String, Int64}()
    function modelinfo()        
        ans = length(findall(x -> x.sex == MALE, humans))
        println("Number of MALE: $ans")
        ans = length(findall(x -> x.sex == FEMALE, humans))
        println("Number of FEMALE: $ans")
        ans = length(findall(x -> x.partner > 0, humans))
        println("Number of people with partners: $ans (distinct pairs: $(div(ans, 2)))")
        ans = length(findall(x -> x.partner > 0 && x.married == true, humans))
        println("Number of people married: $ans (distinct pairs: $(div(ans, 2)))")

    end

    main() = main(1)
    function main(simnumber::Int64) 
        #Random.seed!(simnumber)        
        init_humans(); @info "Population initialized"          
        pair_everyone(); @info "Population paired up"
        marry(); @info "Pairs married"              
    end
    
    function pair_everyone()
        # function assigns partners to non-married people in each age-group
        # an individual is only partnered with someone in their own age group

        # TO DO: add logic in for keeping a certain percentage of partners fixed on top of the marriage.
        # for example, there is a 20% chance the partners remain the same. 

        # before starting, reset everyone's (NON MARRIED) partner. This is important for the "reshuffling" every 6 months.        
        reset = findall(x -> x.partner > 0 && x.married == false, humans)
        @debug "Resetting partners for $(length(reset)) individuals"
        map(x -> humans[x].partner = 0, reset)

        for ag in (1, 2, 3, 4)
            ## get the indices of all the eligible males and females
            malein = findall(x -> x.sex == MALE && get_age_group(x.age) == ag && x.married == false, humans)
            femalein = findall(x -> x.sex == FEMALE && get_age_group(x.age) == ag && x.married == false, humans)
    
            shuffle!(malein)
            shuffle!(femalein)
    
            for (m, f) in zip(malein, femalein)
                @debug "pairing male $m (age: $(humans[m].age)), female $f (age: $(humans[f].age))"
                humans[m].partner = f
                humans[f].partner = m
            end       
        end       
    end

    function marry()        
        h = findall(x -> x.partner > 0 && x.married == false, humans)
        howmany = Int(round(length(h)*P.percent_married))
        @debug "Number of people getting married" howmany
        ctr = 1
        while ctr <= howmany
            rn = rand(h)
            if humans[rn].married == 0 || humans[humans[rn].partner].married == 0
                humans[rn].married = 1
                humans[humans[rn].partner].married = 1
                ctr += 1
            end
        end
    end

    function visualpairs()
        ## allows to run a scatter plot to see the pairings. 
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

end # module
