module thvaccine
    using Distributions
    using Parameters
    using Random

    import Base: show
    include("./parameters.jl")
    include("./functions.jl")
    #main exports.
    export P, humans, main, modelinfo
    
    const P = ModelParameters()
    const humans = Array{Human}(undef, P.num_of_humans)
    const gridsize = P.num_of_humans


    #Base.show(io::IO, ::Type{Human}) = print(io, "this is a Human type")
    function Base.show(io::IO, ::MIME"text/plain", z::Human)
       dump(z)
    end
    
    

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
        
        #initialization stage
        # 1) initialize the population 
        #   -> demographics, 
        #   -> partnerships, 
        #   -> marriage 
        # 2) distribute the initial disease throughout the population.

        init_population()
        partnerup()
        marry()

        # time loop 
        for i = 1:P.sim_time
            ## main simulation loop started. 
            age()
            partnerup()
        end 
    end

    function init_population()    
        @inbounds for i = 1:gridsize       
            humans[i] = Human()   ## create an empty human
            init_human(humans[i]) ## initialize the human
        end
        @debug "population initialized"
    end

    function age()
        ## this increases the age of every individual. If the individual is 49+ we replace with a 15 year old. 
        # if the 49 year old had a partner, that partner is now single and is available for pairing at the next shuffle.
        # if the 49 year old had a partner and they were married, both of them are replaced. 
        howmany = 0
        for i = 1:gridsize
            humans[i].age += 1            
            if humans[i].age > 49                
                exit_population(humans[i])
            end
        end
    end

    function exit_population(h::Human)
        ## this human h is exiting the pool, reset their information
        ## question: How much of the old information is saved? I.e. if a black person leaves, is it a black person coming in?
        ## if a male leaves, is it a male coming back in?
        if h.partner > 0 
            if h.married == true 
                init_human(humans[h.partner])
                humans[h.partner].age = 15   

                ## find another couple to marry. 
                t = findfirst(x -> x.partner > 0 && x.married == false && x.age > 19, humans)
                if t != nothing
                    # ofcourse humans[t] is not married (the filter above). But just double check if the partner of humans[t] is accidently married. this should never happen.
                    humans[humans[t].partner].married == true && error("bug: one partner is married, the other is not")
                    humans[t].married = true 
                    humans[humans[t].partner].married = true
                end

            else 
                humans[h.partner].partner = 0 
                humans[h.partner].married = 0
            end
        end        
        init_human(h)              
        h.age = 15        
        @debug "human exited population successfully."
    end

    ## contact pairings, married pairs        
    function partnerup()
        # function assigns partners to non-married people in each age-group
        # an individual is only partnered with someone in their own age group
        
        # before starting, reset everyone's (NON MARRIED) partner. This is important for the "reshuffling" every 6 months.        
        ## NOT IMPLEMENTED: do some partners tend to stay with each other during the year?
        reset = findall(x -> x.partner > 0 && x.married == false, humans)
        @debug "Resetting partners for $(length(reset)) individuals"
        map(x -> humans[x].partner = 0, reset)

        for eg in (WHITE, BLACK, ASIAN, HIS)
            for ag in (15:19, 20:24, 25:29, 30:34, 35:39, 40:44, 45:49)
                ## get the indices of all the eligible males and females. 
                ## filters: sex, age, ethnic group. 
                ## married = false makes sure we don't reassign partners to married individuals 
                ## NOT IMPLEMENTED: partner > 0 makes sure we don't reassign some of the partners (out of those not married)
                malein = findall(x -> x.sex == MALE && x.age ∈ ag && x.married == false && x.grp == eg, humans)
                femalein = findall(x -> x.sex == FEMALE && x.age ∈ ag && x.married == false && x.grp == eg, humans)
    
                shuffle!(malein)
                shuffle!(femalein)
    
                for (m, f) in zip(malein, femalein)
                    #@debug "pairing male $m (age: $(humans[m].age)), female $f (age: $(humans[f].age))"
                    humans[m].partner = f
                    humans[f].partner = m
                end       
            end       
        end
       
        result = length(findall(x -> x.partner > 0, humans))
        @debug "Number of people with partners: $result (distinct pairs: $(div(result, 2)))"        
    end

    function marry()       
        ## create issue. People only the age of 19+ are married. 
        h = findall(x -> x.partner > 0 && x.married == false && x.age > 19, humans)
        howmany = Int(round(length(h)*P.pct_married))
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
        ans = length(findall(x -> x.partner > 0 && x.married == true, humans))
        @debug "Number of people married: $ans (distinct pairs: $(div(ans, 2)))"
    end

end # module
