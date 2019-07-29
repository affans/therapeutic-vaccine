

module thvaccine
    using Distributions
    using Parameters
    using Random
    include("./parameters.jl")
    include("./functions.jl")
    export P, humans, main, modelinfo
  
    const P = ModelParameters()
    const humans = Array{Human}(undef, P.num_of_humans)
    

    function modelinfo()
        @info "FEMALE = 0, MALE = 1"
        ans = length(findall(x -> x.sex == 1, humans))
        @info "Number of MALE" ans
        ans = length(findall(x -> x.sex == 0, humans))
        @info "Number of FEMALE" ans

        ans = length(findall(x -> x.partner > 0, humans))/2
        @info "Number of people with partners" ans
    end

    function main(simnumber::Int64) 
        #Random.seed!(simnumber)        
        init_humans()
        @info "Population initialized"        
    end

    function calculatesexfrequency(age, sex)
        ## this function calculates sex frequency based on the relevant distribution distribution
        ag = get_age_group(age)  
        mfd, wfd = distribution_sexfrequency()  ## get the distributions
        rn = rand() ## roll a dice
        sexfreq = 0
        if sex == MALE 
            sexfreq = findfirst(x -> rn <= x, mfd[ag]) - 1   #if male, use the male distribution
        else 
            sexfreq = findfirst(x -> rn <= x, wfd[ag]) - 1   #if female, use the female distribution
        end
        return sexfreq
    end
    
    function pair_everyone()
        # function assigns partners to non-married people in each age-group
        # an individual is only partnered with someone in their own age group
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

end # module
