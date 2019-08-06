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
    const tracking = 0

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

    function track(id)
        @debug "Tracking $id"
        tracking = id
    end
    export track  

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
        init_disease()
        # time loop 
        # for i = 1:P.sim_time
        #     ## main simulation loop started. 
        #     age()
        #     partnerup()
        # end 
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
                oldgrp = humans[h.partner].grp
                oldsex = humans[h.partner].sex  
                init_human(humans[h.partner])
                humans[h.partner].grp = oldgrp
                humans[h.partner].age = oldsex
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
        oldgrp = h.grp
        oldsex = h.sex  
        init_human(h)              
        h.age = 15   
        h.grp = oldgrp     
        h.sex = oldsex
        @debug "human exited population successfully."
    end
    export exit_population

    ## contact pairings, married pairs        
    function partnerup()
        # function assigns partners to non-married people in each age-group
        # an individual is only partnered with someone in their own age group
        
        # before starting, reset everyone's (NON MARRIED) partner. This is important for the "reshuffling" every 6 months.        
        ## NOT IMPLEMENTED: do some partners tend to stay with each other during the year?
        ## Use the parameter P.pct_partnerchange to implement this.
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
    export partnerup

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
    export marry

    function calc_prob(age, sex, grp)
        # for grp in ((15:19), (20:29), (30:39), (40:49))
        #     println(length(findall(x -> x.age ∈ grp, humans))/10000)
        # end

        if sex == MALE  
            s = :male
            sd = :male_dis
        else 
            s = :female
            sd = :female_dis
        end

        if grp == WHITE 
            g = :white
            gd = :white_dis
        elseif grp == BLACK 
            g = :black
            gd = :black_dis
        elseif grp == ASIAN
            g = :asian
            gd = :asian_dis
        elseif grp == HIS
            g = :hispanic
            gd = :hispanic_dis
        end

        if age in (15:19)
            a = :ag1
            ad = :ag1_dis
        elseif age in (20:29)
            a = :ag2
            ad = :ag2_dis
        elseif age in (30:39)
            a = :ag3
            ad = :ag3_dis
        elseif age in (40:49)
            a = :ag4
            ad = :ag4_dis
        end

        @debug "$a, $s, $g"
        @debug "$ad, $sd, $gd"

        conds = Dict{Symbol, Float64}()
        push!(conds, :dis => 0.12)
        push!(conds, :male => 0.50)
        push!(conds, :female => 0.50)
        push!(conds, :dis_male => 0.08)
        push!(conds, :dis_female => 0.16)
        push!(conds, :white => 0.65)
        push!(conds, :black => 0.12)
        push!(conds, :asian => 0.06)
        push!(conds, :hispanic => 0.17)
        push!(conds, :dis_white => 0.08)
        push!(conds, :dis_black => 0.346)
        push!(conds, :dis_asian => 0.038)
        push!(conds, :dis_hispanic => 0.094)
        push!(conds, :ag1 => 0.1354)
        push!(conds, :ag2 => 0.2393)
        push!(conds, :ag3 => 0.2434)
        push!(conds, :ag4 => 0.3819)
        push!(conds, :dis_ag1 => 0.008)
        push!(conds, :dis_ag2 => 0.076)
        push!(conds, :dis_ag3 => 0.133)
        push!(conds, :dis_ag4 => 0.212)
        conds[:male_dis] = conds[:dis_male]*conds[:male]/conds[:dis]            
        conds[:female_dis] = conds[:dis_female]*conds[:female]/conds[:dis]
        conds[:white_dis] = conds[:dis_white]*conds[:white]/conds[:dis]
        conds[:black_dis] = conds[:dis_black]*conds[:black]/conds[:dis]
        conds[:asian_dis] = conds[:dis_asian]*conds[:asian]/conds[:dis]
        conds[:hispanic_dis] = conds[:dis_hispanic]*conds[:hispanic]/conds[:dis]
        conds[:ag1_dis] = conds[:dis_ag1]*conds[:ag1]/conds[:dis]
        conds[:ag2_dis] = conds[:dis_ag2]*conds[:ag2]/conds[:dis]
        conds[:ag3_dis] = conds[:dis_ag3]*conds[:ag3]/conds[:dis]
        conds[:ag4_dis] = conds[:dis_ag4]*conds[:ag4]/conds[:dis]
    
        ## calculate the conds
        prob = (conds[ad]*conds[sd]*conds[gd]*conds[:dis])/(conds[a]*conds[s]*conds[g])
        return round(prob, digits = 4)
    end

    function init_disease()
        for i = 1:10000
            rn = rand()
            prb = calc_prob(humans[i].age, humans[i].sex, humans[i].grp)
            if rn < prb
                humans[i].health = INF
            end
        end
    end
    export init_disease

    function get_pairs()
        ## helper function which may not get used. 
        tmp = Array{Tuple{Int64, Int64}, 1}()
        for i = 1:10000
            if humans[i].partner > 0
                push!(tmp, (i, humans[i].partner))
            end            
        end
        tmp2 = unique(x -> Set(x), tmp)
        return tmp2
    end

    function get_sick_pairs()
        ## this function goes through all the pairs and returns only pairs 
        ## that have ONE OR MORE infected.
        ## it returns an array of tuples of the human/partner id. 
        tmp = Array{Tuple{Int64, Int64}, 1}()
        for i = 1:10000
            if humans[i].partner > 0 && humans[i].health == INF
                push!(tmp, (i, humans[i].partner))
            end            
        end
        tmp2 = unique(x -> Set(x), tmp)
        return tmp2
    end

    function get_single_sick_pairs()
        ## this function goes through all the pairs 
        ## and returns only pairs that have ONE infected. 
        ## if both of them are infected or susceptible, it dosn't return. 
        ## it returns an array of tuples of the human/partner id. 
        ## moreover the first element of the tuple is the SICK person and the second element is the SUSC person. 
        tmp = Array{Tuple{Int64, Int64}, 1}()
        for i = 1:10000
            if humans[i].partner > 0 && humans[i].health == INF && humans[humans[i].partner].health == SUSC
                push!(tmp, (i, humans[i].partner))
            end            
        end
        return tmp
    end
    export get_single_sick_pairs

    function get_sexual_encounters()
        ## function returns how many sexual encounters a pair will have
    end

    
    function applydis(pair)
        #go through each pair
        p1 = humans[pair[1]]
        p2 = humans[pair[2]]

        @debug p1
        @debug p2

        diseasetransfer = false

        recur_dist = p1.firstyearinfection ? Categorical(P.num_recur_firstyear) : Categorical(P.num_recur_thereafter)
        numofepisodes = rand(recur_dist)
        println("total symptomatic episodes: $numofepisodes.")
        if numofepisodes > 0 
            if p1.firstyearinfection  ## this really shouldnt make a difference.               
                durationsymp = P.duration_first[p1.sex] + P.duration_recur[p1.sex]*(numofepisodes - 1)
            else 
                durationsymp =  P.duration_recur[p1.sex]*(numofepisodes)
            end
            durationshed =  Int(round(P.pct_days_shed_symp*durationsymp / 7))
            # sample for each week of symtpomatic shedding the number of times they have sex
            numofsex = reduce(+, [calculatesexfrequency(p1.age, p1.sex) for i=1:durationshed])
                
            # check if disease will transfer in these sexual encounters. 
            for i = 1:numofsex 
                if rand() < P.beta
                    diseasetransfer = true
                end
            end
            println("total duration of symptomatic episodes in days: $durationsymp")
            println("total shedding time in weeks: $durationshed")
            println("total sexual encounters in $durationshed weeks: $numofsex")
            println("disease was successfully transferred during symptomatic periods")
        end

        ## now deal with all the days the individual is asymptmatic.
        durationasymp = 360 - (durationsymp)
        durationshed =  Int(round(P.pct_days_shed_asymp*(durationasymp) / 7))
        numofsex = reduce(+, [calculatesexfrequency(p1.age, p1.sex) for i=1:durationshed])
        for i = 1:numofsex 
            if rand() < P.beta*P.asymp_reduction
                diseasetransfer = true
            end
        end 
        println("total duration of symptomatic episodes in days: $durationasymp")
        println("total shedding time in weeks: $durationshed")
        println("total sexual encounters in $durationshed weeks: $numofsex")
        println("disease was successfully transferred during asymptomatic periods")

        ## check how many people don't develop primary infection 
        ## the numofsex while in symptomatic weeks should be less.
        ## what happens to first year infection for the partner status?
        return diseasetransfer
    end
    export applydis

end # module
