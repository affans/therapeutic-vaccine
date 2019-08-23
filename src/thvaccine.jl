### before camping (after push to github so SAVE then delete these comments)
### basically added a new struct containing data frames in parameters.jl
### In this struct add any dataframe for collecting data. 
### Sometimes easier to preallocate, sometimes easier to just push

### In terms of disease dynamics.. In disyn() if INF/SUSC makes susc sick, then I run _dis for the susceptible also for the data collection
### TO DO: FIRST THING.. run _dis for INF/INF pairing as well. This is just data collection. 
### It dosn't change the dynamics (well except that it may turn `firstyearinfection` off)
### MAKE NOTE OF THIS IN GITHUB ISSUES.

## TO DO Remove Legions from human structure

module thvaccine
    using Distributions
    using Parameters
    using Random
    using DataFrames

    import Base: show
    include("./parameters.jl")
    include("./functions.jl")
    #main exports.
    export P, humans, main, modelinfo
    
    const gridsize = 10000
    const P = ModelParameters()
    const humans = Array{Human}(undef, gridsize)
    const verbose = false ## not used
    
    main() = main(1)
    function main(simnumber::Int64) 
        #Random.seed!(simnumber) 
        
        global dat = SimData(P) ## we can't use const here at the global level since each simulation needs to be on its own
        global yr ## make the year available to all functions
        #show(dat.prevalence) #this shows that at every run of the function main(), the data is a new object
        # setup initial simulation
        init_population()
        create_partners()
        marry()
        init_disease()
        
        yr = 1                ## record the initial data
        record_prevalence(yr) ## record prevalence data

        for i = 1:P.sim_time # step through discrete time (add one to sim_time since starting at 2)
            yr = i
            record_prevalence(yr) ## record prevalence data
            transmission(yr)   
            age()                
            create_partners()
        end 
        return dat ## return the data structure.
    end

    function record_prevalence(year)
        ## this just runs some queries for data collection.
        ## make sure dat is initialized as a SimData object.
        total = length(findall(x -> x.health == INF, humans))
        totalpartners = length(findall(x -> x.partner > 0 && x.health == INF, humans))
        ag1 = length(findall(x -> x.health == INF && x.age ∈ 15:19, humans))
        ag2 = length(findall(x -> x.health == INF && x.age ∈ 20:29, humans))
        ag3 = length(findall(x -> x.health == INF && x.age ∈ 30:39, humans))
        ag4 = length(findall(x -> x.health == INF && x.age ∈ 40:49, humans))
       
        wte = length(findall(x -> x.health == INF && x.grp == WHITE, humans))
        blk = length(findall(x -> x.health == INF && x.grp == BLACK, humans))
        asn = length(findall(x -> x.health == INF && x.grp == ASIAN, humans))
        his = length(findall(x -> x.health == INF && x.grp == HIS, humans))
        
        M = length(findall(x -> x.health == INF && x.sex == MALE, humans))
        F = length(findall(x -> x.health == INF && x.sex == FEMALE, humans))
        
        ## enter data in dataframes

        dat.prevalence[year, 1:end] .=  (total, totalpartners, ag1, ag2, ag3, ag4, wte, blk, asn, his, M, F)
        dat.partners[year, 1:end] .=  partner_info()
    end
    export record_prevalence

    function init_population()    
        @inbounds for i = 1:gridsize       
            humans[i] = Human()   ## create an empty human
            init_human(humans[i], i) ## initialize the human
        end
    end
    export init_population

    function age() 
        ## this increases the age of every individual. If the individual is 49+ we replace with a 15 year old. 
        # if the 49 year old had a partner, that partner is now single and is available for pairing at the next shuffle.
        # if the 49 year old had a partner and they were married, both of them are replaced.       
        ct = 0
        for h in humans 
            h.age += 1 
            if h.age > 49 
                ct += 1
                exit_population(h)
            end
        end
        return ct
    end
    export age 

    function exit_population(h::Human)
        ## this human h is exiting the pool, reset their information
        ## if the person is married, their partner leaves as well.
        ## question: How much of the old information is saved? I.e. if a black person leaves, is it a black person coming in?
        ## if a male leaves, is it a male coming back in?
        if h.partner > 0 
            if h.married == true 
                replace_human(humans[h.partner])  
                ## find another couple to marry. 
                t = findfirst(x -> x.partner > 0 && x.married == false && x.age > 19, humans)
                if t != nothing
                    # ofcourse humans[t] is not married. When the reshuffling happens the marrieds are protected.  
                    # But just double check if the partner of humans[t] is accidently married. this should never happen.
                    humans[humans[t].partner].married == true && error("bug: one partner is married, the other is not")
                    humans[t].married = true 
                    humans[humans[t].partner].married = true
                end
            else 
                humans[h.partner].partner = 0 
                humans[h.partner].married = 0
            end
        end      
        replace_human(h)
    end


    ## Partnerships
    struct Partner
        a::Int64
        b::Int64
    end
    Base.hash(p::Partner, h::UInt) = hash(hash(p.a) + hash(p.b), h)
    function Base.isequal(p1::Partner, p2::Partner)
        return (p1.a == p2.a && p1.b == p2.b) || (p1.a == p2.b && p1.b == p2.a)
    end
    export Partner


    function reset_all_partners()
        cnt = 0 
        for x in humans
            if x.partner > 0 && x.married == false
                x.partner = 0
                cnt += 1
            end
        end
        return cnt
    end
        
    function create_partners()
        # function assigns partners to non-married people in each age-group
        # an individual is only partnered with someone in their own age group
        reset_all_partners()
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
        #@debug "Number of people with partners: $result (distinct pairs: $(div(result, 2)))"        
    end
    export create_partners

    function marry()       
        ## marries people everytime this function is called. 
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

    function get_partners(; onlysick = false)
        ## only returns susceptible/infected or infected/infected partners.
        arr = Array{Partner, 1}()
        for (x, i) in zip(humans, 1:gridsize)
            if x.partner > 0 
                if onlysick
                    if x.health == INF
                        push!(arr, Partner(i, x.partner))            
                    end
                else
                    push!(arr, Partner(i, x.partner))            
                end
            end
        end
        return unique(arr)
    end
    export get_partners

    function _resetdemo()
        ## helper function forr esetting the population while debugging.
        init_population()
        create_partners()
        marry()
    end

    function _resetdisease()
        init_disease()
    end

    function _reset()
        _resetdemo()
        _resetdisease()
    end
    export _reset, _resetdemo, _resetdisease

    function init_disease()
        ## initialize the disease based on data distribution
        cnt = 0
        for x in humans
            rn = rand()
            prb = disease_probability(x.age, x.sex, x.grp)
            if rn < prb 
                x.health = INF
                cnt += 1
            end
        end
        return cnt
    end
    export init_disease

    function pair_stats()
        all_pairs = get_partners()  
        sick_pairs = get_partners(onlysick = true)
        infsusc_pairs = 0
        for p in sick_pairs
            p1health = humans[p.a].health 
            p2health = humans[p.b].health
            if xor(p1health == INF, p2health == INF)
                infsusc_pairs += 1
            end
        end
        return length(all_pairs), length(sick_pairs), infsusc_pairs
    end
    export pair_stats


    ## HAVE TO REWRITE THE TRANSMISSION + VACCINE FUNCTIONS.
    function transmission(year = 1)
        ## we are checking disease transmission between each pair.
        ## if it's a SUSC/INF then there is chance of transmission
        ## if it's a INF/INF then it's just a matter of data collection.
        ## if its a SUSC/SUSC, then everything is irrelevant. 

        ## TO DO: In an SUSC/INF pair, if the susc get's infected... run _dis() for them also.
        
        #dfnames = [:year, keys(nt)...]
        #df = DataFrame([Int64 for i = 1:8], dfnames, 0)
        sick_pairs = get_partners(onlysick = true)
              
        ctr_xor = 0   ## total number of SUSC/INF pairs
        ctr_dis = 0   ## total number of the SUSC/INF pairs that got sick. 
        for p in sick_pairs
            p1health = humans[p.a].health 
            p2health = humans[p.b].health
            if xor(p1health == INF, p2health == INF)
                ctr_xor += 1
                if p1health == INF && p2health != INF
                    sick = p.a
                    susc = p.b 
                end
                if p1health != INF && p2health == INF
                    sick = p.b
                    susc = p.a
                end
                dta = _dis(sick, year)
                if dta 
                    ctr_dis += 1
                    humans[susc].health = INF 
                    humans[susc].firstyearinfection = true
                    _dis(susc, year)  ## nothing will happen. Just need to record the data. 
                end
            end
        end
        dat.transmission[year, 1:end] .= (ctr_xor, ctr_dis)
        #println("total: $ctr_pair, xor'ed: $ctr_xor, diseased: $ctr_dis")
        return ctr_xor, ctr_dis
    end
    export transmission

    function _dis(id, year = 1)
        #@debug P.beta
        x = humans[id]
        xpartner = humans[id].partner
        dt = false ## disease transfer flag
  
        if x.health != INF ## should not happen in the main course of the simulations.
            return dt
        end
        recur_dist = x.firstyearinfection ? Categorical(P.num_recur_firstyear) : Categorical(P.num_recur_thereafter)
        numofepisodes = rand(recur_dist) - 1  ## since julia is one based. 

        ## initialize variables 
        durationsymp = 0
        durationshed_symp =  0
        numofsex_symp = 0
        if numofepisodes > 0 
            if x.firstyearinfection  ## this really shouldnt make a difference.               
                durationsymp = P.duration_first[x.sex] + P.duration_recur[x.sex]*(numofepisodes - 1)
            else 
                durationsymp =  P.duration_recur[x.sex]*(numofepisodes)
            end
            durationshed_symp =  Int(round(P.pct_days_shed_symp*durationsymp / 7))
            # sample for each week of symtpomatic shedding the number of times they have sex
            numofsex_perweek_symp = [calculatesexfrequency(x.age, x.sex) for i=1:durationshed_symp]
            numofsex_symp = reduce(+, numofsex_perweek_symp)

            # for this sick person, this has all happened within the span of one year.
            firstyearinfection = false
        end
        ## now deal with all the days the individual is asymptmatic.
        ## this part is repeated from above... can be put in it's own function. 
        durationasymp = 365 - (durationsymp)
        durationshed_asymp =  Int(round(P.pct_days_shed_asymp*(durationasymp) / 7))

        numofsex_perweek_asymp = [calculatesexfrequency(x.age, x.sex) for i=1:durationshed_asymp]
        numofsex_asymp = reduce(+, numofsex_perweek_asymp)
        for i = 1:numofsex_asymp
            if rand() < P.beta*P.asymp_reduction
                dt = true
            end
        end 

         # check if disease will transfer in these sexual encounters. 
        for i = 1:numofsex_symp
            if rand() < P.beta
                dt = true
            end
        end
        
        ## return data for collection. 
        ret = (year, id, xpartner, dt, numofepisodes, 
                durationsymp, durationshed_symp, numofsex_symp, 
                durationasymp, durationshed_asymp, numofsex_asymp)
        push!(dat.episodes, ret)
        return dt
    end 

end # module
