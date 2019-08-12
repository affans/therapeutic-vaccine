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
    
    const P = ModelParameters()
    const humans = Array{Human}(undef, P.num_of_humans)
    const gridsize = P.num_of_humans
    const verbose = false ## not used

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
        println("Number of people with partners: $ans (distinct pairs: $(div(ans, 2)), confirm: $(length(get_partners())) )")
        ans = length(findall(x -> x.partner > 0 && x.married == true, humans))
        println("Number of people married: $ans (distinct pairs: $(div(ans, 2)))")
        a1 = length(findall(x -> x.health == INF, humans))
        a2 = length(findall(x -> x.partner > 0 && x.health == INF, humans))
        println("Number of infected: $a1 (with partners: $a2) ")
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

        ## data frames
        df1names = Symbol.(["Total", "Ag1", "Ag2", "Ag3", "Ag4", "Wte", "Blk", "Asn", "His", "M", "F"])
        df1types = [Int64 for i = 1:length(df1names)]
        df1 = DataFrame(df1types, df1names, P.sim_time)

        init_population()
        partnerup()
        marry()
        init_disease()
        println("simulation setup success")
        #time loop 
        for i = 1:P.sim_time
            #main simulation loop started. 
            disdyn()
            age()
            partnerup()
            modelinfo()
            println("Year $i done...")
            println("\n")
            df1[i, 1:end] .= data_df1()
        end 
        return (df1)
    end

    function data_df1()
        total = length(findall(x -> x.health == INF, humans))
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
        
        return (total, ag1, ag2, ag3, ag4, wte, blk, asn, his, M, F)
    end

    function _reset()
        ## helper function forr esetting the population while debugging.
        init_population()
        partnerup()
        marry()
        init_disease()
    end 
    function init_population()    
        @inbounds for i = 1:gridsize       
            humans[i] = Human()   ## create an empty human
            humans[i].id = i
            init_human(humans[i]) ## initialize the human
        end
        @debug "population initialized"
    end


    function age() 
        ## this increases the age of every individual. If the individual is 49+ we replace with a 15 year old. 
        # if the 49 year old had a partner, that partner is now single and is available for pairing at the next shuffle.
        # if the 49 year old had a partner and they were married, both of them are replaced.       
        for x in humans 
            x.age += 1 
            x.age > 49 && exit(x)
        end
    end
    export age 

    function exit(h::Human)
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
    export exit

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
        
    function partnerup()
        # function assigns partners to non-married people in each age-group
        # an individual is only partnered with someone in their own age group
        
        # before starting, reset everyone's (NON MARRIED) partner. This is important for the "reshuffling" every 6 months.        
        ## NOT IMPLEMENTED: do some partners tend to stay with each other during the year?
        ## Use the parameter P.pct_partnerchange to implement this.
        reset = findall(x -> x.partner > 0 && x.married == false, humans)
       # @debug "Resetting partners for $(length(reset)) non married individuals"
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
        #@debug "Number of people with partners: $result (distinct pairs: $(div(result, 2)))"        
    end
    export partnerup

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
        for (x, i) in zip(humans, 1:10000)
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

    function init_disease()
        ## initialize the disease based on data distribution
        for x in humans
            rn = rand()
            prb = disease_probability(x.age, x.sex, x.grp)
            rn < prb && (x.health = INF)
        end
    end
    export init_disease

    function disdyn()
        ## we are checking disease transmission between each pair.
        ## if it's a SUSC/INF then there is chance of transmission
        ## if it's a INF/INF then it's just a matter of data collection.
        ## if its a SUSC/SUSC, then everything is irrelevant. 
        
        #dfnames = [:year, keys(nt)...]
        #df = DataFrame([Int64 for i = 1:8], dfnames, 0)
        all_pairs = get_partners(onlysick = true)

        ctr_pair = 0
        ctr_xor = 0
        ctr_dis = 0
        for p in all_pairs
            # if humans[p.a].health == INF && humans[p.b].health == INF
            #     dta, nta = _dis(p.a)
            #     dtb, ntb = _dis(p.b)                
            # end
            p1health = humans[p.a].health 
            p2health = humans[p.b].health
            ctr_pair += 1
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
                #@debug "$(p.a) and $(p.b) either or are sick"
                #@debug "p.a=$(p.a) health=$(p1health), p.b=$(p.b) health=$(p2health)"
                #@debug "$sick is sick ID and $susc is susc ID"
                dta, nta = _dis(sick)
                if dta 
                    ctr_dis += 1
                    humans[susc].health = INF 
                    #println("$(p.a) and $(p.b) either or are sick")
                    #println("p.a=$(p.a) health=$(p1health), p.b=$(p.b) health=$(p2health)")
                    #println("$sick is sick ID and $susc is susc ID")
                    #println("susc=$susc is now sick")             
                end
            end
        end
        println("total: $ctr_pair, xor'ed: $ctr_xor, diseased: $ctr_dis")
    end
    export disdyn

    function _dis(id)
        #@debug P.beta
        x = humans[id]
        dt = false ## disease transfer flag
        if x.health != INF 
            ret = (NumSympEpisodes = 0, DurationOfSympDays = 0, 
            NumOfShedWeeksSymp = 0, TotalNumOfSexSymp = 0, 
            DurationOfASympDays = 0, NumOfShedWeeksASymp = 0, NumOfSexAsymp = 0)
            return (dt, ret)
        end
       
        recur_dist = x.firstyearinfection ? Categorical(P.num_recur_firstyear) : Categorical(P.num_recur_thereafter)
        numofepisodes = rand(recur_dist)
        # println("total symptomatic episodes: $numofepisodes.")
        # println("total duration of symptomatic episodes in days: $durationasymp")
        # println("total shedding time in weeks: $durationshed")
        # println("total sexual encounters in $durationshed weeks: $numofsex")
        # println("disease was successfully transferred during asymptomatic periods")

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
        ret = (NumSympEpisodes = numofepisodes, DurationOfSympDays = durationsymp, 
               NumOfShedWeeksSymp = durationshed_symp, TotalNumOfSexSymp = numofsex_symp, 
               DurationOfASympDays = durationasymp, NumOfShedWeeksASymp = durationshed_asymp, NumOfSexAsymp = numofsex_asymp)

        return (dt, ret)
    end 

end # module
