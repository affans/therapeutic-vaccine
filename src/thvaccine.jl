module thvaccine
using Distributions
using Parameters
using Random
using DataFrames
using Distributed

import Base: show
include("./parameters.jl")
include("./functions.jl")
#main exports.
export P, humans, main, modelinfo

const gridsize = 10000
const P = ModelParameters()
const humans = Array{Human}(undef, gridsize)
const verbose = false ## not used

function main(simnumber=1, warmup_beta=0.0, main_beta=0.0, 
    warmup_time=0, eql_time=0, run_time=0, includetreatment=false) 
    #Random.seed!(simnumber) 
 
    ## error checks 
    warmup_beta + main_beta == 0.0 && error("β is set to zero. No disease will happen")
    P.sim_time = warmup_time + eql_time + run_time
    P.sim_time == 0 && error("simulation time is set to zero")

    dat = SimData(P)  # initialize data collection

    ## setup initial distribution. 
    init_population()
    create_partners()
    marry()
    init_disease()
    if includetreatment 
        treatment(2, 1.0) 
    end
    # ## give all the infected people episodic treatment. 

    ## create the right end points for the time llops
    t1 = warmup_time
    t2 = warmup_time + eql_time
    t3 = warmup_time + eql_time + run_time

    _tmpcntarr = zeros(Int64, P.sim_time)## temp array to store data

    P.beta = warmup_beta
    for yr = 1:t1
        transmission(dat, yr)   
        age()                
        create_partners()
        record_data(dat, yr) ## record prevalence data        
    end 
    
    P.beta = main_beta
    for yr = (t1 + 1):t2
        transmission(dat, yr)   
        age()                
        create_partners()
        record_data(dat, yr) 
    end

    ## main scenario run loop
    P.beta = main_beta
    if P.scenario == 1     ## treatment 
        _func = treatment
    elseif P.scenario == 2 ## vaccine
        _func = vaccine
    end
    
    for yr = (t2 + 1):t3        
        transmission(dat, yr)   
        age()                
        create_partners()
        cnt = _func(1, P.treatment_coverage)
        _tmpcntarr[yr] = cnt 
        record_data(dat, yr)  
    end
    dat.gendata[!, :treated] = _tmpcntarr

    return dat ## return the data structure.
end

function record_data(dat, year)
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
    dat.partners[year, 1:end] .=  pair_stats()
    dat.agedist[year, 1] = length(findall(x -> get_age_group(x.age) == 1, humans))
    dat.agedist[year, 2] = length(findall(x -> get_age_group(x.age) == 2, humans))
    dat.agedist[year, 3] = length(findall(x -> get_age_group(x.age) == 3, humans))
    dat.agedist[year, 4] = length(findall(x -> get_age_group(x.age) == 4, humans))
end
export record_data

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

function pair_stats()
    ## this is an information function used to record data in the dataframe. 
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

function init_disease()
    ## initialize the disease based on data distribution from paper.
    ## this might not work for calibration since we are trying to get to 12% after calibration.
    cnt = 0
    a = findall(x -> x.partner > 0 && x.married == false, humans)
    r = sample(a, 100; replace = false)
    for i in r
        humans[i].health = INF
        cnt += 1
    end
    return cnt
end
export init_disease

function _infinf(sick1::Human, sick2::Human)
    sick1.health != INF && error("infected person is not infected.")
    sick2.health != INF && error("infected person is not infected.")

    sick1r = _mynaturalhistory(sick1)
    sick2r = _mynaturalhistory(sick2)

    return (false, sick1r, sick2r) 
    ## return false as first element since disease is not transferred. 
    ## and to keep it consistent with _infsusc() function as well. 
    ## ofcourse it's false... it's a inf/inf scenario.
end
export _infinf

function _infsusc(sick::Human, susc::Human)
    ## error checks
    sick.health != INF  && error("infected person is not infected.")
    susc.health != SUSC && error("susceptible person is not susceptible.")

    sickr = _mynaturalhistory(sick) ## run natural history of disease. 
    dt = check_transfer_disease(sickr.numofsex_symp, sickr.numofsex_asymp)
    if dt 
        ## disease is passed onto the susceptible person
        susc.health = INF 
        susc.firstyearinfection = true
        suscr = _mynaturalhistory(susc) ## run natural history of disease.
        return (dt, sickr, suscr) 
    end
    # if disease is not transferred, return a zero'ed object (with proper id set) even if nothing happens. 
    # initially I was returning Nothing but this way we can have two rows per inf/susc pair
    return (dt, sickr, NaturalHistory(susc.id, sick.id, [0 for i = 1:(length(fieldnames(NaturalHistory)) - 2)]...))    
end
export _infsusc

function check_transfer_disease(symp_t, asymp_t)
    # symp_t:  number of sexual encounters during symptomatic phase.
    # asymp_t: number of sexual encounters during asymptomatic phases. Beta is reduced by 50%
    dt = false 
    # check if disease will transfer in these sexual encounters for symptomatic episodes
    for i = 1:symp_t
        if rand() < 1 - (1 - (1 - P.p)*(1 - P.q)*P.beta)
            dt = true
        end
    end    
    for i = 1:asymp_t
        if rand() < 1 - (1 - (1 - P.p)*(1 - P.q)*P.beta*P.asymp_reduction)   
            dt = true
        end
    end 
    return dt
end
export check_transfer_disease

function unpack_naturalhistory(dt, rtuple)
    ## this is a helper function that unpacks the NaturalHistory() object
    ## and creates a tuple so it can be pushed to a data frame. 
    for i in fieldnames(NaturalHistory)
        rtuple = (rtuple..., getfield(dt, i))
    end
    return rtuple
end
export unpack_naturalhistory

function transmission(dataobj::SimData, year = 1)
    # for each pair
    # if inf/susc -> see if the inf will transfer disease in the whole year 
    # -- record the DiseaseInfo()
    # -- if disease is successfully transferred, the partner is now infected as well. 
    # -- record DiseaseInfo() -- it could very well be  the partner will have episodes and also recieve vaccine. 
    
    # -- these partners may split up and join other susceptibles (but since they are vaccinated nothing will happen)

    ## get the sick pairs. these could be inf/susc or inf/inf
    sick_pairs = get_partners(onlysick = true)
            
    ctr_xor = 0   ## total number of SUSC/INF pairs
    ctr_dis = 0   ## total number of the SUSC/INF pairs that got sick. 
    ctr_inf = 0   ## total number of INF/INF pairs.

    for p in sick_pairs
        p1 = humans[p.a]
        p2 = humans[p.b]
        p1health = p1.health 
        p2health = p2.health

        if p1health == INF && p2health == INF 
            ctr_inf += 1
            rdt =  _infinf(p1, p2) #returns a Tuple{Bool, NaturalHistory, NaturalHistory}
            
            ## record the data for both people. rdt[2] and rdt[3] contains the data objects to be unpacked.
            ## care to make sure p1t, p2t are properly set.
            ## unpack the results (by appending the fields to an already existing tuple)
            p1t = unpack_naturalhistory(rdt[2], (year, 1, Int(p1.sex), 0))
            p2t = unpack_naturalhistory(rdt[3], (year, 1, Int(p2.sex), 0))

            #push!(dataobj.episodes, p1t)
            #push!(dataobj.episodes, p2t)
        end 

        if xor(p1health == INF, p2health == INF)
            ctr_xor += 1
            if p1health == INF && p2health != INF
                sick = p1
                susc = p2
            end
            if p1health != INF && p2health == INF
                sick = p2
                susc = p1
            end
            rdt = _infsusc(sick, susc)
            rdt[1] && (ctr_dis += 1)
            ## record the data for both people. rdt[2] and rdt[3] contains the data objects to be unpacked.
            ## care to make sure p1t, p2t are properly set.
            ## unpack the results (by appending the fields to an already existing tuple)
            p1t = unpack_naturalhistory(rdt[2], (year, 1, Int(p1.sex), rdt[1]))
            p2t = unpack_naturalhistory(rdt[3], (year, 1, Int(p2.sex), 0))
            
            #push!(dataobj.episodes, p1t)
            #push!(dataobj.episodes, p2t)
        end
    end
    return ctr_inf, ctr_xor, ctr_dis
end
export transmission

function _get_shedding_weeks(x::Human) 
    # quick calculation to see how many weeks one will shed in asympotmatic and symptmatic infections during a year. 
    # this is based of the JAMA paper, where about 90% of the swabs/days were subclinical. 
    # i then used the assumption that this will be from 85-95%. 
    # the amount of shedding while asymptomatic is about 12% of those days. Can we put a distribution around this? 

    # to do: still have to write unit tests for this

    # percentage of the days someone is symptomatic/asymptomatic
    pct_asymptomatic = rand(85:95)/100  ## need to justify this assumption. 
    pct_symptomatic = 1 - pct_asymptomatic   
    pct_shed_asymptomatic = 0.122
    pct_shed_symptomatic = 0.689

    days_symptomatic = pct_symptomatic * 365
    days_asymptomatic = 365 - days_symptomatic    
  
    shed_symptomatic = days_symptomatic * pct_shed_symptomatic
    shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic

    # 50% reduction in symptomatic days if individual is vaccinated.
    if x.vaccinated 
        days_symptomatic = days_symptomatic * (1 - 0.50)
        days_asymptomatic = 365 - days_symptomatic
        shed_symptomatic = days_symptomatic * pct_shed_symptomatic * (1 - 0.50)
        shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic * (1 - 0.50)
    end

    # suppressive treatment.
    if x.treated == 1 
        days_symptomatic = days_symptomatic * (1 - 0.80)    ## need to justify
        days_asymptomatic = 365 - days_symptomatic
        shed_symptomatic = days_symptomatic * pct_shed_symptomatic * (1 - 0.50)
        shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic * (1 - 0.50)
    end

    # episodic treatment.
    if x.treated == 2
        days_symptomatic = days_symptomatic * (1 - 0.60)   ## need to justify
        days_asymptomatic = 365 - days_symptomatic
        shed_symptomatic = days_symptomatic * pct_shed_symptomatic * (1 - 0.50)
        shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic
    end

    weeks_asymptomatic = Int(shed_asymptomatic ÷ 7)
    weeks_symptomatic = Int(shed_symptomatic ÷ 7)

    return (weeks_asymptomatic, weeks_symptomatic)
end

function _mynaturalhistory(x::Human)
    # to do: still have to write unit tests for this
    x.health != INF && error("person is not sick") # error check
    weeks_asymptomatic, weeks_symptomatic = _get_shedding_weeks(x)
    
    sex_encounters_asymptomatic = [calculatesexfrequency(x.age, x.sex) for i=1:weeks_asymptomatic]
    sex_encounters_symptomatic = [calculatesexfrequency(x.age, x.sex) for i=1:weeks_symptomatic]
    
    numofsex_asymp = sum(sex_encounters_asymptomatic)  
    numofsex_symp = sum(sex_encounters_symptomatic)  

    res = NaturalHistory(x.id, x.partner, 0, 0, 0, 0, weeks_symptomatic, numofsex_symp, 
    weeks_asymptomatic, numofsex_asymp)
    return res
end
export _mynaturalhistory

function treatment(type::Int64, coverage)
    ## this function is run at the end of year.
    ## it goes through every single human, and if they are infected 
    ## it sets their treatment on.
    ## if an individual gets sick during the year, their treatment really starts right away
    ## but in our model, it would get recorded in the following year. 

    ## this function needs unit testing and more logic for episodic as well
    cnt = 0
    for x in humans
        if x.health == INF && x.treated == 0
            if rand() < coverage
                x.treated = type
                cnt += 1
            end
        end
    end
    return cnt
end

export treatment

function vaccine()
    # this function is run at the end of year 
    # it goes through every single human, and if they are infected
    # it sets their vaccine status on. 
    # (if their age is past the expiry, it sets their vaccine status off)
    
    cnt = 0 ## count of new vaccinated in the year.
    for x in humans
        if x.health == INF
            x.vaccinated = true
            x.vaccineexpiry = x.age + P.vac_waningtime 
            cnt += 1
        end

        if x.age == x.vaccineexpiry 
            x.vaccinated = false
            x.vaccineexpiry = 999
        end
    end
    return cnt
end
export vaccine



end # module
