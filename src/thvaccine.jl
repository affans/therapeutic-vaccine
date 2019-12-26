# add distribut

module thvaccine
using Distributions
using Parameters
using Random
using DataFrames
using Distributed

import Base: show
include("./functions.jl")
#main exports.
export P, humans, main, modelinfo

const gridsize = 10000
const P = ModelParameters()
const humans = Array{Human}(undef, gridsize)

function main(simnumber=1, scenario = 1.0, cov=0.0, efficacy=0.0, warmup_beta=0.016, main_beta=0.07, 
    warmup_time=50, eql_time=100, run_time=10) 
    Random.seed!(simnumber) 
    #println(simnumber)
    ## error checks 
    warmup_beta + main_beta == 0.0 && error("β is set to zero. No disease will happen")
    
    P.sim_time = warmup_time + eql_time + run_time
    P.sim_time == 0 && error("simulation time is set to zero")

    P.treatment_coverage = cov
    P.vaccine_efficacy = efficacy

    dat = SimData(P)  # initialize data collection

    ## setup initial distribution. 
    init_population()
    create_partners()
    marry()
    init_disease()   

    ## create the right end points for the time loops
    t1 = warmup_time
    t2 = warmup_time + eql_time
    t3 = warmup_time + eql_time + run_time

    #println("t1: $t1, t2: $t2, t3: $t3")
    P.beta = warmup_beta
    for yr = 1:(t1-1)
        dat.disease[yr, 1:end] .= transmission(dat, yr)   
        record_prev(dat, yr) ## record prevalence data       
        dat.agedist[yr, [:left, :left_ct]] .= age()
        create_partners()
    end 
    
    P.beta = main_beta
    for yr = t1:(t2-1)
        dat.disease[yr, 1:end] .= transmission(dat, yr)  
        record_prev(dat, yr) ## record prevalence data     
        dat.agedist[yr, [:left, :left_ct]] .= age()
        create_partners()
    end

    ## select the right intervention function
    if scenario == 1     ## treatment 
        _func = suppressive_treatment
    elseif scenario == 2 ## vaccine
        _func = vaccine
    end
    ## give it an extra year so the counting process smoothes itself out
    for yr = (t2):t3     
        dat.disease[yr, 1:end] .= transmission(dat, yr)        
        dat.treatment[yr, :total_treated] = _func(P.treatment_coverage)
        record_prev(dat, yr) ## record data before system changes at end of year  
        dat.agedist[yr, [:left, :left_ct, :left_treated]] .= age()
        create_partners()
    end
    return dat ## return the data structure.
end

function record_prev(dat, year)
    ## this just runs some queries for data collection.
    ## make sure dat is initialized as a SimData object.
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
    
    ## enter data in dataframes
    nt = length(findall(x -> x.health == INF && x.newlyinfected == true, humans))
    dat.prevalence[year, 1:end] .=  (total, nt, ag1, ag2, ag3, ag4, wte, blk, asn, his, M, F)

    #dat.[!, :tt] = _tmpcntarr

    dat.partners[year, 1:end] .=  pair_stats()
    dat.agedist[year, 1] = length(findall(x -> get_age_group(x.age) == 1, humans))
    dat.agedist[year, 2] = length(findall(x -> get_age_group(x.age) == 2, humans))
    dat.agedist[year, 3] = length(findall(x -> get_age_group(x.age) == 3, humans))
    dat.agedist[year, 4] = length(findall(x -> get_age_group(x.age) == 4, humans))
end
export record_prev

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
    # this also set's the newlyinfected parameter back to false.       
    ct = 0          ## number of people leaving
    ct_inf = 0      ## number of people leaving that were infected
    ct_treated = 0  ## number of people leaving that were vaccinated/treated
    for h in humans 
        h.age += 1 
        h.newlyinfected = false
        if h.age > 49 
            ct += 1
            if h.health == INF 
                ct_inf += 1
            end
            if h.treated > 0 || h.vaccinated == true
                ct_treated += 1
            end
            exit_population(h)
        end
    end
    return ct, ct_inf, ct_treated
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
    ## resets all the partners except married couples
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
    sick1r = _mynaturalhistory(sick1)
    sick2r = _mynaturalhistory(sick2)

    ## return false as first element since disease is not transferred. 
    ## and to keep it consistent with _infsusc() function as well.
    return (false, sick1r, sick2r)  
end
export _infinf

function _infsusc(sick::Human, susc::Human)
    sickr = _mynaturalhistory(sick) ## run natural history of disease. 
    dt = check_transfer_disease(sick, sickr.numofsex_symp, sickr.numofsex_asymp)
    if dt 
        ## disease is passed onto the susceptible person
        susc.health = INF 
        susc.newlyinfected = true
        suscr = _mynaturalhistory(susc) ## run natural history of disease.
        return (dt, sickr, suscr) 
    end
    # if disease is not transferred, return a zero'ed object (with proper id set) even if nothing happens. 
    # initially I was returning Nothing but this way we can have two rows per inf/susc pair
    return (dt, sickr, NaturalHistory(susc.id, sick.id, [0 for i = 1:(length(fieldnames(NaturalHistory)) - 2)]...))    
end
export _infsusc

function check_transfer_disease(x, symp_t, asymp_t)
    # symp_t:  number of sexual encounters during symptomatic phase.
    # asymp_t: number of sexual encounters during asymptomatic phases. 
    # this is split incase since beta is different in each phase, dependent on the scenario 
    # oct 5th: during sexual interactions in symptomatic days, the vaccinated individual has no reduction in beta. 
    # ... in other words, their beta during symptomatic is same as calibrated beta. 
    # ... seyed has papers justifying this. 

    dt = false     
    
    sympbeta = P.beta
    asympbeta = P.beta
    ## if under suppressive treatment
    if x.treated == 1 
        sympbeta = sympbeta*(1 - 0.80)
        asympbeta = asympbeta*(1 - 0.80)
    end

    ## if under vaccinated scenario
    if x.vaccinated 
        sympbeta = sympbeta
        asympbeta = asympbeta*(1 - P.vaccine_efficacy)
    end 

    for i = 1:symp_t        
        if rand() < sympbeta
            dt = true
        end
    end    
    for i = 1:asymp_t       
        if rand() < asympbeta 
            dt = true
        end
    end 
    return dt
end
export check_transfer_disease

function unpack_naturalhistory(dt, rtuple)
    ## this function not currently being used
    ## this is a helper function that unpacks the NaturalHistory() object
    ## and creates a tuple so it can be pushed to a data frame. 
    for i in fieldnames(NaturalHistory)
        rtuple = (rtuple..., getfield(dt, i))
    end
    return rtuple
end
export unpack_naturalhistory

function transmission(dataobj::SimData, year = 1)
    ## main transmission dynamics function.    
 
    ## get the sick pairs (i.e. inf/susc or inf/inf)
    sick_pairs = get_partners(onlysick = true)
            
    ctr_xor = 0   ## total number of SUSC/INF pairs
    ctr_dis = 0   ## total number of the SUSC/INF pairs that got sick. 
    ctr_inf = 0   ## total number of INF/INF pairs.

    ## "total" counters for everyone. 
    ds = 0  ## total symptomatic days
    ss = 0  ## total shedding days
    da = 0  ## total asympotmatic days
    sa = 0  ## total shedding days

    ds_notreat = 0 
    ss_notreat = 0

    for p in sick_pairs
        p1 = humans[p.a]
        p2 = humans[p.b]
        p1health = p1.health 
        p2health = p2.health

        if p1health == INF && p2health == INF 
            ctr_inf += 1
            dt, p1nathis, p2nathis =  _infinf(p1, p2) #returns a Tuple{Bool, NaturalHistory, NaturalHistory}
            

            ## add to the statistic variables
            ds += (p1nathis.duration_symp + p2nathis.duration_symp)
            ss += (p1nathis.shedding_symp + p2nathis.shedding_symp)
            da += (p1nathis.duration_asymp + p2nathis.duration_asymp)
            sa += (p1nathis.shedding_asymp + p2nathis.shedding_asymp)

            if p1.treated == 0
                ds_notreat += p1nathis.duration_symp
                ss_notreat += p1nathis.shedding_symp
            end

            if p2.treated == 0
                ds_notreat += p2nathis.duration_symp
                ss_notreat += p2nathis.shedding_symp
            end

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
            ## if disease is not transferred, the natural history of susc is all zeros.
            dt, p1nathis, p2nathis = _infsusc(sick, susc)
            if dt ## disease has transferred
                ctr_dis += 1
            end

            ds += (p1nathis.duration_symp + p2nathis.duration_symp)
            ss += (p1nathis.shedding_symp + p2nathis.shedding_symp)
            da += (p1nathis.duration_asymp + p2nathis.duration_asymp)
            sa += (p1nathis.shedding_asymp + p2nathis.shedding_asymp)
           

            if p1.treated == 0
                ds_notreat += p1nathis.duration_symp
                ss_notreat += p1nathis.shedding_symp
            end

            if p2.treated == 0
                ds_notreat += p2nathis.duration_symp
                ss_notreat += p2nathis.shedding_symp
            end
        end
    end
    return (ctr_inf, ctr_xor, ctr_dis, ds, ss, ds_notreat, ss_notreat, da, sa)
end
export transmission

function _get_shedding_weeks(x::Human) 
    # quick calculation to see how many weeks one will shed in asympotmatic and symptmatic infections during a year. 
    # this is based of the JAMA paper, where about 90% of the swabs/days were subclinical. 
    # I used the assumption that this will be from 85-95%. 
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

    #println("sampled days symptomatic (shedding): $days_symptomatic ($shed_symptomatic)")
    #println("sampled days asymptomatic (shedding): $days_asymptomatic ($shed_asymptomatic) \n")

    ## hacky way to make sure everyone is under episodic treatment... 
    if x.treated == 2 || true  
        ## save the original number of symptomatic days
        r = days_symptomatic 

        ## there is reduction in symptomatic days (which is not used in the calculation for reduction of shedding)
        ## take the original sampled days of symptomatic and reduce by 50%. 
        ## recalculate asymptomatic days as well
        days_symptomatic = r * (1 - 0.50)
        days_asymptomatic = 365 - days_symptomatic

        ## reduce shedding from the original number of symptomatic days (with baseline shedding)
        ## shedding from asymptomatic is not reduced because this is episodic treatment (but have to recalculate because there are more asympotmatic days)
        shed_symptomatic = r * pct_shed_symptomatic #* (1 - 0.50)
        shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic

        #println("treated days (epis) symptomatic (shedding): $days_symptomatic ($shed_symptomatic)")
        #println("treated days (epis) asymptomatic (shedding): $days_asymptomatic ($shed_asymptomatic) \n")
    end

    # suppressive treatment.
    if x.treated == 1 
        ## save the original number of symptomatic days
        r = days_symptomatic 

        ## there is also reduction in symptomatic days (which is not used in the calculation for reduction of shedding)
        ## take the original sampled days of symptomatic and reduce by 50%, recalculate asymptomatic days as well
        days_symptomatic = r * (1 - 0.50)
        days_asymptomatic = 365 - days_symptomatic

        ## reduce shedding from the original number of symptomatic days (with baseline shedding)
        ## shedding from asymptomatic is reduced now by 80% as well since this is suppressive
        shed_symptomatic = r * pct_shed_symptomatic #* (1 - 0.80)
        shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic #* (1 - 0.80)
    end

    # 50% reduction in symptomatic days if individual is vaccinated.
    if x.vaccinated 
        r = days_symptomatic

        ## there is also reduction in symptomatic days (which is not used in the calculation for reduction of shedding)
        ## take the original sampled days of symptomatic and reduce by 50%, recalculate asymptomatic days as well
        days_symptomatic = r * (1 - 0.50)
        days_asymptomatic = 365 - days_symptomatic

        shed_symptomatic = r * pct_shed_symptomatic #* (1 - 0.7)
        shed_asymptomatic = days_asymptomatic * pct_shed_asymptomatic #* (1 - 0.7)
    end
  
    return round.((days_symptomatic, shed_symptomatic, days_asymptomatic, shed_asymptomatic); digits=2)
end

# 22 days -> convert to number of weeks, 4 weeks -> 2/sex per week = 8 interactions. 
# 

function _mynaturalhistory(x::Human)
    ds, ss, da, sa  = _get_shedding_weeks(x)
    
    ## using the number of shedding days, calculate the number of weeks they are shedding.
    ## if they are shedding for atleast 1 day, make it default to 1 week. 
    wa = sa > 0 ? max(1, Int(ceil(sa/7))) : 0
    ws =  ss > 0 ? max(1, Int(ceil(ss/7))) : 0 
    #println("weeks symptomatic (shedding):  $weeks_symptomatic")
    #println("weeks asymptomatic (shedding): $weeks_asymptomatic")
    
    ## calculate the number of sexual encounters during shedding weeks
    sex_encounters_asymptomatic = [calculatesexfrequency(x.age, x.sex) for i=1:wa]
    sex_encounters_symptomatic = [calculatesexfrequency(x.age, x.sex) for i=1:ws]
    
    numofsex_asymp = sum(sex_encounters_asymptomatic)  
    numofsex_symp = sum(sex_encounters_symptomatic)  

    res = NaturalHistory(x.id, x.partner, ds, ss, da, sa, numofsex_symp, numofsex_asymp)
    return res
end
export _mynaturalhistory

function suppressive_treatment(coverage)    
    ## if an individual gets sick during the year, their treatment really starts right away
    ## but in our model, it would get recorded in the following year.  
    cnt = 0
    days = 0 
    for x in humans
        if x.health == INF
            if (rand() < coverage && x.newlyinfected == true) || x.treated == 1
                x.treated = 1
                cnt += 1
            end
        end
    end
    return cnt
end
export suppressive_treatment

function vaccine(coverage)
    # this function is run at the end of year 
    # it goes through every single human, and if they are infected, 
    # they get a vaccine dose based on coverage. 
    # they do not get a dose everyyear though, only when the vaccine has waned. 

    cnt = 0 ## count of new vaccinated in the year.
    for x in humans
        if x.health == INF 
            if (rand() < coverage && x.newlyinfected == true) || x.vaccinated == true ## or this person has been vaccinated before
                x.vaccinated = true               
                cnt += 1
            end           
        end
    end
    return cnt  
end
export vaccine

end # module
