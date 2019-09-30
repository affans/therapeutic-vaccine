## Enums
@enum HEALTH SUSC=1 INF=2 ASYMP=3 SYMP=4 VAC=5 REC=6 
@enum SEX MALE=1 FEMALE=2
@enum GRP WHITE=1 BLACK=2 ASIAN=3 HIS=4

## MAIN SYSTEM PARAMETER
@with_kw mutable struct ModelParameters @deftype Float64
    # general parameters
    sim_time::Int64 = 20  ## ten years in 6 month intervals. 
    
    ## demographic information https://factfinder.census.gov/bkmk/table/1.0/en/ACS/17_5YR/DP05
    grp_hispanic = 0.17
    grp_white = 0.65
    grp_black = 0.12
    grp_asian = 0.06
    pct_married = 0.20       ## percentage of people married at the start of sims
    pct_partnerchange = 0.50 ## not implemented yet
    beta = 1.0 #0.01
    asymp_reduction = 1.0    ## make sure interpretation is right away... is it (1 - reduction)*beta or (reduction*beta)
    vac_waningtime::Int64 = 5  ## how long vaccine provides efficacy. 
    scenario = 1 ## 1 = treatment, 2 = vaccine
    treatment_coverage = 1.0
end

mutable struct Human
    id::Int64
    health::HEALTH 
    ## demographics
    age::Int64
    sex::SEX # 0: female, 1:male   
    grp::GRP
    partner::Int64      ## partnership and pairings
    married::Bool
    firstyearinfection::Bool # infection first year
    vaccinated::Bool  # would get vaccinated after first episode
    vaccineexpiry::Int64
    treated::Int64 # 0 = no treatment, 1 = suppressive treatment, 2 = episodic treatment
    newlyinfected::Bool ## needed for coverage scenarios.
    Human() = new()
end

struct NaturalHistory
    id::Int64
    pd::Int64    
    duration_symp::Float64
    shedding_symp::Float64
    duration_asymp::Float64
    shedding_asymp::Float64
    numofsex_symp::Int64
    numofsex_asymp::Int64
end

struct SimData
    # data frames
    prevalence::DataFrame
    agedist::DataFrame
    partners::DataFrame
    episodes::DataFrame # don't use this as of now  
    disease::DataFrame  
    treatment::DataFrame
 
    function SimData(P)
        ## set up dataframes. when setting up data frames for yearly level data, add 1 to sim_time for the initial year

        ## setup prevalence dataframe
        _names = Symbol.(["Total", "NewInfections", "Ag1", "Ag2", "Ag3", "Ag4", "Wte", "Blk", "Asn", "His", "M", "F"])
        prev = DataFrame([Int64 for i = 1:length(_names)], _names, P.sim_time)
        prev .= 0

        partners = DataFrame([Int64 for i = 1:3], [:partners, :partners_sick, :ctr_xor], P.sim_time)
        partners .= 0
        
        episodes = DataFrame([Int64 for i = 1:12], [:year, :type, :sex, :dt, fieldnames(NaturalHistory)...])
        
        agedist = DataFrame([Int64 for i = 1:7], [:gr1, :gr2, :gr3, :gr4, :left, :left_ct, :left_treated], P.sim_time)
        agedist .= 0

        treatment = DataFrame([Int64], [:total_treated], P.sim_time)
        treatment .= 0

        disease = DataFrame([Float64 for _ = 1:9], [:ctr_inf, :ctr_xor, :ctr_dis, :ds, :ss, :ds_nt, :ss_nt, :da, :sa], P.sim_time)
        disease .= 0
        new(prev, agedist, partners, episodes, disease, treatment)
    end
end


#Base.show(io::IO, ::Type{Human}) = print(io, "this is a Human type")
Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

function init_human(h::Human, id)
    init_human(h)
    h.id = id
end

function init_human(h::Human)   
    ## all demographic data https://factfinder.census.gov/bkmk/table/1.0/en/ACS/10_5YR/DP05
    ## hispanic is not a race. All the race categories add to 1, so how to distribute the hispanics?
    
    ## everything is adjusted for the races we have. 
    ## For example, the census reports other races as well. 
    ## to calculate our distribution, I summed up the population sizes only from the races we have, and used that as the denominator. 
    #agedist = Categorical([0.13, 0.12, 0.24, 0.25, 0.26])
    agedist = Categorical([0.14, 0.14, 0.30, 0.28, 0.14])  ## based off internet pyramid https://www.populationpyramid.net/united-states-of-america/2019/
    agebraks = [15:19, 20:24, 25:34, 35:44, 45:49]
    grpdist = Categorical([P.grp_white, P.grp_black, P.grp_asian, P.grp_hispanic])
  
    h.health = SUSC
   
    # demographics -- this is where most of the allocations happen
    h.age = rand(agebraks[rand(agedist)])
    h.sex = rand() < 0.5 ? MALE : FEMALE
    h.grp = GRP(rand(grpdist))

    # partners
    h.partner = 0
    h.married = false
    h.firstyearinfection = true
    h.vaccinated = false
    h.vaccineexpiry = 999
    h.treated = 0
    h.newlyinfected = false
    return h
end

function replace_human(h::Human)
    oldid =  h.id
    oldgrp = h.grp
    oldsex = h.sex  
    init_human(h)
    h.grp = oldgrp
    h.sex = oldsex
    h.age = 15
    h.id = oldid
end

## reset population to default functions.
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