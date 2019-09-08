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
    pct_married = 0.20 ## percentage of people married at the start of sims
    pct_partnerchange = 0.50 ## not implemented yet

    # disease parameters.
    
    beta = 1.0 #0.01
    asymp_reduction = 1.0 # 0.50
    incubation = 4.3   ## average incubation days.

    ## if these change to a distribution, we have to make the changes in runtests.jl as well otherwise those tests will fail
    duration_first::Dict{SEX, Int64} = Dict(MALE => 17, FEMALE => 20)
    duration_recur::Dict{SEX, Int64} = Dict(MALE => 10, FEMALE => 12)

    ## number of recurrances: 0 = 11%, 1-6: 51% (so 51/6 = 8.5%)
    num_recur_firstyear::Array{Float64} = [0.11, 0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.095, 0.095, 0.095, 0.095]
    num_recur_thereafter::Array{Float64} = [0.20, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.025, 0.025, 0.025, 0.025, 0.04]
    ## the 0.04 at the end there added by me.... 

    ## percentage of shedding in symptomatic episodes with/without legions.
    pct_legions = 0.30 ## percentage of an episode having legions
    pct_shed_legions = 0.69  
    pct_shed_nolegions = 0.12


    vaccine_on::Bool = false
    vac_efficacy = 1.0
    vac_waningtime::Int64 = 5  ## how long vaccine provides efficacy. 

    initialsymptomatic = 0.82 ## percentage of people that have had an initial symptomatic period at start of infection. from JAMA paper according to their definitions. 
    initialepisode = 0.20  ## what is the chance of developing the first episode.
end

mutable struct Human
    id::Int64
    health::HEALTH
    
    ## demographics
    age::Int64
    sex::SEX # 0: female, 1:male   
    grp::GRP

    ## partnership and pairings
    partner::Int64
    married::Bool

    ## whether infection happens in first year.
    firstyearinfection::Bool # infection first year
    vaccinated::Bool  # would get vaccinated after first episode
    vaccineexpiry::Int64

    #
    hadfirstepisode::Bool


    Human() = new()
end

struct NaturalHistory
    id::Int64
    pd::Int64
    
    vaccinated::Int64     ## whether the individual was vaccinated
    
    numofepisodes::Int64  ## if vaccine is turned on, the "numofepisodes" and "numoflegions" woulnd't change.
    numoflegions::Int64   ## -- however the effect of vaccine is reflected in "duration_symp" 
    duration_symp::Int64
    shedding_symp::Int64
    numofsex_symp::Int64

    #duration_asymp::Int64
    shedding_asymp::Int64
    numofsex_asymp::Int64
end

struct SimData
    # data frames
    prevalence::DataFrame
    partners::DataFrame
    episodes::DataFrame

    function SimData(P)
        ## set up dataframes. when setting up data frames for yearly level data, add 1 to sim_time for the initial year

        ## setup prevalence dataframe
        _names = Symbol.(["Total","TotalPartners", "Ag1", "Ag2", "Ag3", "Ag4", "Wte", "Blk", "Asn", "His", "M", "F"])
        prev = DataFrame([Int64 for i = 1:length(_names)], _names, P.sim_time)
        prev .= 0

        partners = DataFrame([Int64 for i = 1:3], [:partners, :partners_sick, :ctr_xor], P.sim_time)
        partners .= 0
        
        episodes = DataFrame([Int64 for i = 1:14], [:year, :type, :sex, :dt, fieldnames(NaturalHistory)...])
        
        new(prev, partners, episodes)
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
    agedist = Categorical([0.13, 0.12, 0.24, 0.25, 0.26])
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

    ## if they get infected, then it's going to be their first year of infection. 
    h.firstyearinfection = true
    h.vaccinated = false
    h.vaccineexpiry = 999

    h.hadfirstepisode = false
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