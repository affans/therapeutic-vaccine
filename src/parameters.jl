## Enums
@enum HEALTH SUSC=1 INF=2 ASYMP=3 SYMP=4 VAC=5 REC=6 
@enum SEX MALE=1 FEMALE=2
@enum GRP WHITE=1 BLACK=2 ASIAN=3 HIS=4

## MAIN SYSTEM PARAMETER
@with_kw struct ModelParameters @deftype Float64
    # general parameters
    sim_time::Int64 = 20  ## ten years in 6 month intervals. 
    num_of_humans::Int64 = 10000
   

    ## demographic information https://factfinder.census.gov/bkmk/table/1.0/en/ACS/17_5YR/DP05
    grp_hispanic = 0.17
    grp_white = 0.65
    grp_black = 0.12
    grp_asian = 0.06
    pct_married = 0.20 ## percentage of people married at the start of sims
    pct_partnerchange = 0.50 ## not implemented yet

    # disease parameters.
    pct_legions = 0.30 ## percentage of people married at the start of sims
    beta = 0.01
    asymp_reduction = 0.50
    incubation = 4.3   ## average incubation days.

    duration_first::Dict{SEX, Int64} = Dict(MALE => 17, FEMALE => 20)
    duration_recur::Dict{SEX, Int64} = Dict(MALE => 10, FEMALE => 12)

    ## number of recurrances: 0 = 11%, 1-6: 51% (so 51/6 = 8.5%)
    num_recur_firstyear::Array{Float64} = [0.11, 0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.095, 0.095, 0.095, 0.095]
    num_recur_thereafter::Array{Float64} = [0.20, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.025, 0.025, 0.025, 0.025]

    pct_days_shed_symp = 0.69
    pct_days_shed_asymp = 0.10

    
end

mutable struct Human
    health::HEALTH
    legions::Bool
    recur::Int64 ## number of recurrances in a time step (~1 year?)

    ## demographics
    age::Int64
    sex::SEX # 0: female, 1:male   
    grp::GRP

    partner::Int64
    married::Bool

    firstyearinfection::Bool # infection first year 

    Human() = new()
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
    h.legions = rand() < P.pct_legions ? true : false

    # demographics
    h.age = rand(agebraks[rand(agedist)])
    h.sex = rand() < 0.5 ? MALE : FEMALE
    h.grp = GRP(rand(grpdist))

    # partners
    h.partner = 0
    h.married = false

    h.firstyearinfection = true

    return h
end
