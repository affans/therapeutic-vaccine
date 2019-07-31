## Enums
@enum HEALTH SUSC=1 ASYMP=3 INF=4 VAC=5 REC=6 
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
    pct_legions = 0.30 ## percentage of people married at the start of sims
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

    return h
end
