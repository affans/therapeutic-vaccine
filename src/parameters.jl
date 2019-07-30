## Enums
@enum HEALTH SUSC=1 ASYMP=3 INF=4 VAC=5 REC=6 
@enum SEX MALE=1 FEMALE=2
@enum GRP WHITE=1 BLACK=2 ASIAN=3 HIS=4

## MAIN SYSTEM PARAMETER
@with_kw struct ModelParameters @deftype Float64
    # general parameters
    sim_time::Int64 = 20  ## ten years in 6 month intervals. 
    num_of_humans::Int64 = 10000
    percent_married = 0.20 ## percentage of people married at the start of sims

    ## demographic information https://factfinder.census.gov/bkmk/table/1.0/en/ACS/17_5YR/DP05
    grp_hispanic = 0.17
    grp_white = 0.63
    grp_black = 0.12
    grp_asian = 0.06
end

mutable struct Human
    health::HEALTH
    
    ## demographics
    age::Int64
    sex::SEX # 0: female, 1:male   
    grp::ETH

    partner::Int64
    married::Bool
    Human() = new(SUSC, 0, MALE, WHITE, NONHIS, 0, false)
end

function init_humans()
    
    ## all demographic data https://factfinder.census.gov/bkmk/table/1.0/en/ACS/10_5YR/DP05
    ## hispanic is not a race. All the race categories add to 1, so how to distribute the hispanics?
    
    ## everything is adjusted for the races we have. 
    ## For example, the census reports other races as well. 
    ## to calculate our distribution, I summed up the population sizes only from the races we have, and used that as the denominator. 
    
    agedist = Categorical([0.130156576, 0.124754688, 0.236352852, 0.248203293, 0.260532591])
    agebraks = [15:19, 20:24, 25:34, 35:44, 45:49]
   
    grpdist = Categorical([P.grp_white, P.grp_black, P.grp_asian, P.grp_hispanic])
    @inbounds for i = 1:length(humans)        
        humans[i] = Human()
        humans[i].age = rand(agebraks[rand(agedist)])
        humans[i].sex = rand() < 0.5 ? MALE : FEMALE
        humans[i].race = ETH(rand(grpdist))
    end
end
