## Enums
@enum HEALTH SUSC=1 LAT=2 ASYMP=3 SYMP=4 VAC=5 REC=6 
@enum SEX MALE=1 FEMALE=2
@enum ETH WHITE=1 BLACK=2 ASIAN=3 HIS=4

## MAIN SYSTEM PARAMETER
@with_kw struct ModelParameters @deftype Int64
    # general parameters
    sim_time = 20  ## ten years in 6 month intervals. 
    num_of_humans = 10000
    percent_married::Float64 = 0.20 ## percentage of people married at the start of sims
end

mutable struct Human
    health::HEALTH
    
    ## demographics
    age::Int64
    sex::SEX # 0: female, 1:male
    eth::ETH

    numofsex::Int64
    partner::Int64
    married::Bool
    Human() = new(SUSC, rand(15:49), SEX(rand(1:2)), ETH(rand(1:4)), 0, 0, false)
end

function init_humans()
    @inbounds for i = 1:length(humans)        
        humans[i] = Human()
    end
end