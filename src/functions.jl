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
    treatment_coverage = 0.5
    vaccine_efficacy = 0.5  
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
    vaccinated::Bool  # would get vaccinated after first episode  
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
    h.vaccinated = false
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

function disease_probability(age, sex, grp)
    ## this function calculates the initial disease probabilities
    ## it follows the initial distribution spread in the article "Prevalence of Herpes Simplex Virus Type 1 and Type 2" by Geraldine McQuillan.
    ## for a given age, sex, group in the population, it calculates the probability this person will have disease. 
    ## the math is in the notebook. 
    if sex == MALE  
        s = :male
        sd = :male_dis
    else 
        s = :female
        sd = :female_dis
    end

    if grp == WHITE 
        g = :white
        gd = :white_dis
    elseif grp == BLACK 
        g = :black
        gd = :black_dis
    elseif grp == ASIAN
        g = :asian
        gd = :asian_dis
    elseif grp == HIS
        g = :hispanic
        gd = :hispanic_dis
    end

    if age in (15:19)
        a = :ag1
        ad = :ag1_dis
    elseif age in (20:29)
        a = :ag2
        ad = :ag2_dis
    elseif age in (30:39)
        a = :ag3
        ad = :ag3_dis
    elseif age in (40:49)
        a = :ag4
        ad = :ag4_dis
    end
    
    conds = Dict{Symbol, Float64}()
    push!(conds, :dis => 0.12)
    push!(conds, :male => 0.50)
    push!(conds, :female => 0.50)
    push!(conds, :dis_male => 0.08)
    push!(conds, :dis_female => 0.16)
    push!(conds, :white => 0.65)
    push!(conds, :black => 0.12)
    push!(conds, :asian => 0.06)
    push!(conds, :hispanic => 0.17)
    push!(conds, :dis_white => 0.08)
    push!(conds, :dis_black => 0.346)
    push!(conds, :dis_asian => 0.038)
    push!(conds, :dis_hispanic => 0.094)
    push!(conds, :ag1 => 0.1354)
    push!(conds, :ag2 => 0.2393)
    push!(conds, :ag3 => 0.2434)
    push!(conds, :ag4 => 0.3819)
    push!(conds, :dis_ag1 => 0.008)
    push!(conds, :dis_ag2 => 0.076)
    push!(conds, :dis_ag3 => 0.133)
    push!(conds, :dis_ag4 => 0.212)
    conds[:male_dis] = conds[:dis_male]*conds[:male]/conds[:dis]            
    conds[:female_dis] = conds[:dis_female]*conds[:female]/conds[:dis]
    conds[:white_dis] = conds[:dis_white]*conds[:white]/conds[:dis]
    conds[:black_dis] = conds[:dis_black]*conds[:black]/conds[:dis]
    conds[:asian_dis] = conds[:dis_asian]*conds[:asian]/conds[:dis]
    conds[:hispanic_dis] = conds[:dis_hispanic]*conds[:hispanic]/conds[:dis]
    conds[:ag1_dis] = conds[:dis_ag1]*conds[:ag1]/conds[:dis]
    conds[:ag2_dis] = conds[:dis_ag2]*conds[:ag2]/conds[:dis]
    conds[:ag3_dis] = conds[:dis_ag3]*conds[:ag3]/conds[:dis]
    conds[:ag4_dis] = conds[:dis_ag4]*conds[:ag4]/conds[:dis]

    ## calculate the conds
    prob = (conds[ad]*conds[sd]*conds[gd]*conds[:dis])/(conds[a]*conds[s]*conds[g])
    return round(prob, digits = 4)
end

function calculatesexfrequency(age, sex)
    ## this function calculates sex frequency based on the relevant distribution distribution
    ag = get_age_group(age)  
    dist_men = [[0.167, 0.167, 0.229, 0.229, 0.104, 0.104],
        [0.109, 0.463, 0.1855, 0.1855, 0.0295, 0.0275],
        [0.201, 0.473, 0.134, 0.134, 0.029, 0.029],
        [0.254, 0.51, 0.0995, 0.0995, 0.0185, 0.0185]]
            
    dist_women = [[0.265, 0.147, 0.1765, 0.1765, 0.1175, 0.1175],
        [0.151, 0.477, 0.176, 0.176, 0.01, 0.01],
        [0.228, 0.502, 0.1095, 0.1095, 0.0255, 0.0255],
        [0.298, 0.466, 0.104, 0.104, 0.0135, 0.0145]]
    
    if sex == MALE 
        sexfreq = rand(Categorical(dist_men[ag])) - 1
    else 
        sexfreq = rand(Categorical(dist_women[ag])) - 1   #if female, use the female distribution
    end
    return sexfreq
end


@inline function get_age_group(age::Int64)
    ## this is mainly used for sex frequency 
    if age >= 15 && age < 25
        agegroup = 1
    elseif age >= 25 && age < 30
        agegroup = 2
    elseif age >= 30 && age < 40
        agegroup = 3
    elseif age >= 40 && age < 50 
        agegroup = 4 
    else 
        error("in get_age_group(): age out of bounds")
    end 
    return agegroup
end

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

function fs()
    findfirst(x -> x.health == INF, humans)
end
export fs