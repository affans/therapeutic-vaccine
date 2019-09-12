

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