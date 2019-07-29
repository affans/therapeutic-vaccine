 ## male/female sexual frequency
 function distribution_sexfrequency()
    dist_men = [    [0.167, 0.334, 0.563, 0.792, 0.896, 1],     # 15 - 24    
                    [0.109, 0.572, 0.7575, 0.943, 0.9725, 1], # 25 - 29                   
                    [0.201, 0.674, 0.808, 0.942, 0.971, 1], # 30 - 39
                    [0.254, 0.764, 0.8635, 0.963, 0.9815, 1]] # 40 - 49                                       
    dist_women = [  [0.265, 0.412, 0.5885, 0.765, 0.8825, 1],     # 15 - 24    
                    [0.151, 0.628, 0.804, 0.98, 0.99, 1], # 25 - 29                   
                    [0.228, 0.73, 0.8395, 0.949, 0.9745, 1], # 30 - 39
                    [0.298, 0.764, 0.868, 0.972, 0.9855, 1]] # 40 - 49                   
    return dist_men, dist_women
end

function calculatesexfrequency(age, sex)
    ## this function calculates sex frequency based on the relevant distribution distribution
    ag = get_age_group(age)  
    mfd, wfd = distribution_sexfrequency()  ## get the distributions
    rn = rand() ## roll a dice
    sexfreq = 0
    if sex == MALE 
        sexfreq = findfirst(x -> rn <= x, mfd[ag]) - 1   #if male, use the male distribution
    else 
        sexfreq = findfirst(x -> rn <= x, wfd[ag]) - 1   #if female, use the female distribution
    end
    return sexfreq
end

function get_age_group(age::Int64)
    ## this is mainly used for sex frequency 
    if age >= 15 && age < 25
        agegroup = 1
    elseif age >= 25 && age < 30
        agegroup = 2
    elseif age >= 30 && age < 40
        agegroup = 3
    elseif age >= 40 && age < 50 
        agegroup = 4 
    end 
    return agegroup
end