using Test, thvaccine, Statistics

const th = thvaccine
const gridsize = length(thvaccine.humans)

@testset "DEMO" begin
    ## init_human test
    x = th.Human()
    th.init_human(x, 1)  ## init human with ID
    @test x.id == 1 && x.health == th.SUSC
    @test x.age > 0
    @test x.age >= 15 && x.age <= 49
    @test x.partner == 0 && x.married == false
    @test x.vaccinated == false && x.treated == false && x.newlyinfected == false
    
    ## replace the human and check if certain properties are set
    oldid =  x.id
    oldgrp = x.grp
    oldsex = x.sex  
    th.replace_human(x)  
    @test x.id == oldid && x.health == th.SUSC
    @test x.age > 0 
    @test x.age == 15 && x.sex == oldsex && x.grp == oldgrp 
    @test x.partner == 0 && x.married == false
    @test x.vaccinated == false && x.treated == false && x.newlyinfected == false
    

    ## check change of ID
    th.init_human(x, 2)  ## init human with ID
    @test x.id == 2

    ## now test the array of population
    th.init_population()
    @test length(findall(x -> x == undef, humans)) == 0 ## check for proper initialization
    for i = 1:gridsize
        x = humans[i]
        @test x.id == i && x.health == th.SUSC
        @test x.age > 0
        @test x.age >= 15 && x.age <= 49
        @test x.partner == 0 && x.married == false
        @test x.vaccinated == false && x.treated == false && x.newlyinfected == false
    end
    
    ## test the distribution of demographics
    # m = length(findall(x -> x.sex == MALE, humans))/10000
    # f = length(findall(x -> x.sex == FEMALE, humans))/10000
    # @test abs(m - f) <= 1.0 

    ## testing age() and exit_population() functions
    ## this equires the create_partners()/marry() functions (=> tested in the next section)
    th.init_population()
    humans[1].age = 15
    th.age()
    @test humans[1].age == 16
    ## internally tests exit_population
    humans[1].age = 49
    th.age()
    @test humans[1].age == 15 ## if the human leaves, a 15 year old is added. 

    ## test exit_population() for partners
    th.init_population()
    
    # test when partnerd, by not married. human[1] should be new. human[2] should remain
    humans[1].partner = 2
    humans[1].married = false
    humans[2].partner = 1
    humans[2].partner = false 
  
    th.exit_population(humans[1]) 
    @test humans[1].age == 15
    @test humans[1].partner == 0
    @test humans[1].married == 0
    
    @test humans[2].partner == 0
    @test humans[2].married == 0

    ## test when partnered and  married.
    humans[3].partner = 4
    humans[3].married = true
    humans[4].partner = 3
    humans[4].married = true 
    th.exit_population(humans[4])
    @test humans[3].age == 15
    @test humans[3].partner == 0
    @test humans[3].married == 0    
    @test humans[4].age == 15
    @test humans[4].partner == 0
    @test humans[4].married == 0    
end

@testset "PART" begin
    ## check whether the struct for unique partners works.
    a = th.Partner(1, 2)
    b = th.Partner(2, 1)
    @test isequal(a, b) == true
   
    ## check if there are an even number of partners, and if all partners are mutual (a1)
    th.init_population()
    howmany = th.create_partners()
    allpartners = findall(x -> x.partner > 0, humans)
    @test length(allpartners) > 0   ## check if the function even worked
    @test howmany == length(allpartners)
    @test mod(length(allpartners), 2) == 0 ## test if even number of partners
    # check if pairings are logically okay
    a1 = findall(a -> humans[humans[a].partner].partner != a, allpartners)
    @test length(a1) == 0 
    
    # check if all partners are within their age groups
    a2 = findall(allpartners) do h
        age1 = humans[h].age
        age2 = humans[humans[h].partner].age
        diffage = abs(age1 - age2) 
        diffage > 5  
    end
    @test length(a2) == 0 

    # check if they are all within their ethnicites/subgroups
    a2 = findall(allpartners) do x
        grp1 = humans[x].grp
        grp2 = humans[humans[x].partner].grp
        grp1 != grp2
    end
    @test length(a2) == 0

    # no one should be married since we havn't run the marry function yet.
    a3 = findall(x -> x.married == true, humans)
    @test length(a3) == 0

    # check if the reset_all_partners() function works
    ## all partners should be 0 unless they are married
    for x in humans  ## turn of marriage (incase marry() was already run)
        x.married = false
    end
    cw = th.reset_all_partners()
    @test cw == length(allpartners)
    @test length(findall(x -> x.partner > 0, humans)) == 0

    ## do the opposite where everyone is married and try to reset
    ## no one should be reset since they are all married
    for x in humans
        x.partner = 99
        x.married = true
    end
    cw = th.reset_all_partners()
    @test cw == 0
    @test length(findall(x -> x.partner == 0, humans)) == 0

    ## check the marriage function  -- reset the population 
    th.init_population()
    howmany = th.create_partners()
    th.marry()

    a3 = findall(x -> x.married == true, humans)
    @test length(a3) > 0
    @test length(a3) <= howmany ## cant have more married people than partnered people
    
    # check if there is an even number of married pairs
    a4 = findall(x -> x.married == true, humans)
    @test mod(length(a4), 2) == 0

    # check if there are any young people getting married
    a5 = findall(x -> x.married == true && x.age < 19, humans)    
    @test length(a5) == 0

    ## check if someone is married, make sure they have a partner.
    a6 = findall(x -> x.married == true && x.partner == 0, humans)
    @test length(a6) == 0

    ## go through all humans. If human is married, check if their partner is also married
    a7 = findall(humans) do h
        h.married == true && humans[h.partner].married != true
    end
    @test length(a7) == 0

    ## check that after a re-shuffle, the married pairs remain the same.
    th.init_population()
    th.create_partners()
    th.marry()
    a = findall(x -> x.married == true, humans)
    ## do the shuffle
    th.create_partners()
    b = findall(x -> x.married == true, humans)
    @test a == b

    ## testing  get_partners() function
    ## check for uniqueness in the returned pairs from get_partners() 
    th.init_population()
    th.create_partners()
    th.marry()    
    
    pairs = th.get_partners()
    a = [p.a for p in pairs]
    b = [p.b for p in pairs]
    @test a ∉ b && b ∉ a
end

@testset "IDIS" begin
    ## this tests the initial disease distribution. 
    ## Basically the function init_disease
    th._resetdemo() ## reset the demographics (init_pop, partners, marry)
 
    ## we havn't run the init_disease function yet. 
    ## this test is order dependent
    sickpartners = th.get_partners(onlysick = true)
    @test length(sickpartners) == 0 

    ## check init_disease() function.
    cnt = th.init_disease()
    a1 = findall(x -> x.health == th.INF, humans)
    @test length(a1) == cnt

    ### see if we can find a pair that's not sick. 
    ## also check for uniqueness in the sick only returned pairs      
    pairs_sick = th.get_partners(onlysick = true)
    @test length(pairs_sick) > 0 ## we should have positive lenght now 
    a = [p.a for p in pairs_sick]
    b = [p.b for p in pairs_sick]
    @test a ∉ b && b ∉ a
    r = findall(pairs_sick) do p
        p.a == th.INF || p.b == th.INF
    end
    @test length(r) == 0

    ## test init_disease() probabilities. 
    ## test whether the initialized diseased is within the confidence limits of the data paper.
    ## TO DO


    ## here we are checking that the total number of infected is properply distributed in the partners
    ## this is tested using the pair_stats() function. 
    th._reset() ## reset the population. 
    ## prevalence counts with partners. these are total # of individuals with partners (some infected individuals may not be without partners)
    prevcnt = length(findall(x -> x.partner > 0 && x.health == th.INF, humans))
    allpartnercnt = length(th.get_partners())
    sickpartnercnt = length(th.get_partners(onlysick = true)) ## get the total number of pairs with either inf/susc or inf/inf
    ## run transmission. the last ctr_dis is obviously dependent on beta but we don't care about that right now. 
    _a, _b, _c = th.pair_stats()
    @test _a == allpartnercnt
    @test _b == sickpartnercnt
    # _c contains all the pairs with only INF/SUSC
    infinf = _b - _c
    @test _c + infinf*2 == prevcnt
end

@testset "SHED" begin
    th._reset() ## reset the model with infected people
    ## no one is vaccinated or treated ## check this just incase
    for x in humans
        @test x.vaccinated == false && x.treated == false && x.newlyinfected == false
    end
    
    ## this testset only tests _get_shedding_weeks() and builds a distribution 
    ## can change the size from `gridsize` to anything else if needed a higher resolution distribution
    ds = zeros(Float64, gridsize)
    ss = zeros(Float64, gridsize)
    da = zeros(Float64, gridsize)
    sa = zeros(Float64, gridsize)
    for i = 1:gridsize        
        ds[i], ss[i], da[i], sa[i]  = th._get_shedding_weeks(humans[1])
    end
    println("mean ds: $(mean(ds)), ss: $(mean(ss)), da: $(mean(da)), sa: $(mean(sa))")

    ## TO DO: Test after person is vaccinated/treated that symptomatic days/shedding days are reduced
    humans[1].treated = true
    ds = zeros(Float64, gridsize)
    ss = zeros(Float64, gridsize)
    da = zeros(Float64, gridsize)
    sa = zeros(Float64, gridsize)
    for i = 1:gridsize        
        ds[i], ss[i], da[i], sa[i]  = th._get_shedding_weeks(humans[1])
    end
    println("after treatment: mean ds: $(mean(ds)), ss: $(mean(ss)), da: $(mean(da)), sa: $(mean(sa))")
end

# @testset "TREAT" begin
#     ## test the treatment function
#     ## test that with and without treatment what the average number of symptomatic days there are. 
#     _reset() ## reset the population with infected.
#     cnt = length(findall(x -> x.health == th.INF, humans))
 
#     ## test with 100% episodic treatment. this is usually at the start of warmup period and for calibration.
#     treatment(2, 1.0) ## give everyone episodic treatment
#     for x in humans
#         if x.health == th.INF
#             @test x.treated == 2
#         end
#     end
 
#     ## test with 0% episodic treatment. this is usually at the start of warmup period and for calibration.
#     _reset() ## reset the population with infected.
#     cnt = length(findall(x -> x.health == th.INF, humans))
#     treatment(2, 0.0)
#     for x in humans
#         if x.health == th.INF
#             @test x.treated == 0
#         end 
#     end

#     ## test whether shedding is reduced as expected. 
#     _reset() ## reset the population with infected. 
#     nta = Array{Int64, 1}()
#     nts = Array{Int64, 1}()
#     for x in humans
#         if x.health == th.INF
#             a, s = th._get_shedding_weeks(x) ## get number of weeks shedding asymptomatic/symptomatic
#             push!(nta, a)
#             push!(nts, s)
#         end 
#     end

#     _reset() ## reset the population with infected. 
#     treatment(2, 1.0)  #give everyone episodic treatment
#     ta = Array{Int64, 1}()
#     ts = Array{Int64, 1}()
#     for x in humans
#         if x.health == th.INF
#             a, s = th._get_shedding_weeks(x) ## get number of weeks shedding asymptomatic/symptomatic
#             push!(ta, a)
#             push!(ts, s)
#         end 
#     end
    
#     ## test how different nts (no treatment symptomatic periods) and ts (treatment symptomatic periods) differ. 


#     _reset() ## reset the population with infected. 
#     treatment(1, 1.0)  #give everyone suppressive treatment
#     ta = Array{Int64, 1}()
#     ts = Array{Int64, 1}()
#     for x in humans
#         if x.health == th.INF
#             a, s = th._get_shedding_weeks(x) ## get number of weeks shedding asymptomatic/symptomatic
#             push!(ta, a)
#             push!(ts, s)
#         end 
#     end

# end
