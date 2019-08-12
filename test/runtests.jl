using Test, thvaccine

const th = thvaccine
@testset "Demographics" begin
    ## init_human test
    x = thvaccine.Human()
    thvaccine.init_human(x, 1)
    @test x.age >= 15 && x.age <= 49
    @test x.partner == 0 && x.married == false
    @test x.firstyearinfection == true
    thvaccine.replace_human(x)  ## replace the human
    @test x.partner == 0
    @test x.age == 15

    thvaccine.init_population()
    a1 = [humans[i].age for i = 1:gridsize]
    @test length(findall(x -> x < 15 || x > 49, a1)) == 0
end

@testset "Sexual Partners" begin
    ## check if there are an even number of partners
    ## AND  check if all partners are mutual (a1)
    thvaccine.init_population()
    thvaccine.partnerup()
    allpartners = findall(x -> x.partner > 0, humans)
    a1 = findall(a -> humans[humans[a].partner].partner != a, allpartners)
    @test mod(length(allpartners), 2) == 0
    @test length(a1) == 0
    @test length(allpartners) > 0
    # check if they are all within their age groups
    a2 = map(allpartners) do h
        age1 = humans[h].age
        age2 = humans[humans[h].partner].age
        diffage = abs(age1 - age2) <= 4     
    end
    @test false ∉ a2

    # check if they are all within their ethnicites/subgroups
    a2 = findall(allpartners) do x
        grp1 = humans[x].grp
        grp2 = humans[humans[x].partner].grp
        grp1 != grp2
    end
    @test length(a2) == 0

    # no one should be married without partnering up
    thvaccine.init_population()
    thvaccine.marry()
    a3 = findall(x -> x.married == true, humans)
    @test length(a3) == 0

    # check if there are an even number of people married after partnering up
    thvaccine.init_population()
    thvaccine.partnerup()
    thvaccine.marry()
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
    thvaccine.init_population()
    thvaccine.partnerup()
    thvaccine.marry()
    a = findall(x -> x.married == true, humans)
    ## do the shuffle
    thvaccine.partnerup()
    b = findall(x -> x.married == true, humans)
    @test a == b

    # if a married person leaves the population, they should reset
    thvaccine.init_population()
    thvaccine.partnerup()
    thvaccine.marry()
    id = findfirst(x -> x.married == true && x.age == 49, humans) ## in the test, this could return 0 but low chance
    partner = humans[id].partner
    totalbefore = findall(x -> x.married == true, humans)
    if id != nothing
        thvaccine.exit(humans[id])
        @test humans[id].partner == 0
        @test humans[id].married == 0
        @test humans[id].age == 15
        @test humans[partner].partner == 0
        @test humans[partner].married == 0
        @test humans[partner].age == 15       
    else 
        println("could not test `exit_population` due to stochastic effects")
    end
    totalafter = findall(x -> x.married == true, humans)
    @test length(totalafter) == length(totalbefore)

end

@testset "Disease functions" begin
    ## these tests include the testing of get_partners function as well. 
    thvaccine.init_population()
    thvaccine.partnerup()

    ## test basic lengths
    allpartners = thvaccine.get_partners()
    @test length(allpartners) == length(findall(x -> x.partner > 0, humans))/2
    

    sickpartners = thvaccine.get_partners(onlysick = true)
    @test length(sickpartners) == 0 ## we havn't run the init_disease function yet. 

    thvaccine.init_disease()
    sickpartners = thvaccine.get_partners(onlysick = true)
    @test length(sickpartners) > 0 ## we should have positive lenght now
    @test length(sickpartners) <= length(allpartners) ## basic bug check

    ## see if we can find a pair that's not sick. 
    a = findall(sickpartners) do p 
        humans[p.a].health != INF && humans[p.b].health != INF 
    end
    @test length(a) == 0

    ## test init_disease() probabilities. 
    ## test whether the initialized diseased is within the confidence limits of the data paper.
   
end

@testset "Misc/Generic Functions" begin
    ## test the age groups.. should be from 1-4
    a = [thvaccine.get_age_group(humans[i].age) for i = 1:gridsize]
    @test length(findall(x -> x ∉ (1, 2, 3, 4), a)) == 0

    # ## test model parameters -- distribution of ethnicities
    # c = [P.eth_white, P.eth_black, P.eth_asian, P.eth_hispanic]
    # @test sum(c) == 1
end
