using Test, thvaccine

@testset "Demographics" begin
    thvaccine.init_humans()
    a1 = [humans[i].age for i = 1:P.num_of_humans]
    @test length(findall(x -> x < 15 || x > 49, a1)) == 0
end

@testset "Sexual Partners" begin
    thvaccine.pair_everyone()
    allpartners = findall(x -> x.partner > 0, humans)
    @test mod(length(allpartners), 2) == 0 ## check if even number of partners
  
    # check if all partners are mutual
    a1 = findall(a -> humans[humans[a].partner].partner != a, allpartners)
    @test length(a1) == 0

    # check if they are all within their age groups
    a2 = findall(allpartners) do h
        age1 = thvaccine.get_age_group(humans[h].age)
        age2 = thvaccine.get_age_group(humans[humans[h].partner].age) 
        age1 != age2       
    end
    @test length(a2) == 0
    
    thvaccine.marry()
    a3 = findall(x -> x.married == true, humans)
    @test mod(length(a3), 2) == 0

    ## go through all humans. If human is married, check if their partner is also married
    a4 = findall(humans) do h
        h.married == true && humans[h.partner].married != true
    end
    @test length(a4) == 0

    ## check if someone is married, make sure they have a partner. 
    a = findall(x -> x.married == true && x.partner == 0, humans)
    @test length(a) == 0

    ## check that after a re-shuffle, the married pairs remain the same. 
    a = findall(x -> x.married == true, humans)
    thvaccine.pair_everyone()
    b = findall(x -> x.married == true, humans)
    @test a == b

    ## check that after a re-shuffle, the married pairs remain the same. 
    ## this is the second test.
    a = findall(x -> x.married == true, humans)
    thvaccine.pair_everyone()
    b = findall(x -> x.married == true, humans)
    @test a == b
end

@testset "Misc/Generic Functions" begin
    ## test the age groups.. should be from 1-4
    a = [thvaccine.get_age_group(humans[i].age) for i = 1:P.num_of_humans]
    @test length(findall(x -> x ∉ (1, 2, 3, 4), a)) == 0
end
