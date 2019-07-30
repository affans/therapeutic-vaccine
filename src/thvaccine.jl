module thvaccine
    using Distributions
    using Parameters
    using Random
    include("./parameters.jl")
    include("./functions.jl")
    export P, humans, main, modelinfo
  
    const P = ModelParameters()
    const humans = Array{Human}(undef, P.num_of_humans)
    #const runsummary = Dict{String, Int64}()
    function modelinfo()        
        ans = length(findall(x -> x.sex == MALE, humans))
        println("Number of MALE: $ans")
        ans = length(findall(x -> x.sex == FEMALE, humans))
        println("Number of FEMALE: $ans")
        ans = length(findall(x -> x.partner > 0, humans))
        println("Number of people with partners: $ans (distinct pairs: $(div(ans, 2)))")
        ans = length(findall(x -> x.partner > 0 && x.married == true, humans))
        println("Number of people married: $ans (distinct pairs: $(div(ans, 2)))")

    end

    main() = main(1)
    function main(simnumber::Int64) 
        #Random.seed!(simnumber)        
        init_humans(); @info "Population initialized"          
        pair_everyone(); @info "Population paired up"
        marry(); @info "Pairs married"              
    end

    function 
    
end # module
