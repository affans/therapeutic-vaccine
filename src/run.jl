## main run file for HSV project. 
using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using DataFrames
using CSV

## multi core. 
addprocs(4; exeflags="--project=.")
@everywhere using thvaccine
@everywhere using ProgressMeter

function calibration(numofsims=0, beta=0.0)    
    sims = RemoteChannel(()->Channel{Tuple}(numofsims));
    res = @showprogress pmap(1:numofsims) do x
        main(x, false, beta, sims)
    end 

    avgprev =  zeros(Int64, thvaccine.P.sim_time, numofsims)   
    for i = 1:numofsims
        avgprev[:, i] = res[i].prevalence.Total
    end

    #round.(mean(avgprev, dims=2), digits = 3)
    av = dropdims(round.(mean(avgprev, dims=2), digits = 2), dims=2)

    return av
end
#@time res = map(x -> main(x), 1:2)

# for i = 1:5
#     dts = res[i]
#     insertcols!(dts.prevalence, 1, :sim => i)
#     insertcols!(dts.partners, 1, :sim => i)
#     insertcols!(dts.episodes, 1, :sim => i)
# end

# # prevalence
# p = vcat([res[i].prevalence for i = 1:5]...)
# e = vcat([res[i].episodes for i = 1:5]...)
# s = vcat([res[i].partners for i = 1:5]...)

function process_prev(res)
    avg_prev = DataFrame([Int64 for i = 1:5], [Symbol("sim$i") for i = 1:5], 20)
    #insertcols!(avg_prev, 6, :avg => 0)
    for i = 1:5
        avg_prev[!, Symbol("sim$i")] .= res[i].prevalence[:, :Total]
    end
    c = convert(Matrix, avg_prev[:, 1:5])
    m = dropdims(mean(c; dims=2), dims=2) 
    m = round.(m; digits=2)
    avg_prev[!, :avg] = m
    
    
    ## ways of taking the average row-wise
    # df |> @mutate(d=mean(_)) |> DataFrame
    # mean(eachcol(df))
    # mean.(eachrow(df))
    # map(mean, eachrow(df));
end


function __calibration(numofsims)
    println("running calibration with total sims = $numofsims")
    betas = round.([0.01 + 0.005i for i in 0:15]; digits = 3)
    dt = DataFrame([Float64, Float64], [:betas, :average])
    for b in betas
        println("Testing β = $b")
        res = @showprogress pmap(1:numofsims) do x
            main(x, false, b)
        end
        arr = zeros(Float64, numofsims)
        for i in 1:numofsims
            arr[i] = res[i].prevalence[20, :Total]
        end
        ap = mean(arr)
        println("average prevalence at 20 years = $ap")
        push!(dt, (b, ap))
    end  
    return dt
end




## this example shows that even though agents is defined at the global scope and is available all the time 
## the fact that each worker is "indepedent" and dosn't share memory, 
## and the fact that each worker is running one instance of main() at any time.. 
## means we are two processes can not share the array of agents. 
module tp
export agents, mymain
agents = Array{Int64, 1}(undef, 100)
function mymain(sim)
    oldvalue  = maximum(agents)
    sleep(0.1)
    fill!(agents, sim)
    newvalue = maximum(agents)
    println("welcome from $(myid()), old value of agents: $oldvalue, new value of agents: $newvalue")
    return maximum(agents)
end

end