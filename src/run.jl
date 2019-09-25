## main run file for HSV project. 
using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using DataFrames
using CSV
using ClusterManagers

## multi core. 
##addprocs(4; exeflags="--project=.")
## a = SlurmManager(200, N=17)

#addprocs(SlurmManager(10*32), N=32) ## this throws an error because N = 32. Create PR. 
addprocs(SlurmManager(512), N=16) 

@everywhere using thvaccine
@everywhere using ProgressMeter


function calibration(numofsims=0, warmup_beta=0.0, main_beta=0.0, 
    warmup_time=0, eql_time=0, run_time=0, includetreatment=false)    
    ## this function runs a certain number of simulations for a given warmup beta/sim_beta.
    # calibration results (sept 19) 
    # P.sim_time = 100 or above. Use the last 10 years. 
    # cd = calibration(500, 0.02, 0.25, 50) 
    # Reulst are cyclic because people come in and out. 
    # changing beta to 1 dosn't make prev 100%. 
    # the bump of beta after warmup period puts it in some equilibrium. 
    # sims = RemoteChannel(()->Channel{Tuple}(numofsims));
    res = @showprogress pmap(1:numofsims) do x
        main(x, warmup_beta, main_beta, warmup_time, eql_time, run_time, includetreatment)
    end 
    avgprev =  zeros(Int64, warmup_time+eql_time+run_time, numofsims)   
    for i = 1:numofsims
        avgprev[:, i] = res[i].prevalence.Total
    end
    #hcat((ai.prevalence.Total for ai in a))...)
    #reduce(hcat, [ai.prevalence.Total for ai in a]) to avoid splatting
    #round.(mean(avgprev, dims=2), digits = 3)
    av = dropdims(round.(mean(avgprev, dims=2), digits = 2), dims=2)
    return av
end

function loopoverbetas() 
    ## this function loops over betas and calls the calibration function for each beta.
    ## if the average prevalence from start to end is less than 50, we consider that beta good.
    betas = round.([0.035 + 0.0001i for i in 0:25]; digits = 5)
    avgs = zeros(Float64, length(betas))
    @everywhere @eval thvaccine P.sim_time = 100
    println("Simulation Time: $(thvaccine.P.sim_time)")
    for i in 1:length(betas)
        println("starting calibration for: $(betas[i])")
        cd = calibration(100, betas[i])
        avgs[i] = cd[end]
        println("...average prevalence: $(cd[end])")
        println("...starting next sim \n")
        if abs(cd[end] - cd[1]) <= 50 
            println("found suitable beta: $(beta[i])")
            break
        end
    end
    return avgs
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