## main run file for HSV project. 
## this file loads all the packages required for analysis (the packages required for the simulations 
## are added automatically when adding the packing using `add thvaccine`).

## this file also sets up the parallelism: the code you see is specific to ABM's high performance cluster.  
## you can change this to connect to another cluster or just run `addprocs()` 
## to connect to your local PC. See Julia's documentation on `addprocs`.

using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using DataFrames
using CSV
using ClusterManagers
using Query
using Statistics
using thvaccine
using UnicodePlots
using Dates
using DependentBootstrap


#addprocs(SlurmManager(10*32), N=32) ## this throws an error because N = 32. Create PR. 
addprocs(SlurmManager(512), N=16) 
@everywhere using thvaccine
@everywhere using ProgressMeter

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