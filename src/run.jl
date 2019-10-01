## main run file for HSV project. 
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

## multi core. 
##addprocs(4; exeflags="--project=.")
## a = SlurmManager(200, N=17)

#addprocs(SlurmManager(10*32), N=32) ## this throws an error because N = 32. Create PR. 
addprocs(SlurmManager(512), N=16) 

@everywhere using thvaccine
@everywhere using ProgressMeter

@everywhere const numofsims = 500
@everywhere const warmup_beta=0.016
@everywhere const main_beta=0.08
@everywhere const warmup_time=40
@everywhere const eql_time=100 
@everywhere const run_time=10
@everywhere const totaltime = warmup_time+eql_time+run_time


## description of the model logic.
# model is at equilibrium with some beta and 12% prevalence. 
# this is with episodic treatment, automatically (see _get_shedding_weeks()). 

# at the start of the year, a few things will happen. 
# -> transmission is run. this looks at all infinf/infsusc. 
# -> for infinf, we simply have to run the natural history of disease. this is important because we would like to count the total number of symptom days (reduced) under different scenarios. 
# -> for infsusc, we run the natural history of disease for the infected. the symptomatic days/shedding days are recorded  
# -> -> in addition for infsusc, disease transfer can take place. if it happens, run the natural history of this newly infected person as well. 
# -> total number of symptomatic days/shedding days in year is recorded 
# -> total number of new infections is also recorded

# -> at the end of the year
# -> if a person is infected (i.e. newly infected or past infected) and suppressive treatment is on, 
# -> make this person x.treated = 1. This means they are under suppressive treatment. 
# -> their shedding is going to be significantly lowered. 
# -> how to use for cost-effectiveness? record their age. they will be under suppressive until they hit 49

# -> age function. 
# -> increase everyone's age by one. this is how the population is refreshed every year. 
# -> everyone over 49 leaves and is replaced by a 15 year old. 
# -> this is captured in the agedist dataframe

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