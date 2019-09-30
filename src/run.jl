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
@everywhere const run_time=20
@everywhere const totaltime = warmup_time+eql_time+run_time

function single(scenario = 1.0, cov = 1.0; write=false, dirname=".", showprogress=true) 
    ## this function runs the main scenario of the model,
    ## it processes the results by taking the average of simulations
    ## and returns a named tuple with the results. 

    @everywhere @eval thvaccine P.scenario = $scenario
    @everywhere @eval thvaccine P.treatment_coverage = $cov

    if showprogress
        cd = @showprogress pmap(1:numofsims) do x
            main(x, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    else
        cd = pmap(1:numofsims) do x
            main(x, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    end


    ## filename setup
    if scenario == 1.0 
        fname = "supp_cov$(replace(string(cov), "." => ""))"
    else
        fname = "vacc_cov$(replace(string(cov), "." => ""))"
    end

    ## first let's get the average prevalence
    avgprev =  zeros(Int64, warmup_time+eql_time+run_time, numofsims)   
    for i = 1:numofsims
        avgprev[:, i] = cd[i].prevalence.Total
    end
    av = dropdims(round.(mean(avgprev, dims=2), digits = 2), dims=2)
    # ways of taking the average row-wise
    # df |> @mutate(d=mean(_)) |> DataFrame
    # mean(eachcol(df))
    # mean.(eachrow(df))
    # map(mean, eachrow(df));

    ## add some extra columns for processing]
    for i = 1:numofsims
        dts = cd[i]    
        insertcols!(dts.prevalence, 1, :sim => i)    
        insertcols!(dts.disease, 1, :sim => i)    
        insertcols!(dts.treatment, 1, :sim => i)
        insertcols!(dts.agedist, 1, :sim => i)    
        insertcols!(dts.partners, 1, :sim => i)
        insertcols!(dts.prevalence, 1, :year => 1:totaltime)
        insertcols!(dts.disease, 1, :year => 1:totaltime)
        insertcols!(dts.treatment, 1, :year => 1:totaltime)
        insertcols!(dts.agedist, 1, :year => 1:totaltime)
        insertcols!(dts.partners, 1, :year => 1:totaltime)
    end
    
    ## process disease file
    d = vcat([cd[i].disease for i = 1:length(cd)]...)
    dd = d |> @groupby(_.year) |> @map({year=key(_), avg_ds=mean(_.ds), avg_ss=mean(_.ss), 
                                                     avg_ds_nt=mean(_.ds_nt), avg_ss_nt=mean(_.ss_nt),
                                                     avg_da=mean(_.da), avg_sa=mean(_.sa)}) |> DataFrame
    round.(dd; digits=2)
   
    # prevalence
    p = vcat([cd[i].prevalence for i = 1:length(cd)]...)
    pp = p |> @groupby(_.year) |> @map({year=key(_), total_avg_prevalence=mean(_.Total), new_avg_infections=mean(_.NewInfections)}) |> DataFrame
    round.(pp; digits=2)
    
    # age dist
    a = vcat([cd[i].agedist for i = 1:length(cd)]...)
    aa = a |> @groupby(_.year) |> @map({year=key(_), avg_gr1=mean(_.gr1), avg_gr2=mean(_.gr2), avg_gr3=mean(_.gr3), avg_gr4=mean(_.gr4), avg_left=mean(_.left), avg_left_ct=mean(_.left_ct)}) |> DataFrame
    round.(aa; digits=2)
   
    # treatment
    t = vcat([cd[i].treatment for i = 1:length(cd)]...)
    tt = t |> @groupby(_.year) |> @map({year=key(_), avg_treated=mean(_.total_treated)}) |> DataFrame
    # add cost of treatment
    sk = Gamma(1.94, 1.42)
    if scenario == 2.0
        tt[!, :treatment_cost] .= tt |> @map(_.avg_treated*50) |> collect 
        tt[!, :episodic_cost] .= dd |> @map(_.avg_ds_nt*rand(sk)) |> collect
    end
    if scenario == 1.0
        tt[!, :treatment_cost] .= tt |> @map((totaltime + 1 - _.year)*365*_.avg_treated*rand(sk)) |> collect 
        tt[!, :episodic_cost] .= dd |> @map(_.avg_ds_nt*rand(sk)) |> collect
    end
    
    # partners
    x = vcat([cd[i].partners for i = 1:length(cd)]...)
    
    if write
        CSV.write("$dirname/modeloutput/raw_partners_$fname.dat", x)    
        CSV.write("$dirname/modeloutput/$(fname)_raw_disease.dat", d) 
        CSV.write("$dirname/modeloutput/$(fname)_avg_disease.dat", dd)         
        CSV.write("$dirname/modeloutput/$(fname)_raw_treatments.dat", t)    
        CSV.write("$dirname/modeloutput/$(fname)_avg_treatments.dat", tt)    
        CSV.write("$dirname/modeloutput/$(fname)_raw_agedist.dat", a)
        CSV.write("$dirname/modeloutput/$(fname)_avg_agedist.dat", aa)    
        CSV.write("$dirname/modeloutput/$(fname)_raw_prevalence.dat", p)
        CSV.write("$dirname/modeloutput/$(fname)_avg_prevalence.dat", pp)
        writedlm("$dirname/modeloutput/$(fname)_avg_prev.dat", av)      
    end

    return (av=av, dd=dd, pp=pp, tt=tt, aa=aa)
end

function scenarios()
    ## 
    dn = Dates.format(Dates.now(), dateformat"mmdd_HHMM")
    mkpath("$dn/modeloutput")
    println("created directory: $dn")
    for cov in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) #, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        fp = string(Int(cov*100))
        println("running coverage: $fp%")
        t1 = single(1.0, cov; write=true, dirname=dn, showprogress=false);
        v1 = single(2.0, cov; write=true, dirname=dn, showprogress=false);        
        process(t1, v1, "$dn/m$fp.dat")
    end
end

function process(t, v, fname)
    ## this function takes the results of two "single()" runs and puts together a 
    ## dataframe that combines the two results. 
    year = [i for i = 1:totaltime]
    
    ct_t = t.tt.treatment_cost 
    ce_t = t.tt.episodic_cost 
    
    ct_v = v.tt.treatment_cost 
    ce_v = v.tt.episodic_cost 

    ds_t = t.dd.avg_ds
    ds_v = v.dd.avg_ds

    t_t = t.tt.avg_treated
    t_v = v.tt.avg_treated

    i_t = t.pp.new_avg_infections
    i_v = v.pp.new_avg_infections

    p_t = t.pp.total_avg_prevalence
    p_v = v.pp.total_avg_prevalence

    idf = DataFrame(yr=year, cost_treatment_supp = ct_t, cost_episodic_supp = ce_t,
                             cost_treatment_vacc = ct_v, cost_episodic_vacc = ce_v,
                             symp_days_supp = ds_t, symp_days_vacc = ds_v,
                             num_treated_supp = t_t, num_treated_vacc = t_v, 
                             new_infect_supp=i_t, new_infect_vacc=i_v,
                             prev_supp = p_t, prev_vacc = p_v)

    CSV.write(fname, idf)
    return idf
end


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