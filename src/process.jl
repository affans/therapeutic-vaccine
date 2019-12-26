## this file runs the main scenarios and processes the results. 
## 

function single(scenario = 1.0, cov = 1.0, efficacy = 0.8, vcpi = 50, 
    warmup_beta=0.016, main_beta=0.07, warmup_time=50, eql_time=100, run_time=20; showprogress=true) 
    ## returns: 
    ## `rd` the raw data: this is all the data collected over sims and years. It will have sims x years rows. 
    ## `ya` the yearly averages: the averages of ALL years taken over all simulations. 
    ## `sa` a vector of 10 DataFrames, each element is the raw data of 10 years post warm-up time for 500 simulations 
    numofsims = 500
    totaltime = warmup_time+eql_time+run_time
    println("single simulation details: scenario = $scenario, sims: $numofsims, coverage = $cov, efficacy = $efficacy, vcpi = $vcpi")
    println("... warmup_beta: $(warmup_beta), main_beta: $(main_beta), total time = $totaltime")

    if showprogress
        cd = @showprogress pmap(1:numofsims) do x
            main(x, scenario, cov, efficacy, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    else
        cd = pmap(1:numofsims) do x
            main(x, scenario, cov, efficacy, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    end

    ## add some extra columns for processing, ideally this should be done at the simulation level but it dosn't matter
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
    
    ## first let's get the average prevalence, this is done later on in the yearly avgs as well, 
    # avgprev =  zeros(Int64, warmup_time+eql_time+run_time, numofsims)   
    # for i = 1:numofsims
    #     avgprev[:, i] = cd[i].prevalence.Total
    # end
    #av = dropdims(round.(mean(avgprev, dims=2), digits = 2), dims=2)
    # ways of taking the average row-wise
    # df |> @mutate(d=mean(_)) |> DataFrame
    # mean(eachcol(df))
    # mean.(eachrow(df))
    # map(mean, eachrow(df));
    
    ## vcat all the simulation data together (then merge in big dataframe)
    d = vcat([cd[i].disease for i = 1:length(cd)]...)
    p = vcat([cd[i].prevalence for i = 1:length(cd)]...)
    a = vcat([cd[i].agedist for i = 1:length(cd)]...)
    t = vcat([cd[i].treatment for i = 1:length(cd)]...)
    # add the episodic cost. this is the number of episodic days from the non treated individuals. 
    # this is because individuals on suppressive treatment need not get episodic treatment. 
    sk = Gamma(1.94, 1.42)
    t[!, :ecost] .= d |> @map(_.ds_nt*rand(sk)) |> collect
    # treatment cost 
    if scenario == 1  ## treatment scenario
        t[!, :tcost] .= t |> @map(_.total_treated*365*rand(sk))
    else ## vaccine scenario
        t[!, :tcost] .= t |> @map(_.total_treated*vcpi)
    end
    t[!, :totalcost] .= t.ecost + t.tcost

    ## merge all the data frames together
    rd = join(d, p, on=[:year, :sim])
    rd = join(rd, t, on=[:year, :sim])
    rd = join(rd, a, on=[:year, :sim])

    ## take the raw data and compute yearly averages
    ya  = rd |> @groupby(_.year) |> 
                @map({
                year = key(_), 
                avg_prevalence=mean(_.Total), 
                avg_new_infections=mean(_.NewInfections), 
                avg_ds=mean(_.ds), avg_ss=mean(_.ss), 
                avg_ds_nt=mean(_.ds_nt), avg_ss_nt=mean(_.ss_nt), 
                avg_da=mean(_.da), avg_sa=mean(_.sa), 
                avg_left_system=mean(_.left), 
                avg_left_system_infected=mean(_.left_ct), 
                avg_left_system_treated=mean(_.left_treated),
                avg_treated=mean(_.total_treated), 
                tcost=mean(_.totalcost) }) |> DataFrame

     ## take the raw data and sum up everything at the simulation level PER YEAR
    stime = warmup_time+eql_time+1
    yrsleft = totaltime - stime
    simavgs = []
    for i in 0:yrsleft
        df_temp = rd |> @filter(stime <= _.year <= stime+i) |> @groupby(_.sim) |>
        @map({sim=key(_), 
        sum_ds=sum(_.ds), sum_ss=sum(_.ss), 
        sum_ds_nt=sum(_.ds_nt), sum_ss_nt=sum(_.ss_nt), 
        sum_da=sum(_.da), sum_sa=sum(_.sa), 
        sum_prevalence=sum(_.Total), 
        sum_new_infections=sum(_.NewInfections), 
        sum_left=sum(_.left), 
        sum_left_infected=sum(_.left_ct), 
        sum_left_treated=sum(_.left_treated),
        sum_treated=sum(_.total_treated), 
        sum_cost = sum(_.totalcost)}) |> DataFrame
        push!(simavgs, df_temp)
    end 
   
    return (rd=rd, ya=ya, sa=simavgs)
end

function scenarios()
    ## this is a helper function that runs all the scenarios we want. 
    ## it calls the single() function for each scenario. 
    ## it calls `process` functions to extract relevant information from ya/sa retruend from single 
    ## and put it in its own data files. 
    savepathprefix = "/data/hsvvaccine/"
    dn = "$savepathprefix/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
    mkpath("$dn")
    println("saving results to folder: $dn")

    rt = 20 ## number of years post processing to save the q files.

    ## run baseline scenario
    baseline = single(1.0, 0.0, 0.0, 0.0; showprogress=false)
    CSV.write("$dn/baseline_raw_disease.dat", baseline.rd)     
    CSV.write("$dn/baseline_yearavg_disease.dat", baseline.ya) 
    
    ctr = 1
    for eff in (0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
        for vcpi in (50, 100, 150, 200, 250, 300)
            for cov in (0.2, 0.4, 0.5, 0.6, 0.8) 
                println("working scenario number: $ctr/180")
                ctr += 1                
                fp = "eff$(Int(eff*100))_vcpi$(vcpi)_cov$(Int(cov*100))"
                mkpath("$dn/$fp/modeloutput")     
                t1 = single(1.0, cov, eff, vcpi; showprogress=false);
                v1 = single(2.0, cov, eff, vcpi; showprogress=false);  

                CSV.write("$dn/$fp/modeloutput/vac_raw_disease.dat", v1.rd) 
                CSV.write("$dn/$fp/modeloutput/sup_raw_disease.dat", t1.rd) 
                CSV.write("$dn/$fp/modeloutput/vac_yearavg_disease.dat", v1.ya) 
                CSV.write("$dn/$fp/modeloutput/sup_yearavg_disease.dat", t1.ya)
                
                # process the individial scenarios into one dataframe
                idf = process_yearly_averages(t1, v1, baseline)        
                CSV.write("$dn/$fp/m$(Int(cov*100)).dat", idf)

                for i = 1:rt ## for the 10 or 20 years post warm-up simulation sums
                    adf = process_sim_sums(i, t1, v1, baseline)
                    CSV.write("$dn/$fp/q$(Int(cov*100))_yr$i.dat", adf)
                end
            end
        end
    end
end

function process_yearly_averages(t, v, b)
    ## this function takes the results of two "single()" runs and puts together a 
    ## dataframe that combines the two results. 
    year = [i for i = 1:maximum(b.ya.year)]
    
    cst_b = b.ya.tcost
    cst_t = t.ya.tcost 
    cst_v = v.ya.tcost 
   
    ds_b = b.ya.avg_ds
    ds_t = t.ya.avg_ds
    ds_v = v.ya.avg_ds

    t_b = b.ya.avg_treated
    t_t = t.ya.avg_treated
    t_v = v.ya.avg_treated

    i_b = b.ya.avg_new_infections
    i_t = t.ya.avg_new_infections
    i_v = v.ya.avg_new_infections

    p_b = b.ya.avg_prevalence
    p_t = t.ya.avg_prevalence
    p_v = v.ya.avg_prevalence

    idf = DataFrame(yr=year, cost_supp = cst_t, cost_vacc = cst_v, cost_base = cst_b,
                             symp_days_supp = ds_t, symp_days_vacc = ds_v, symp_days_base = ds_b,
                             num_treated_supp = t_t, num_treated_vacc = t_v, num_treated_based = t_b, 
                             new_infect_supp = i_t, new_infect_vacc = i_v, new_infect_base = i_b,
                             prev_supp = p_t, prev_vacc = p_v, prev_base = p_b)


    return idf
end


function process_sim_sums(idx, t, v, b)
    ## this function takes the results of two "single()" runs and puts together a 
    ## dataframe that combines the two results. 
    sim = [i for i = 1:500]
    
    cst_b = b.sa[idx].sum_cost
    cst_t = t.sa[idx].sum_cost
    cst_v = v.sa[idx].sum_cost

    ds_b = b.sa[idx].sum_ds
    ds_t = t.sa[idx].sum_ds
    ds_v = v.sa[idx].sum_ds

    t_b = b.sa[idx].sum_treated
    t_t = t.sa[idx].sum_treated
    t_v = v.sa[idx].sum_treated

    i_b = b.sa[idx].sum_new_infections
    i_t = t.sa[idx].sum_new_infections
    i_v = v.sa[idx].sum_new_infections

    p_b = b.sa[idx].sum_prevalence
    p_t = t.sa[idx].sum_prevalence
    p_v = v.sa[idx].sum_prevalence
    
    idf = DataFrame(sim=sim, cost_supp = cst_t, cost_vacc = cst_v, cost_base = cst_b,                             
                             symp_days_supp = ds_t, symp_days_vacc = ds_v, symp_days_base = ds_b,
                             num_treated_supp = t_t, num_treated_vacc = t_v, num_treated_base = t_b, 
                             new_infect_supp=i_t, new_infect_vacc=i_v, new_infect_base=i_b,
                             prev_supp = p_t, prev_vacc = p_v, prev_base = p_b)
    return idf
end

## recursive file moving in linux 
## https://stackoverflow.com/questions/8798153/recursively-move-files-of-certain-type-and-keep-their-directory-structure
