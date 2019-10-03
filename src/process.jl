## this file runs the main scenarios and processes the results. 

## constants that don't need to change a lot once fixed. 
@everywhere  numofsims = 500
@everywhere  warmup_beta=0.016
@everywhere  main_beta=0.07
@everywhere  warmup_time=50
@everywhere  eql_time=100 
@everywhere  run_time=10
@everywhere  totaltime = warmup_time+eql_time+run_time


function single(scenario = 1.0, cov = 1.0, vcpi = 50, efficacy = 0.8; showprogress=true) 
    ## this function runs the main scenario of the model,
    ## it processes the results by taking the average of simulations
    ## and returns a named tuple with the results. 
    @everywhere @eval thvaccine P.scenario = $scenario
    @everywhere @eval thvaccine P.treatment_coverage = $cov
    @everywhere @eval thvaccine P.vaccine_efficacy = $efficacy

    println("single simulation details: beta: $(main_beta), total time = $totaltime")
    println("parameters: scenario = $scenario, coverage = $cov, efficacy = $efficacy, vcpi = $vcpi")
    if showprogress
        cd = @showprogress pmap(1:numofsims) do x
            main(x, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    else
        cd = pmap(1:numofsims) do x
            main(x, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    end

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

    ## first let's get the average prevalence -- this is done later on now. 
    # avgprev =  zeros(Int64, warmup_time+eql_time+run_time, numofsims)   
    # for i = 1:numofsims
    #     avgprev[:, i] = cd[i].prevalence.Total
    # end
    # av = dropdims(round.(mean(avgprev, dims=2), digits = 2), dims=2)
    # ways of taking the average row-wise
    # df |> @mutate(d=mean(_)) |> DataFrame
    # mean(eachcol(df))
    # mean.(eachrow(df))
    # map(mean, eachrow(df));
    
    ## process disease file
    d = vcat([cd[i].disease for i = 1:length(cd)]...)
    ## create yearly average
    dd = d |> @groupby(_.year) |>
        @map({year=key(_), avg_ds=mean(_.ds), avg_ss=mean(_.ss), avg_ds_nt=mean(_.ds_nt), avg_ss_nt=mean(_.ss_nt), avg_da=mean(_.da), avg_sa=mean(_.sa)}) |> DataFrame

    ## create cumalative sums of the first ten years. 
    stime = warmup_time+eql_time+1
    println("stime is: $stime \n")

    dfs = []
    for i in 0:9
        #q(row) = 141 <= row.year <= 141+i
        #df_temp = filter(q, d)
        df_temp = d |> @filter(stime <= _.year <= stime+i) |> @groupby(_.sim) |>
            @map({sim=key(_), sum_ds=sum(_.ds), sum_ss=sum(_.ss), sum_ds_nt=sum(_.ds_nt), sum_ss_nt=sum(_.ss_nt), sum_da=sum(_.da), sum_sa=sum(_.sa)}) |> DataFrame
        push!(dfs, df_temp)
    end
   
    # process prevalence file
    p = vcat([cd[i].prevalence for i = 1:length(cd)]...)
    pp = p |> @groupby(_.year) |> 
         @map({year=key(_), avg_prevalence=mean(_.Total), avg_new_infections=mean(_.NewInfections)}) |> DataFrame
    pfs = []
    for i in 0:9
        df_temp = p |> @filter(stime <= _.year <= stime+i) |> @groupby(_.sim) |>                  
                  @map({sim=key(_), sum_prevalence=sum(_.Total), sum_new_infections=sum(_.NewInfections)}) |> DataFrame
        push!(pfs, df_temp)
    end 
    
    # process age distribution file
    a = vcat([cd[i].agedist for i = 1:length(cd)]...)
    aa = a |> @groupby(_.year) |> 
         @map({year=key(_), avg_left=mean(_.left), avg_left_infected=mean(_.left_ct), avg_left_treated=mean(_.left_treated)}) |> DataFrame
    afs = []
    for i in 0:9
        df_temp = a |> @filter(stime <= _.year <= stime+i) |> @groupby(_.sim) |>                  
                  @map({sim=key(_), sum_left=sum(_.left), sum_left_infected=sum(_.left_ct), sum_left_treated=sum(_.left_treated)}) |> DataFrame
        push!(afs, df_temp)
    end 
    
    # treatment and cost analysis
    t = vcat([cd[i].treatment for i = 1:length(cd)]...)
    # add the episodic cost. this is the number of episodic days from the non treated individuals. 
    # this is because agents on suppressive treatment need not get episodic treatment. 
    sk = Gamma(1.94, 1.42)
    t[!, :ecost] .= d |> @map(_.ds_nt*rand(sk)) |> collect
    # treatment cost 
    if scenario == 1  ## treatment scenario
        t[!, :tcost] .= t |> @map(_.total_treated*365*rand(sk))
    else ## vaccine scenario
        t[!, :tcost] .= t |> @map(_.total_treated*vcpi)
    end
    t[!, :totalcost] .= t.ecost + t.tcost
    
    tt = t |> @groupby(_.year) |> @map({year=key(_), avg_treated=mean(_.total_treated), tcost=mean(_.tcost)}) |> DataFrame 

    tfs = []
    for i in 0:9
        df_temp = t |> @filter(stime <= _.year <= stime+i) |> @groupby(_.sim) |>                  
                  @map({sim=key(_), sum_treated=mean(_.total_treated), sum_cost = sum(_.totalcost)}) |> DataFrame
        push!(tfs, df_temp)
    end 

    ## merge all the dataframes together.
    rawdata = join(d, p, on=[:year, :sim])
    rawdata = join(rawdata, t, on=[:year, :sim])
    rawdata = join(rawdata, a, on=[:year, :sim])

    yearlyaverages = join(dd, pp, on=:year)
    yearlyaverages = join(yearlyaverages, aa, on=:year)
    yearlyaverages = join(yearlyaverages, tt, on=:year)

    simavgs = []
    for i in 1:10 
        tempdf = join(dfs[i], pfs[i], on=:sim)
        tempdf = join(tempdf, afs[i], on=:sim)
        tempdf = join(tempdf, tfs[i], on=:sim)
        push!(simavgs, tempdf)
    end
    return (rd=rawdata, ya=yearlyaverages, sa=simavgs)
end


function scenarios()
    ## 
    dn = "/data/hsvvaccine/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
   
    println("created directory: $dn")

    for eff in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
        for vcpi in (50, 75, 100, 150, 200, 250)
            for cov in (0.2, 0.4, 0.5, 0.6, 0.8) 
                fp = "eff$(Int(eff*100))_vcpi$(vcpi)_cov$(Int(cov*100))"
                mkpath("$dn/$fp/modeloutput")
                t1 = single(1.0, cov, vcpi, eff; showprogress=false);
                v1 = single(2.0, cov, vcpi, eff; showprogress=false);  

                CSV.write("$dn/$fp/modeloutput/vac_raw_disease.dat", v1.rd) 
                CSV.write("$dn/$fp/modeloutput/sup_raw_disease.dat", t1.rd) 
                CSV.write("$dn/$fp/modeloutput/vac_yearavg_disease.dat", v1.ya) 
                CSV.write("$dn/$fp/modeloutput/sup_yearavg_disease.dat", t1.ya)
                
                # process the individial scenarios into one dataframe
                idf = process_yearly_averages(t1, v1)        
                CSV.write("$dn/$fp/m$(Int(cov*100)).dat", idf)

                for i = 1:10 ## for the 10 years post warm-up simulation sums
                    adf = process_sim_sums(i, t1, v1)
                    CSV.write("$dn/$fp/q$(Int(cov*100))_yr$i.dat", adf)
                end
            end
        end
    end
end

function process_yearly_averages(t, v)
    ## this function takes the results of two "single()" runs and puts together a 
    ## dataframe that combines the two results. 
    year = [i for i = 1:totaltime]
    
    cst_t = t.ya.tcost 
    cst_v = v.ya.tcost 
   
    ds_t = t.ya.avg_ds
    ds_v = v.ya.avg_ds

    t_t = t.ya.avg_treated
    t_v = v.ya.avg_treated

    i_t = t.ya.avg_new_infections
    i_v = v.ya.avg_new_infections

    p_t = t.ya.avg_prevalence
    p_v = v.ya.avg_prevalence

    idf = DataFrame(yr=year, cost_supp = cst_t, cost_vacc = cst_v,
                             symp_days_supp = ds_t, symp_days_vacc = ds_v,
                             num_treated_supp = t_t, num_treated_vacc = t_v, 
                             new_infect_supp=i_t, new_infect_vacc=i_v,
                             prev_supp = p_t, prev_vacc = p_v)


    return idf
end


function process_sim_sums(idx, t, v)
    ## this function takes the results of two "single()" runs and puts together a 
    ## dataframe that combines the two results. 
    sim = [i for i = 1:500]
    
    cst_t = t.sa[idx].sum_cost
    cst_v = v.sa[idx].sum_cost

    ds_t = t.sa[idx].sum_ds
    ds_v = v.sa[idx].sum_ds

    t_t = t.sa[idx].sum_treated
    t_v = v.sa[idx].sum_treated

    i_t = t.sa[idx].sum_new_infections
    i_v = v.sa[idx].sum_new_infections

    p_t = t.sa[idx].sum_prevalence
    p_v = v.sa[idx].sum_prevalence
    
    idf = DataFrame(sim=sim, cost_supp = cst_t, cost_vacc = cst_v,                             
                             symp_days_supp = ds_t, symp_days_vacc = ds_v,
                             num_treated_supp = t_t, num_treated_vacc = t_v, 
                             new_infect_supp=i_t, new_infect_vacc=i_v,
                             prev_supp = p_t, prev_vacc = p_v)
    return idf
end