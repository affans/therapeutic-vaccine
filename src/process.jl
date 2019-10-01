## this file runs the main scenarios and processes the results. 
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

    ## create cumalative sums of year 141 to year 150
    dfs = []
    for i in 0:9
        #q(row) = 141 <= row.year <= 141+i
        #df_temp = filter(q, d)
        df_temp = d |> @filter(141 <= _.year <= 141+i) |> @groupby(_.sim) |>
            @map({sim=key(_), sum_ds=sum(_.ds), sum_ss=sum(_.ss), sum_ds_nt=sum(_.ds_nt), sum_ss_nt=sum(_.ss_nt), sum_da=sum(_.da), sum_sa=sum(_.sa)}) |> DataFrame
        push!(dfs, df_temp)
    end
   
    # process prevalence file
    p = vcat([cd[i].prevalence for i = 1:length(cd)]...)
    pp = p |> @groupby(_.year) |> 
         @map({year=key(_), avg_prevalence=mean(_.Total), avg_new_infections=mean(_.NewInfections)}) |> DataFrame
    pfs = []
    for i in 0:9
        df_temp = p |> @filter(141 <= _.year <= 141+i) |> @groupby(_.sim) |>                  
                  @map({sim=key(_), sum_prevalence=sum(_.Total), sum_new_infections=sum(_.NewInfections)}) |> DataFrame
        push!(pfs, df_temp)
    end 
    
    # process age distribution file
    a = vcat([cd[i].agedist for i = 1:length(cd)]...)
    aa = a |> @groupby(_.year) |> 
         @map({year=key(_), avg_left=mean(_.left), avg_left_infected=mean(_.left_ct), avg_left_treated=mean(_.left_treated)}) |> DataFrame
    afs = []
    for i in 0:9
        df_temp = a |> @filter(141 <= _.year <= 141+i) |> @groupby(_.sim) |>                  
                  @map({sim=key(_), sum_left=sum(_.left), sum_left_infected=sum(_.left_ct), sum_left_treated=sum(_.left_treated)}) |> DataFrame
        push!(afs, df_temp)
    end 
    
    # treatment
    t = vcat([cd[i].treatment for i = 1:length(cd)]...)
    # add the episodic cost
    t[!, :ecost] .= d |> @map(_.ds_nt*3) |> collect
    sk = Gamma(1.94, 1.42)
    # treatment cost 
    if scenario == 1
        t[!, :tcost] .= t |> @map(_.total_treated*365*3)
    else 
        t[!, :tcost] .= t |> @map(_.total_treated*50)
    end
    t[!, :totalcost] .= t.ecost + t.tcost
    
    tt = t |> @groupby(_.year) |> @map({year=key(_), avg_treated=mean(_.total_treated)}) |> DataFrame
    # add cost of treatment
    if scenario == 2.0
        tt[!, :treatment_cost] .= tt |> @map(_.avg_treated*50) |> collect 
        tt[!, :episodic_cost] .= dd |> @map(_.avg_ds_nt*3) |> collect
    end
    if scenario == 1.0
        tt[!, :treatment_cost] .= tt |> @map(_.avg_treated*365*3) |> collect 
        tt[!, :episodic_cost] .= dd |> @map(_.avg_ds_nt*3) |> collect
    end
    tt[!, :totalcost] .= tt.treatment_cost + tt.episodic_cost
   
    tfs = []
    for i in 0:9
        df_temp = t |> @filter(141 <= _.year <= 141+i) |> @groupby(_.sim) |>                  
                  @map({sim=key(_), sum_treated=mean(_.total_treated), sum_cost = sum(_.totalcost)}) |> DataFrame
        push!(tfs, df_temp)
    end 

    # check the original code using the new dataframe
    # ttt = t |> @groupby(_.year) |> @map({year=key(_), tcost=mean(_.tcost)}) |> DataFrame 

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
    mkpath("$dn/modeloutput")
    println("created directory: $dn")
    for cov in ( 0.5) #, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        fp = string(Int(cov*100))
        println("running coverage: $fp%")
        t1 = single(1.0, cov; write=true, dirname=dn, showprogress=false);
        v1 = single(2.0, cov; write=true, dirname=dn, showprogress=false);  
        
        CSV.write("$dn/modeloutput/($fp)_vac_raw_disease.dat", v1.rd) 
        CSV.write("$dn/modeloutput/($fp)_sup_raw_disease.dat", t1.rd) 
        CSV.write("$dn/modeloutput/($fp)_vac_yearavg_disease.dat", v1.ya) 
        CSV.write("$dn/modeloutput/($fp)_sup_yearavg_disease.dat", t1.ya)
        
        idf = process_yearly_averages(t1, v1)        
        CSV.write("$dn/m$fp.dat", idf)
        
        for i = 1:10 ## for the yearly simulation sums
            adf = process_sim_sums(i, t1, v1)
            CSV.write("$dn/q$(fp)_yr$i.dat", idf)
        end
    end
end

function process_yearly_averages(t, v)
    ## this function takes the results of two "single()" runs and puts together a 
    ## dataframe that combines the two results. 
    year = [i for i = 1:totaltime]
    
    cst_t = t.ya.treatment_cost 
    cst_v = v.ya.episodic_cost 
   
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