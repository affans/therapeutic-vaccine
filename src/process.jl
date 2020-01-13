## this file runs the main scenarios and processes the results. 
## you may use this file for reproducibility or running your own scenarios.
struct output
    ## all data frames are 500 rows (the number of bootstrap replicates) and 20 columns (the number of years)
    v_newinfects::DataFrame ## vaccine only, bootstrapped icer: cost/number of new infections
    v_sympdays::DataFrame   ## vaccine only, bootstrapped icer: cost/number of symp days
    s_newinfects::DataFrame ## suppressive only, bootstrapped icer: cost/number of new infections
    s_sympdays::DataFrame   ## suppressive only, bootstrapped icer: cost/number of symp days
    c_newinfects::DataFrame ## combination of suppressive/vaccine
    c_sympdays::DataFrame   ## combination of suppressive/vaccine
    function output(num_bootstrap)
        _names = Symbol.(["yr$i" for i = 1:20])
        df = DataFrame([Float64 for i = 1:20], _names, num_bootstrap)
        df .= 0     
        new(copy(df), copy(df), copy(df), copy(df), copy(df), copy(df))
    end
end

struct outputtwo
    ## all data frames are 500 rows (the number of bootstrap replicates) and 20 columns (the number of years)
    newinfects::DataFrame ##  bootstrapped icer: cost/number of new infections
    sympdays::DataFrame   ##  bootstrapped, bootstrapped icer: cost/number of symp days
    function outputtwo()
        _names = Symbol.(["yr$i" for i = 1:20])
        df = DataFrame([Float64 for i = 1:20], _names, 500)
        df .= 0     
        new(copy(df), copy(df))
    end
end

function single(scenario = :none, supp_cov=0.0, vacc_cov=0.0, efficacy = 0.8, vcpi = 50, 
    warmup_beta=0.016, main_beta=0.07, warmup_time=50, eql_time=100, run_time=20; showprogress=true) 
    ## returns: 
    ## `rd` the raw data: this is all the data collected over sims and years. It will have sims x years rows. 
    ## `ya` the yearly averages: the averages of ALL years taken over all simulations. 
    ## `sa` a vector of 10 DataFrames, each element is the raw data of 10 years post warm-up time for 500 simulations 
    numofsims = 500
    totaltime = warmup_time+eql_time+run_time
    #println("single simulation details: scen=$scenario, sc=$supp_cov, vc=$vacc_cov, eff=$efficacy, vcpi=$vcpi")
    #println("... warmup_beta: $(warmup_beta), main_beta: $(main_beta), total time = $totaltime")

    if showprogress
        cd = @showprogress pmap(1:numofsims) do x
            main(x, scenario, supp_cov, vacc_cov, efficacy, warmup_beta, main_beta, warmup_time, eql_time, run_time)
        end 
    else
        cd = pmap(1:numofsims) do x
            main(x, scenario, supp_cov, vacc_cov, efficacy, warmup_beta, main_beta, warmup_time, eql_time, run_time)
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
    t[!, :tcost] .= t |> @map(_.total_treated * 365 * rand(sk))
    t[!, :vcost] .= t |> @map(_.total_vaccinated * vcpi)
    t[!, :totalcost] .= t.ecost + t.tcost + t.vcost

    ## merge all the data frames together
    rd = join(d, p, on=[:year, :sim])
    rd = join(rd, t, on=[:year, :sim])
    rd = join(rd, a, on=[:year, :sim])

    ## take the raw data and compute yearly averages
    ## this takes the averages per year (and not at a cumulative level). 
    ## eg: "tcost" is the average of costs incurred in a particular year only. 
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
                avg_vaccinated=mean(_.total_vaccinated),  
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
        sum_vaccinated=sum(_.total_vaccinated),
        sum_cost = sum(_.totalcost)}) |> DataFrame
        insertcols!(df_temp, 1, :postyr => (i+1))
        push!(simavgs, df_temp)
    end 

    simavgvcat = vcat([simavgs[i] for i = 1:length(simavgs)]...)
   
    return (rd=rd, ya=ya, sa=simavgvcat)
end

function scenarios()
    ## this is a helper function that runs all the scenarios we want. 
    ## it calls the single() function for each scenario. 
    ## it calls `process` functions to extract relevant information from ya/sa retruend from single 
    ## and put it in its own data files. 
    savepathprefix = "/data/hsvvaccine"
    dn = "$savepathprefix/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
    mkpath("$dn")
    println("saving results to folder: $dn")


    ## run baseline scenario
    baseline = single(:none, 0.0, 0.0, 0.0; showprogress=false)
    #CSV.write("$dn/baseline_raw_disease.dat", baseline.rd)     
    #CSV.write("$dn/baseline_yearavg_disease.dat", baseline.ya) 
    
    alleffs = (0.4, 0.5, 0.6, 0.7, 0.8)
    allvcpi = (50, 100, 150, 200, 250)
    allsc = (0.2, 0.4, 0.6, 0.8)
    allvc = collect(0.10:0.10:1.0)
    tot_scenarios = length(alleffs) * length(allvcpi) * length(allsc) * length(allvc)
    pgs = Progress(tot_scenarios, 1)   # minimum update interval: 1 second

    for sc in allsc  ## loop over fixed coverage levels
        scstr = Int(round(sc*100)) 
        ## run suppressive only scenario. supp scenario dosn't depend on vc, eff, vcpi.
        t1 = single(:suppressive, sc, 0, 0, 0; showprogress=false);
        booted = bootstrap_icer(baseline.sa, t1.sa)              
        CSV.write("$dn/bt_supp_newinfects_cov$(scstr).dat", booted.newinfects)
        CSV.write("$dn/bt_supp_sympdays_cov$(scstr).dat", booted.sympdays)
        for eff in alleffs
            for vcpi in allvcpi                       
                efstr = Int(round(eff*100))               
                v1 = single(:vaccine, sc, sc, eff, vcpi; showprogress=false); 
                booted = bootstrap_icer(baseline.sa, v1.sa)
                CSV.write("$dn/bt_vacc_newinfects_eff$(efstr)_cov$(scstr)_vcpi$(vcpi).dat", booted.newinfects)
                CSV.write("$dn/bt_vacc_sympdays_eff$(efstr)_cov$(scstr)_vcpi$(vcpi).dat", booted.sympdays)
                for vc in allvc  ## the varying vaccine coverage with the suppressive coverage fixed.                    
                    vcstr = Int(round(vc*100)) 
                    b1 = single(:scenA, sc, vc, eff, vcpi; showprogress=false);  
                    booted = bootstrap_icer(baseline.sa, b1.sa)
                    CSV.write("$dn/bt_mix_newinfects_eff$(efstr)_cov$(scstr)_vcov$(vcstr)_vcpi$(vcpi).dat", booted.newinfects)
                    CSV.write("$dn/bt_mix_sympdays_eff$(efstr)_cov$(scstr)_vcoc$(vcstr)_vcpi$(vcpi).dat", booted.sympdays)
                    #idf = process_ya(t1.ya, v1.ya, b1.ya, baseline.ya) 
                    #adf = process_sa(t1.sa, v1.sa, b1.sa, baseline.sa) 
                    
                    # write the m file
                    #CSV.write("$dn/m_eff$(efstr)_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", idf)
                    #CSV.write("$dn/q_eff$(efstr)_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", adf)
                    # calculate the icer values by bootstrapping the data. 
                    #bt_icers = iceronefunction(adf) 

                    # # write the bootstraped files 
                    #CSV.write("$dn/b_eff$(efstr)_vacconly_newinfects_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", bt_icers.v_newinfects)
                    #CSV.write("$dn/b_eff$(efstr)_vacconly_sympdays_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", bt_icers.v_sympdays)
                    ##CSV.write("$dn/b_eff$(efstr)_supponly_newinfects_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", bt_icers.s_newinfects)
                    #CSV.write("$dn/b_eff$(efstr)_supponly_sympdays_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", bt_icers.s_sympdays)
                    #CSV.write("$dn/b_eff$(efstr)_combinat_newinfects_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", bt_icers.c_newinfects)
                    #CSV.write("$dn/b_eff$(efstr)_combinat_sympdays_sv$(scstr)_mv$(vcstr)_vcpi$(vcpi).dat", bt_icers.c_sympdays)

                    next!(pgs; showvalues = [(:eff, eff), (:vcpi,vcpi), (:sc, sc), (:vc, vc)]) #update the progress bar
                end 
            end
        end
    end
end


function merge_sa_columns(df1, df2; names=["base", "supp"])
    ## WIP: refactoring of process_yearly_averages
    newdf = DataFrame()
    colstoextract = [:postyr, :sim, :sum_cost, :sum_ds, :sum_treated, :sum_vaccinated, :sum_new_infections, :sum_prevalence]
    for colname in reverse(colstoextract)
        insertcols!(newdf, 1, Symbol("$(names[1])_$colname") => df1[:, colname])
        insertcols!(newdf, 1, Symbol("$(names[2])_$colname") => df2[:, colname])                
    end
    return newdf
end

function bootstrap_icer(base, intervention)
    ## calculates the icer based on IID bootstrapping
    ## input requires two DataFrames

    ## join the columns needed for cost-effectiveness analysis
    adf = merge_sa_columns(base, intervention, names=["base", "inter"])
    btreps = 500      ## number of bootstrap replicates
    out = outputtwo() ## setup the output dataframe

    for y in 1:20  ## for each of the postyears past the warm up phase         
        bt = adf |> @filter(_.base_postyr == y) |> DataFrame
        btdat = dbootdata(bt, numresample=btreps, bootmethod=:iid)

        mcost = [-mean(btdat[i].inter_sum_cost - btdat[i].base_sum_cost) for i = 1:btreps]
        mnewinfects = [mean(btdat[i].inter_sum_new_infections - btdat[i].base_sum_new_infections) for i = 1:btreps]
        msympdays = [mean(btdat[i].inter_sum_ds - btdat[i].base_sum_ds) for i = 1:btreps]
        out.newinfects[!, y] = mcost ./ mnewinfects
        out.sympdays[!, y] = mcost./ msympdays  
    end
    return out
end

# function process_ya(df_supp, df_vacc, df_comb, df_base)
#     ## WIP: refactoring of process_yearly_averages
#     newdf = DataFrame()
#     year = [i for i = 1:maximum(df_base.year)]
#     insertcols!(newdf, 1, :year => year)
#     colstoextract = [:tcost, :avg_ds, :avg_treated, :avg_vaccinated, :avg_new_infections, :avg_prevalence]
#     for colname in colstoextract
#         insertcols!(newdf, 2, Symbol("supp_$colname") => df_supp[:, colname])
#         insertcols!(newdf, 2, Symbol("vacc_$colname") => df_vacc[:, colname])
#         insertcols!(newdf, 2, Symbol("comb_$colname") => df_comb[:, colname])
#         insertcols!(newdf, 2, Symbol("base_$colname") => df_base[:, colname])        
#     end
#     #df = DataFrame(cst_b, cst_t, cst_v, cst_c)
#     return newdf
# end

# function process_sa(df_supp, df_vacc, df_comb, df_base)
#     ## WIP: refactoring of process_yearly_averages
#     newdf = DataFrame()
#     postyear = df_base.postyr
#     sims = df_base.sim
#     insertcols!(newdf, 1, :postyear => postyear)
#     insertcols!(newdf, 2, :sims => sims)

#     colstoextract = [:sum_cost, :sum_ds, :sum_treated, :sum_vaccinated, :sum_new_infections, :sum_prevalence]
#     for colname in colstoextract
#         insertcols!(newdf, 3, Symbol("supp_$colname") => df_supp[:, colname])
#         insertcols!(newdf, 3, Symbol("vacc_$colname") => df_vacc[:, colname])
#         insertcols!(newdf, 3, Symbol("comb_$colname") => df_comb[:, colname])
#         insertcols!(newdf, 3, Symbol("base_$colname") => df_base[:, colname])        
#     end
#     return newdf
# end

# ## Bootstrapping code 
# function iceronefunction(psa)
#     # adf = [CSV.read("$(qf)_yr$i") for i = 1:20]
#     bt_reps = 500 ## number of bootstrap replicates
#     out = output(bt_reps) ## initialize data structure to save
 
#     for y in unique(psa.postyear) 
#         adf = psa |> @filter(_.postyear == y) |> DataFrame
#         btdat = dbootdata(adf, numresample=bt_reps, bootmethod=:iid)

#         # calculate vaccine only ICER values from bootstrapped data
#         mcost = [-mean(btdat[i].vacc_sum_cost - btdat[i].base_sum_cost) for i = 1:bt_reps]
#         mnewinfects = [mean(btdat[i].vacc_sum_new_infections - btdat[i].base_sum_new_infections) for i = 1:bt_reps]
#         msympdays = [mean(btdat[i].vacc_sum_ds - btdat[i].base_sum_ds) for i = 1:bt_reps]
#         icer_vaccine_newinfect = mcost ./ mnewinfects
#         icer_vaccine_sympdays = mcost./ msympdays

#         # calculate suppressive ICER values from bootstrapped data
#         mcost = [-mean(btdat[i].supp_sum_cost - btdat[i].base_sum_cost) for i = 1:bt_reps]
#         mnewinfects = [mean(btdat[i].supp_sum_new_infections - btdat[i].base_sum_new_infections) for i = 1:bt_reps]
#         msympdays = [mean(btdat[i].supp_sum_ds - btdat[i].base_sum_ds) for i = 1:bt_reps]
#         icer_supp_newinfect = mcost ./ mnewinfects
#         icer_supp_sympdays = mcost./ msympdays

#         # calculate suppressive/vaccine combination ICER values from bootstrapped data
#         mcost = [-mean(btdat[i].comb_sum_cost - btdat[i].base_sum_cost) for i = 1:bt_reps]
#         mnewinfects = [mean(btdat[i].comb_sum_new_infections - btdat[i].base_sum_new_infections) for i = 1:bt_reps]
#         msympdays = [mean(btdat[i].comb_sum_ds - btdat[i].base_sum_ds) for i = 1:bt_reps]
#         icer_comb_newinfect = mcost ./ mnewinfects
#         icer_comb_sympdays = mcost./ msympdays

#         out.v_newinfects[:, y] = icer_vaccine_newinfect
#         out.v_sympdays[:, y] = icer_vaccine_sympdays
#         out.s_newinfects[:, y] = icer_supp_newinfect
#         out.s_sympdays[:, y] = icer_supp_sympdays
#         out.c_newinfects[:, y] = icer_comb_newinfect
#         out.c_sympdays[:, y] = icer_comb_sympdays
#     end
#     return out
# end

# # check if all columns are the same in two dataframes a and b
# function test_supp_scenarios() 
#     ## this tests that suppressive scenarios are fixed regarldess of vaccine coverage, efficacy, vcpi. 
#     ## this will eventually need to be moved to the tests. 
#     b = single(:suppressive, 40, 0, 0, 0; showprogress=false);
#     bnames = names(b.rd)
#     for eff in (0:0.1:0.5)
#         for vcpi in (50, 100, 150)
#             for vc in (0:0.10:0.5)
#                 t1 = single(:suppressive, 40, vc, eff, vcpi; showprogress=false);           
#                 for n in bnames
#                     if !(t1.rd[:, n] == b.rd[:, n]) 
#                         println("colname: $n")
#                     end
#                 end             
#                 println("next sim")   
#             end
#         end
#     end
# end


# function process_yearly_averages(t, v, c, b)
#     ## this function takes the results of two "single()" runs and puts together a 
#     ## dataframe that combines the two results. 
#     year = [i for i = 1:maximum(b.ya.year)]
    
#     cst_b = b.ya.tcost
#     cst_t = t.ya.tcost 
#     cst_v = v.ya.tcost 
#     cst_c = c.ya.tcost
   
#     ds_b = b.ya.avg_ds
#     ds_t = t.ya.avg_ds
#     ds_v = v.ya.avg_ds
#     ds_c = c.ya.avg_ds

#     t_b = b.ya.avg_treated
#     t_t = t.ya.avg_treated
#     t_v = v.ya.avg_treated
#     t_c = c.ya.avg_treated

#     v_b = b.ya.avg_vaccinated
#     v_t = t.ya.avg_vaccinated
#     v_v = v.ya.avg_vaccinated
#     v_c = c.ya.avg_vaccinated

#     i_b = b.ya.avg_new_infections
#     i_t = t.ya.avg_new_infections
#     i_v = v.ya.avg_new_infections
#     i_c = c.ya.avg_new_infections

#     p_b = b.ya.avg_prevalence
#     p_t = t.ya.avg_prevalence
#     p_v = v.ya.avg_prevalence
#     p_c = c.ya.avg_prevalence

#     idf = DataFrame(yr=year, cost_supp = cst_t, cost_vacc = cst_v, cost_both = cst_c, cost_base = cst_b,
#                              symp_days_supp = ds_t, symp_days_vacc = ds_v, symp_days_both=ds_c, symp_days_base = ds_b,
#                              num_treated_supp = t_t, num_treated_vacc = t_v, num_treated_both = t_c, num_treated_base = t_b,
#                              num_vaccinated_supp = v_t, num_vaccinated_vacc = v_v, num_vaccinated_both = v_c, num_vaccinated_base = v_b, 
#                              new_infect_supp = i_t, new_infect_vacc = i_v, new_infect_both = i_c, new_infect_base = i_b,
#                              prev_supp = p_t, prev_vacc = p_v, prev_both = p_c, prev_base = p_b)


#     return idf
# end

# function process_sim_sums(idx, t, v, c, b)
#     ## this function takes the results of two "single()" runs and puts together a 
#     ## dataframe that combines the two results. 
#     sim = [i for i = 1:500]
    
#     cst_b = b.sa[idx].sum_cost
#     cst_t = t.sa[idx].sum_cost
#     cst_v = v.sa[idx].sum_cost
#     cst_c = c.sa[idx].sum_cost

#     ds_b = b.sa[idx].sum_ds
#     ds_t = t.sa[idx].sum_ds
#     ds_v = v.sa[idx].sum_ds
#     ds_c = c.sa[idx].sum_ds

#     t_b = b.sa[idx].sum_treated     ## column will be zero since no one is treated in baseline
#     t_t = t.sa[idx].sum_treated  
#     t_v = v.sa[idx].sum_treated     ## column will be zero since vaccine scenario dosn't have treatment
#     t_c = c.sa[idx].sum_treated

#     v_b = b.sa[idx].sum_vaccinated  ## column will be zero since no one is vaccinated in baseline
#     v_t = t.sa[idx].sum_vaccinated  ## column will be zero since treatment scenario dosn't have vaccine
#     v_v = v.sa[idx].sum_vaccinated
#     v_c = c.sa[idx].sum_vaccinated

#     i_b = b.sa[idx].sum_new_infections
#     i_t = t.sa[idx].sum_new_infections
#     i_v = v.sa[idx].sum_new_infections
#     i_c = c.sa[idx].sum_new_infections

#     p_b = b.sa[idx].sum_prevalence
#     p_t = t.sa[idx].sum_prevalence
#     p_v = v.sa[idx].sum_prevalence
#     p_c = c.sa[idx].sum_prevalence
    
#     idf = DataFrame(sim=sim, cost_supp = cst_t, cost_vacc = cst_v, cost_both = cst_c, cost_base = cst_b,                             
#                              symp_days_supp = ds_t, symp_days_vacc = ds_v, symp_days_both = ds_c, symp_days_base = ds_b,
#                              num_treated_supp = t_t, num_treated_vacc = t_v, num_treated_both = t_c, num_treated_base = t_b,
#                              num_vaccinated_supp = v_t, num_vaccinated_vacc = v_v, num_vaccinated_both = v_c, num_vaccinated_base = v_b, 
#                              new_infect_supp=i_t, new_infect_vacc=i_v, new_infect_both=i_c, new_infect_base=i_b,
#                              prev_supp = p_t, prev_vacc = p_v, prev_both = p_c, prev_base = p_b)
#     return idf
# end

## recursive file moving in linux 
## https://stackoverflow.com/questions/8798153/recursively-move-files-of-certain-type-and-keep-their-directory-structure
