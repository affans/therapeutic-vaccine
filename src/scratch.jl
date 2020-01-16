function tt() 
    lof = readdir("/data/hsvvaccine/0106_1150")
    filter!(lof) do x 
        spx = split(x, "_")
        if spx[1] == "q" 
            return true
        else 
            return false 
        end
    end
    return lof
end

function overview()
    qfiles = tt()
    suppvals = Array{Int64, 1}(undef, 0)
    vaccvals = Array{Int64, 1}(undef, 0)
    bothvals = Array{Int64, 1}(undef, 0)
    for x in qfiles
        spx = split(x, "_")
        sv = parse(Int, SubString(spx[2], 5))
        vv = parse(Int, SubString(spx[3], 5))
        cv = parse(Int, SubString(spx[4], 5))
        push!(suppvals, sv)
        push!(vaccvals, vv)
        push!(bothvals, cv)
    end
    println("Suppressive Scenarios: $(map(x -> "$x%", unique(suppvals)))")
    println("Vaccine Scenarios: $(map(x -> "$x%", unique(vaccvals)))")
    println("For each suppressive, comibation of vaccine: $(map(x -> "$x%", unique(bothvals)))")

    ## check if the supp/vacc columns are the same in all the "combination" files. 
     
    for sp in unique(suppvals)
        prefix = Regex("^(q_supp$sp)(_vacc$sp)(_\\w+)(_yr1)(.\\w+)+")  ## only really need to check for 'yr1'.. could also check for other years if needed
        println("finding all files with prefix $prefix")
        fq = filter(qfiles) do x
            m = match(prefix, x)
            m !== nothing && return true
        end
        println("filtered $(length(fq)) files, running sanity check")
        lastval = missing
        for x in fq
            df = CSV.File("/data/hsvvaccine/0106_1150/$x") 
            pval = sum(df.cost_supp) + sum(df.cost_vacc) + sum(df.symp_days_supp) + sum(df.symp_days_supp) + 
                sum(df.new_infect_supp) + sum(df.new_infect_vacc) + sum(df.prev_supp) + sum(df.prev_vacc)
            check = !ismissing(lastval) && pval ≠ lastval
            check && error("sanity check failed")
        end
    end
end

function scheck()
    qfiles = tt()
end