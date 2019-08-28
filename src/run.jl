## main run file for HSV project. 
using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using JSON 


## multi core. 
addprocs(4; exeflags="--project")
@everywhere using thvaccine

using DataFrames
using CSV


res = map(x -> main(x), 1:5)

for i = 1:5
    dts = res[i]
    insertcols!(dts.prevalence, 1, :sim => i)
    insertcols!(dts.partners, 1, :sim => i)
    insertcols!(dts.episodes, 1, :sim => i)
end

# prevalence
p = vcat([res[i].prevalence for i = 1:5]...)
e = vcat([res[i].episodes for i = 1:5]...)
s = vcat([res[i].partners for i = 1:5]...)

