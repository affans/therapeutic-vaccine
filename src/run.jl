## main run file for HSV project. 
using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using JSON 

addprocs(4; exeflags="--project")
@everywhere using thvaccine


res = pmap(x -> main(x), 1:5)