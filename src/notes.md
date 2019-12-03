vaccine effect:

vaccination effect on symptomatic periods:
capture the number of episodes
capture the duration of symp/shedding days/num of sex... everything. 
but make sure there is no disease transfer that happens. 

vaccination effect on asymptomatic period. 
if the efficacy is 100%, there is no shedding at all. 
if the efficacy is 50%, there is 50% less shedding (in days). 

with 100% efficacy, the person will not have any shedding.
with 50% efficacy, the person will shed 

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