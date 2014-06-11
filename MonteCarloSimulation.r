#Monte Carlo Simulation for brown stingray bioenergetics model
#published in:
#Dale, J.J., Drazen, J.C., Holland, K.N. (2013) Stingray life history trade-offs associated with nursery habitat
#use inferred from a bioenergetics model. Marine Biology 160(12): 3181 - 3192.

rm(list=ls())

library(VGAM)
library(MASS)

source('BioenergeticsFunctions v2.R')
load('BioenergeticsData.Rdata') #daily average temperature data, mean and sd of weight at age by sex (temp, weightsF, weightsM)

#2000 simulations for males and females
simRunsF<-consumption(weightsF, temp, 1, 2000)
simRunsM<-consumption(weightsM, temp, 2, 2000)
#average males and females
simAvg<-simAverage(simRunsF, simRunsM)
#extract quantiles for all variables
quantiles<-vector('list', dim(simAvg[[1]])[2])
for(i in 1: dim(simAvg[[1]])[2]){
	quantiles[[i]]<-extractQuantiles(simAvg, i)
}



