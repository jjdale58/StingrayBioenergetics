#run functions seperately for males and females then combine in simaverage
#sex: 1 = female, 2 = male, calculates mean weight at age from logistic function with sex specific parameters and pairs with age specific sd based on observed sd of weight at age data
#data = weights
WeightAtAge<-function(data,sex) {
	data<-data.frame(data)
	if(sex==1) {
		Wi<-63.7
		k<-0.23
		to<-14.61 } else
	{
		Wi<-35.6
		k<-0.25
		to<-11.65
	}
	mweight<-Wi/(1+exp(-k*(data$Age-to)))
	mweight<-cbind(data$Age,mweight,data$SD)
	colnames(mweight)<-c("Age","MeanWeight","SDWeight")
	return(mweight)
}
#Draws one size for each age class from a random distribution with a mean = logistic fit and sd = sd of observed weight at age data
popgen<-function(weights) { 
	mpop<-rnorm(nrow(weights),mean=weights[,2],sd=weights[,3])
	mpop<-cbind(weights[,1],mpop)
	return(mpop)
}
                     
#list returned from popgen goes into simulate
#run simulation model over one year
simulate<-function(population,temp,sex) { 
	ageArray<-array(dim=c(1,length(population[,1])))
	weightArray<-array(dim=c(365,length(population[,1])))
	dGRArray<-array(dim=c(365,ncol(weightArray)))
	dRMRArray<-array(dim=c(365,ncol(weightArray)))
	cfood<-array(dim=c(365,ncol(ageArray)))
		
	if(sex==1) {
		mu<-c(63700,0.23) #mean Wi (in grams) and k
		sigma<-matrix(c(1.9219,-0.00779,-0.00779,0.0000538),2,2) #variance-covariance matrix for Wi and k
	} else
	{
		mu<-c(35600,0.25)
		sigma<-matrix(c(0.85068,-0.00755,-0.00755,0.000102),2,2)
	}
	#create age labels for each column
	for(k in 1:length(population[,1])) ageArray[1,k]<-population[k,1]
	
	#generates a matrix of temperatures, one for each day of the year for each age class, from a random distribution bounded by the min and max for that day.
	#uses inverse-CDF method
	temp$lower<-pnorm(temp$min,mean=temp$Temp,sd=temp$SD) #probablility that temp<=min (Cumulative Distribution Function)
	temp$upper<-pnorm(temp$max,mean=temp$Temp,sd=temp$SD) #probability that temp>=max
	devProb<-t(mapply(runif,nrow(population),temp$lower,temp$upper)) #generate a random probability between the probability for min and max
	tempArray<-apply(devProb,2,qnorm,temp$Temp,temp$SD) #inverse of pnorm, given a probability, returns the number whose cumulative distribution matches the probability
	
	#generate 365 days of data for the following variables for each age class, varies by day
	getRtriangle<-function(mid,lower,upper,ncol){
		m<-matrix(nrow=365,ncol=ncol)
		for(i in 1:ncol) m[,i]<-rtriangle(365,mid,lower,upper)
		return(m)
	}
	SMRact<-getRtriangle(1.11,1.06,1.16,nrow(population)) #randomly choose activity multiplier from a triangular distribution, requires package VGAM
	WArray<-getRtriangle(0.27,0.24,0.30,nrow(population))
	SDArray<-getRtriangle(0.1,0.06,0.17,nrow(population))
	
	#randomly assign growth parameters to each individual, requires Package MASS, individual parameters
	grPArray<-mvrnorm(nrow(population),mu,sigma)
	WiArray<-t(grPArray[,1]) #create vector of Wi and k values. Since they are individual traits, no need to create matrix
	kArray<-t(grPArray[,2])
	
	#generate 365 days of MO2 coefficients for each age class, individual parameters
	SMRcovMatrix<-matrix(c(0.120093242,-0.003784565,-0.075403475,-0.0037845652,0.0004460741,0.0015116337,-0.075403475,0.001511634,0.049644247),3,3)
	SMRcoefMeans<-c(-2.3775254,0.77578034,1.48740566)
	randomSMR<-mvrnorm(nrow(population),SMRcoefMeans,SMRcovMatrix)
	SMRa<-matrix(nrow=365,ncol=nrow(population))
	SMRb<-matrix(nrow=365,ncol=nrow(population))
	SMRc<-matrix(nrow=365,ncol=nrow(population))
	SMRa<-t(randomSMR[,1])
	SMRb<-t(randomSMR[,2])
	SMRc<-t(randomSMR[,3])

	ed<-matrix(rnorm(nrow(population),mean=6030,sd=0.18),1,nrow(population))
	
	WeightArray<-matrix(nrow=365,ncol=nrow(population))
	dGRArray<-matrix(nrow=365,ncol=nrow(population))
	WeightArray[1,]<-t(population[,2]*1000)
	dGRArray[1,]<-(((kArray*WeightArray[1,]*(WiArray-WeightArray[1,]))/WiArray)/365) #currently in g/day, need to convert to J/day by multiplying (g/d) by the energetic content of ray (J/g) for consumption equation
	for(k in 2:365) {
		WeightArray[k,]<-WeightArray[(k-1),]+dGRArray[(k-1),]
		dGRArray[k,]<-(((kArray*WeightArray[k,]*(WiArray-WeightArray[k,]))/WiArray)/365)
	}
	
	m<-log10(WeightArray)
	t<-log10(tempArray)
	dRMRArray<-matrix(nrow=365,ncol=nrow(population))
	for(k in 1:365){
		dRMRArray[k,]<-((10^(SMRa+m[k,]*SMRb+t[k,]*SMRc))*SMRact[1,]*24)*13.59 #13.59 (oxycalorific coefficient (J/mg O2)) used to convert RMR to daily energetic consumption (J/day)
	}
	#calculates consumption in J/day
	for(k in 1:365) dGRArray[k,]<-dGRArray[k,]*ed #convert from g/day to J/day
	
	loss<-1-(WArray+SDArray)
	cfood<-(dGRArray+dRMRArray)/loss #formula for consumption from Dowd et al. 2006, dGRArray multiplied by J/g (energy content of ray) to convert growth rate into J/d, currently JSH value
	
	#calculate means for each age class for each parameter
	means<-array(dim=c(length(population[,1]),14))
	means[,1]<-colMeans(cfood)/1000
	means[,2]<-colMeans(dGRArray)/1000
	means[,3]<-colMeans(dRMRArray)/1000
	means[,4]<-colMeans(WeightArray)/1000
	means[,5]<-colMeans(WArray)
	means[,6]<-colMeans(SDArray)
	means[,7]<-colMeans(SMRa)
	means[,8]<-colMeans(SMRb)
	means[,9]<-colMeans(SMRc)
	means[,10]<-colMeans(SMRact)
	means[,11]<-colMeans(ed)/1000
	means[,12]<-colMeans(WiArray)
	means[,13]<-colMeans(kArray)
	means[,14]<-colMeans(tempArray)
	#puts age of each individual as column title for each array
	colnames(means)<-c("Daily Ration", "Growth","RMR","weight","Waste","SDA","SMRa","SMRb","SMRc","SMRact","Edensity","Wi","k","temp")
	
	return(means)
}
#weights = mean weight at age with sd, temp = daily average temp data, sex; 1=female, 2=male, rep = number of iterations
consumption<-function(weights,temp,sex,rep) {
	mweight<-WeightAtAge(weights,sex)
	AgeMeans<-vector("list",15) #one list for each age class
	for(i in 1:length(AgeMeans)) AgeMeans[[i]]<-matrix(nrow=rep,ncol=14,data=0) #fill matrices with 0
	count<-1 #row counter,equal to rep
	for(k in 1:rep) {
		pop<-popgen(mweight)
		sim<-simulate(pop,temp,sex)
		for(i in 1:nrow(sim)) AgeMeans[[i]][count,]<-sim[i,] #put means for each variable into row corresponding to rep for each age class(e.g. rep 1 means in row 1)
		count<-count+1
	}
	for(i in 1:length(AgeMeans)) colnames(AgeMeans[[i]])<-colnames(sim) 
	print("Done with simulation")
	AgeMeans
}
#takes output from consumption as inputs, one for each sex.  Averages male and female data for all reps of simulate
simAverage<-function(simF,simM) {
	means<-vector("list",length(simF))
	#average means for age classes 1 to 9, all females after age 9
	for(i in 1:10){
		x<-list(a=simF[[i]],b=simM[[i]])
		means[[i]]<-Reduce("+",x)/length(x)
	}
	#fill means for age classes 10 to 14 with female means
	for(i in 11:15) means[[i]]<-simF[[i]]
	means
}

#colnum = column index of variable of interest
#calculates 2.5, 50 and 97.5% quantiles for a specified variable for each age class, returns matrix
extractQuantiles<-function(simAll,colnum){
	quantiles<-matrix(nrow=15,ncol=3)
	for(i in 1:length(simAll)) quantiles[i,]<-t(as.matrix(quantile(simAll[[i]][,colnum],c(0.025,0.5,0.975))))
	name.start<-colnames(simAll[[1]])[colnum]
	colnames(quantiles)<-c(paste(name.start,"0.025",sep="_"),paste(name.start,"0.5",sep="_"),paste(name.start,"0.975",sep="_"))
	quantiles
}