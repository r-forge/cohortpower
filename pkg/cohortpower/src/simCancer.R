simCohort = function(parameters,  populationData, 
                      deathFile, cancerType,
                      cancerRatesByType
	 )	
 {

# parameters is a list with elements
# geno env geno:env size probGeno probEnv genoMiss envMiss oneEnvPerCommunity
# cancerMissRate falseCancerRate sdGroup, sdPerson, Enrollmentprobs, 


if(!all(c("geno", "env", "geno:env", "size", "probGeno", "probEnv", "genoMiss", 
			"envMiss", "oneEnvPerCommunity", "cancerMissRate", "falseCancerRate", 
			"sdGroup", "sdPerson", "Enrollmentprobs","Ncommunity") %in% 
		names(parameters)))
warning("some parameters missing")

#coefs = c(parameters[["geno"]], parameters[["env"]], parameters[["geno:env"]])

thecommunities = rownames(populationData$all)[sort(sample(dim(populationData$F)[1],
parameters$Ncommunity))]


# subset the communities
for(D in 1:length(populationData))
	populationData[[D]] = populationData[[D]][thecommunities,]
      
# find number of subjects per CSD
if (parameters$EqualCommunity==T) 
  thesample=rep(thecommunities, parameters$size/parameters$Ncommunity) 

if (parameters$EqualCommunity==F)  {

#{community = rownames(populationData$all)
population=populationData$all[,"total"]
thesample= sample(thecommunities, parameters$size, replace=T, prob=population/sum(population))
}

thepeople=data.frame(total=as.matrix(table(thesample)))

thepeople[["F"]] = rbinom(dim(thepeople)[1], thepeople$total, 
		populationData$all[rownames(thepeople), "propFemale"])
		
	thepeople$M = thepeople$total - thepeople$F 
                 
datamat = data.frame(CensusDivision = rep(rownames(thepeople), thepeople$total),
	Gender = rep( rep(c("F", "M"), length(thecommunities)), 
					c(t(thepeople[,c("F", "M")])))) 

#simulate age
theages = as.integer(colnames(populationData$F))

# assign an age group



ages = NULL
for(Dcommunity in thecommunities) {
	for(Dsex in c("F", "M")) {
		 	ages =  c(ages,  sample(theages, 
			 		thepeople[Dcommunity, Dsex],
		 		prob=populationData[[Dsex]][Dcommunity,],replace=T) )
		}
}




# ages are uniformly distributed within an age group
datamat$Age = ages + runif(parameters$size, 0, 5)

#simulate enrollment date
Enrollmentcdf=cumsum(parameters$Enrollment)/sum(parameters$Enrollment)
datamat$Enrollment=approx(c(0,Enrollmentcdf), 0:4, runif(parameters$size))$y

# assign gene and environment groups
datamat$geno = rbinom(parameters$size, 1,parameters$probGeno)

if(parameters$oneEnvPerCommunity) {
	uniqueCommunities = as.character(unique(datamat$CensusDivision))
	env = sample(uniqueCommunities, round(length(uniqueCommunities)*parameters$probEnv))
	datamat$env=0
#	return(list(datamat, env))
	datamat[datamat$CensusDivision %in% env, "env"] = 1
} else {
	datamat$env = rbinom(parameters$size, 1,parameters$probEnv) 	
}


# individual and group random effects
personRandomEffect = rnorm(parameters$size, 0, parameters$sdPerson)
groupRandomEffect = rnorm(length(thecommunities), 0, parameters$sdGroup)

datamat$rr = personRandomEffect + groupRandomEffect[datamat$CensusDivision] + 
      parameters$geno * datamat$geno +
      parameters$env * datamat$env + 
      parameters[["geno:env"]] * datamat$geno *datamat$env -
      # subtract off one half the variance to correct the mean of the log-normal
      0.5*(parameters$sdPerson^2 + parameters$sdGroup^2)
datamat$rr = exp(datamat$rr)


thefemales = datamat$Gender== "F"

dataFemales = datamat[thefemales,]
dataMales = datamat[!thefemales,]

# simulate non-cancer deaths 

dataFemales$nonCancerDeath = 
  getDeath(deathFile$F, 
    dataFemales$Age, dataSim=dataFemales)

dataMales$nonCancerDeath =
  getDeath(deathFile$M, 
    dataMales$Age, dataSim=dataMales)
    

# rates for cancers which are and aren't of interest

desiredRates = otherRates = list()

for(D in c("M","F")) {
	allrate = cancerRatesByType[[D]]
	otherRates[[D]] = data.frame(x=allrate$age,
		y=apply(allrate[,-grep(paste("^age$|^",cancerType,"$", sep=""), colnames(allrate), ignore.case=T)], 1, sum))
	desiredRates[[D]] = data.frame(x=allrate$age,
			y=allrate[,grep(paste("^",cancerType,"$", sep=""), colnames(allrate), ignore.case=T)] )
}

# simulate cancer events that we aren't interested in

OtherCancerIncidence = simCancer(rep(1, dim(dataMales)[1]), dataMales$Age,
    	dataMales$nonCancerDeath, otherRates$M)

dataMales$nonCancerDeath = pmin(dataMales$nonCancerDeath, OtherCancerIncidence)


OtherCancerIncidence = simCancer(rep(1, dim(dataFemales)[1]), dataFemales$Age,
    	dataFemales$nonCancerDeath, otherRates$F)

dataFemales$nonCancerDeath = pmin(dataFemales$nonCancerDeath, OtherCancerIncidence)

# simulate from cancer we're interested in

dataMales$CancerIncidence = simCancer(dataMales$rr, dataMales$Age,
   dataMales$nonCancerDeath, desiredRates$M)
   

dataFemales$CancerIncidence = simCancer(dataFemales$rr, dataFemales$Age,
   dataFemales$nonCancerDeath, desiredRates$F)

datamat = rbind(dataMales, dataFemales)

noCancer = datamat$CancerIncidence>199

# misclassify cancers	
missedCancer = rbinom(parameters$size, 1, parameters$cancerMissRate) & !noCancer

falseRateThisCancer = parameters$falseCancerRate * (sum(desiredRates$M$y) + sum(desiredRates$F$y)) /(sum(otherRates$M$y) + sum(otherRates$F$y))

falseCancer = rbinom(parameters$size, 1, falseRateThisCancer) & noCancer

datamat$CancerIncidence[missedCancer] = 200
# missed cancers have uniform age distributions.
datamat$CancerIncidence[falseCancer] = runif(sum(falseCancer),  
	datamat$Age[falseCancer], datamat$nonCancerDeath[falseCancer])


# missclassify env and geno
whichGenoMiss = as.logical(rbinom(parameters$size, 1, parameters$genoMiss))
datamat$geno[whichGenoMiss] = !datamat$geno[whichGenoMiss] 

if(parameters$oneEnvPerCommunity) {
	toswitch = uniqueCommunities[as.logical(rbinom(length(uniqueCommunities), 1, parameters$envMiss))]
    toswitch = datamat$CensusDivision %in% toswitch
	datamat[toswitch,"env"]=as.integer(!datamat[toswitch, "env"])
} else {
	whichEnvMiss = as.logical(rbinom(parameters$size, 1, parameters$envMiss) )
	datamat$env[whichEnvMiss] = !datamat$env[whichEnvMiss] 
}


# compute observed events
datamat$CancerIncidence[is.na(datamat$CancerIncidence)] = 200
datamat$Cancer = datamat$nonCancerDeath > datamat$CancerIncidence
datamat$Event = pmin(datamat$nonCancerDeath, datamat$CancerIncidence)


datamat
				
									
}

simCancer = function(Ri, start, death, lambda) {
# simulate cancer incidence
# lambda is a list with x being the age sequence and y being the rates
# the other arguments are vectors, one element per individual
Nindiv =  length(Ri)
result = rep(0, Nindiv)

      
# add age zero with zero population, 
# and age 200 with zero population
#  in case an individual has 
# age below the lowest cutoff                          
ageSeq = c(0, lambda$x, 200)
lambdaSeq = c(0,lambda$y)                                
Nlambda = length(lambdaSeq)
diffLambda = diff(ageSeq)

# these will contain the cumulative distribution for each individual
thisAgeSeq = thisLambda = rep(0, Nlambda)

if(!is.loaded("simCancer")) 
  dyn.load(paste("../src/simCancer", .Platform$dynlib.ext, sep=""))
result = .C("simCancer" , as.double(Ri), as.double(start), 
  as.double(death), result=as.double(result), 
  as.integer(Nindiv), as.double(ageSeq), as.double(lambdaSeq),
  as.double(diffLambda), as.double(thisAgeSeq), as.double(thisLambda),
  as.integer(Nlambda))$result
  

  result
}

getDeath = function(mortalityFile, age, dataSim) {
# simulate non-cancer death

mydata = mortalityFile        

library(survival)

mort.out<-survreg(Surv(time2,aalive)~1,data=mydata,weights=weights,dist=c("weibull"))
shape<-mort.out$coefficients
scale<-mort.out$scale


# simulate non-cancer death
nonCancerDeath = rweibull(dim(dataSim)[1], scale=exp(shape), shape=1/scale)
smaller = nonCancerDeath < age
while(any(smaller)) {
 nonCancerDeath[smaller] = 
   rweibull(sum(smaller), scale=exp(shape), shape=1/scale)
 smaller = nonCancerDeath < age
}

return(nonCancerDeath )
}


####fuction of calculating different cancer lambda

lambdaDiffCancer=function(myfile)   {
require(foreign)
thefile = read.table(file=myfile, sep=",", header=T, as.is=T, quote="\"")

# get rid of rows with totals in them
thefile = thefile[grep("^[[:digit:]]", thefile[,2]),]


cancerRatesByType=unique(thefile[,1])
cancerRatesByType= cancerRatesByType[as.logical(nchar(cancerRatesByType))]

cancerRatesByType = cancerRatesByType[cancerRatesByType!="Total"]
cancerRatesByType = gsub(" [[:print:]]+$", "", cancerRatesByType)

Ncancer = length(cancerRatesByType)
Nages = length(unique(thefile[,2]))

cancerRatesByType = thefile[,1] = rep(cancerRatesByType, rep(Nages, Ncancer))

# reformat the dataframe, with nice column names and       
mydata=data.frame(cancerRatesByType = cancerRatesByType, age=as.integer(substr(thefile[,2],1,2)), 
	Incidence=as.integer(gsub(",","", thefile[,4])), Rate= as.numeric(thefile[,5]),
	Population = as.numeric(thefile[,6]) )

lambda = reshape(mydata[,c("cancerRatesByType","age","Rate")], direction="wide", timevar="cancerRatesByType", idvar="age")

# convert rates to per person rather than per 100,000.
lambda[,grep("^Rate", colnames(lambda), ignore.case=T)] = lambda[,grep( "^Rate", colnames(lambda),ignore.case=T)]/100000 

colnames(lambda) = gsub("^Rate\\.", "", colnames(lambda))



lambda
}

lambdaDiffCVD=function(myfile,agegroup)   {
require(foreign)
thefile = read.table(file=myfile, sep=",", header=T, as.is=T, quote="\"")

# get rid of rows with totals in them

cancerRatesByType=c(unique(thefile[,1]),"others")
Ncancer = length(cancerRatesByType)

Nages = length(unique(thefile[,2]))

# reformat the dataframe, with nice column names and       
mydata=data.frame(cancerRatesByType = thefile[,1] , age=as.integer(substr(thefile[,2],1,2)), 
	Incidence=as.integer(gsub(",","", thefile[,3])), Rate= as.numeric(thefile[,4]),
	Population = as.numeric(thefile[,5]) )


theages = unique(mydata$age)

others =as.numeric( mydata[mydata$cancerRatesByType=="Total",]$Incidence -  mydata[mydata$cancerRatesByType=="hf",]$Incidence-
mydata[mydata$cancerRatesByType=="ihd",]$Incidence-
mydata[mydata$cancerRatesByType=="stroke",]$Incidence )

otherspop = mydata[mydata$cancerRatesByType=="Total",]$Population
othersrate = others/otherspop*100000
others = cbind(cancerRatesByType=rep("others", length(theages)), age=as.numeric(theages), Incidence=others, Rate=othersrate, Population=otherspop)

mydata = rbind(mydata, others) 
lambda = matrix(as.numeric(mydata$Incidence), nrow=length(theages), ncol=Ncancer, 
	dimnames = list(as.character(theages), cancerRatesByType)) 

CancerTotal = lambda[,1]

lambda = lambda / matrix(CancerTotal, nrow=length(theages), ncol=Ncancer)

return(lambda)
}
