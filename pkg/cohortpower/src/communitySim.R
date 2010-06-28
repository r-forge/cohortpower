

#power function                               
# yesnorate means the cancer has been missed
# noyesrate means a false cancer diagnosis
Power=function(parameters, populationData,  CommonCancer,  
 	 SimulationTime, deathFile, cancerRatesByType, verbose=F, saveData=F) {
  # calculate the power, given parameters, rates

  # signifncance levels
  p = c(0.05, 0.01, 0.001, 0.0001)

	coefs = c("geno", "env", "geno:env")
	allData=list()
  
  # get rid of 'all' bit of common cancer because it's included automatically

  #CommonCancer = c(colnames(cancerRatesByType[[1]])[! CommonCancer %in% c("Total", "[[:alpha:]]")]

  
  if(is.null(CommonCancer)){
CommonCancer=	c(colnames(cancerRatesByType[[1]]) [ !colnames(cancerRatesByType[[1]])== "Total" ],"All")
}  
# else {
#CommonCancer = c(CommonCancer, "All")
#}
    
	parray=array(NA, c(length(parameters$Followup), length(coefs), length(CommonCancer),SimulationTime))
	dimnames(parray) = list(as.character(parameters$Followup), coefs, CommonCancer, NULL)

	carray=varray=parray
	retries=array(NA, c(length(parameters$Followup), length(CommonCancer),SimulationTime,3))
	dimnames(retries) = list(as.character(parameters$Followup),  CommonCancer, NULL,c("95pct", "cc4","cc1"))

	# loop through simulations, get all the p values
	for (i in 1:SimulationTime){
    res =Pvalue(parameters, populationData,  CommonCancer,deathFile, cancerRatesByType, verbose,
				saveData=saveData)
    parray[,,,i] = res
	carray[,,,i]= attributes(res)$coef
	varray[,,,i]= attributes(res)$var
#	return(list(parray, retries, res))
	retries[,,i,]=attributes(res)$retries
	if(saveData) {
		allData[[i]] = attributes(res)$data
	}
	
	cat(i)                                                   
    if(round(i/30)==i) cat("\n")
  }
  cat("\n")
                            
  power = array(NA, 
		c(length(parameters$Followup), length(coefs), length(CommonCancer),length(p)))
	for(Dsignif in 1:length(p)){
		 power[,,,Dsignif] = apply(parray, c(1,2,3), function(qq) mean(qq<p[Dsignif]))
	}

	dimnames(power) = list(as.character(parameters$Followup), coefs, CommonCancer,as.character(p))
	
	attributes(power)$var =varray
	attributes(power)$coef =carray
	attributes(power)$retries = retries
	attributes(power)$data = allData
	
	return(power)
}

######################################
# p value function starts from here
#####################################
Pvalue = function(parameters, populationData, 
	CommonCancer, deathFile, cancerRatesByType, verbose=F, saveData=F) { 
# if oneEnvPerCommunity=T, each community has only one environment effect.  
#  otherwise there's one effect per person.									

	if(length(CommonCancer)>1)
		warning("Only one cancer at a time is allowed, was given", CommonCancer)
	
	require(foreign)
	library(survival)
	

    controlProp=parameters$controlProp
    size=parameters$size
    coefs = unlist(parameters[c("geno", "env", "geno:env")])
    
# create array for p values
    
    estimates = standarderrors = pvalue = 
			array(1, c(length(parameters$Followup), length(coefs), length(CommonCancer)),
dimnames = list(as.character(parameters$Followup), 
		names(coefs), CommonCancer) )   


retries = 
		array(NA, c(length(parameters$Followup),  length(CommonCancer), 3),
				dimnames = list(as.character(parameters$Followup), 
						 CommonCancer, c("95pct", "cc8","cc1")) )   

# simulate a cohort


for(Dcancer in CommonCancer) {  
    
    
    if(verbose) cat(Dcancer)
    
    datamat = simCohort(parameters, populationData, deathData, Dcancer, cancerRatesByType) 
	
# if there's no one in the interaction groups, set p values to a half
	if (sum(datamat$geno&datamat$env) == 0) {
        pvalue[,,Dcancer]=0.5
	} else {
		
        for(Dfollowup in parameters$Followup) {
            if(verbose) cat(Dfollowup)
            # set thisEventTime equal to age at end of followoup period
            ageEndFollowup = datamat$Age + Dfollowup - datamat$Enrollment
            hadEvent = datamat$Event < ageEndFollowup
            datamat$thisEventType = hadEvent &  datamat$Cancer
            
            if(!any(datamat$thisEventType)) {
# if there is nobody with this event,  the pvalues to 1
                pvalue[as.character(Dfollowup),,Dcancer]= 1
                cat("n")
            } else {        
				
                # eventType=T if the cancer of interest occurred during the followup period
                
            	
                # censored failure time
                datamat$thisEventTime = datamat$Event
                datamat$thisEventTime[!hadEvent] = ageEndFollowup[!hadEvent]
				
				
    # for case-control study
  if ( !is.null(parameters$controlProp)   ) {
    controlCohort = datamat[datamat$thisEventType==F,]
    ncase = sum(datamat$thisEventType)
    controlgroup = controlCohort [sample(dim(controlCohort)[1], round(ncase*parameters$controlProp)),]
    # add the maximum as the whole cohort
    datamatFull = datamat                                                   
    datamat = rbind(datamat[datamat$thisEventType==T,], controlgroup)
  } 
# and add left truncation
  response = Surv(datamat$Age, datamat$thisEventTime, datamat$thisEventType)
	
		if(Dcancer == "Breast") {
		coxfitAfterYears=try(coxph(response~geno*env + frailty(CensusDivision),
  		data=datamat, subset=datamat$Gender=="F", 
				control=coxph.control(iter.max=30)))  
		} else if (Dcancer=="Prostate") {
		coxfitAfterYears=try(coxph(response~geno*env + frailty(CensusDivision),
  		data=datamat, subset=datamat$Gender=="M", 
				control=coxph.control(iter.max=30)) ) 		
		}  else { 
# add frailty


		coxfitAfterYears=try(coxph(response~geno*env+
			strata(Gender,na.group=TRUE) + frailty(CensusDivision), 
  		data=datamat, control=coxph.control(iter.max=30,toler.chol=.Machine$double.eps)) )

		# check if a parameter not estimated because of singular matrix
	 # if so, re-run with random 95% of the non-event data until it works
		Ntries <- 0
		if(class(coxfitAfterYears)[1]=="try-error")
				coxfitAfterYears = list(coxfitAfterYears, coef=NA) 
			theNoEvent = which(!datamat$thisEventType)
			theEvents = which(datamat$thisEventType)
			
			while(any(is.na(coxfitAfterYears$coef)) & Ntries <= 10) {
				NtoKeep = round(length(theNoEvent)*0.95)
				cat("s")
			coxfitAfterYears=try(coxph(response ~ geno*env +	strata(Gender,na.group=TRUE) + frailty(CensusDivision), 
							subset = sort(c(theEvents, sample(theNoEvent, NtoKeep))),
  							data=datamat, control=coxph.control(iter.max=30,toler.chol=.Machine$double.eps)) ) 
	if(class(coxfitAfterYears)[1]=="try-error")
		coxfitAfterYears = list(coxfitAfterYears, coef=NA) 
			Ntries = Ntries + 1
		}
		Ntries8 = 0
		# if still singular matrix, keep twice the number of cases as controls
	
	while(any(is.na(coxfitAfterYears$coef)) & Ntries8 <= 40) {
		NtoKeep = round(length(theEvents)*8)
		
			cat("s8")
			coxfitAfterYears=try(coxph(response ~ geno*env +	strata(Gender,na.group=TRUE) + frailty(CensusDivision), 
							subset = sort(c(theEvents, sample(theNoEvent, NtoKeep))),
  							data=datamat, control=coxph.control(iter.max=30,toler.chol=.Machine$double.eps)) ) 
			if(class(coxfitAfterYears)[1]=="try-error")
				coxfitAfterYears = list(coxfitAfterYears, coef=NA) 
			Ntries8 = Ntries8 + 1
		}
		Ntries1 = 0
		# if still singular matrix, keep equal numbers of cases as controls

while(any(is.na(coxfitAfterYears$coef)) & Ntries1 <= 100) {
	NtoKeep = round(length(theEvents)*1)
	
			cat("s1")
			coxfitAfterYears=try(coxph(response ~ geno*env +	strata(Gender,na.group=TRUE) + frailty(CensusDivision), 
							subset = sort(c(theEvents, sample(theNoEvent, NtoKeep))),
  							data=datamat, control=coxph.control(iter.max=30,toler.chol=.Machine$double.eps)) ) 
			if(class(coxfitAfterYears)[1]=="try-error")
				coxfitAfterYears = list(coxfitAfterYears, coef=NA) 
			Ntries1 = Ntries1 + 1
		}
		
		
		
		
		}
	retries[as.character(Dfollowup), Dcancer,] = c(Ntries, Ntries8, Ntries1)
		
			# p value is probability of z being above the observed z scores
			# 1- probability of z being below the observed
			 
    
      if("coefficients" %in% names(coxfitAfterYears)) {
      	pvalue[as.character(Dfollowup),,Dcancer]= 1 - pchisq((coxfitAfterYears$coefficients^2)/diag(coxfitAfterYears$var), 1)

      	estimates[as.character(Dfollowup),,Dcancer] <-
    			coxfitAfterYears$coefficients
    	
		standarderrors[as.character(Dfollowup),,Dcancer] = 
				diag(coxfitAfterYears$var)
		
		
	} else {
cat("\n warning, null coefficieint, probably too many ties\n")
pvalue[as.character(Dfollowup),,Dcancer]=0.5
	}
	 	
	 	} #end if there are events
  if ( !is.null(parameters$controlProp)   ) {
    datamat=datamatFull
  }
	} # end loop through followup
	
	if(verbose) cat("\n")

	} # end if there's someone with the interactions
	
	
} # end loop through cancers



# replace NA's with 1?
pvalue[is.na(pvalue)] = 1 

#dimnames(pvalue)[[3]][length(CommonCancer)] = "All"

attributes(pvalue)$coef = estimates
attributes(pvalue)$var = standarderrors
attributes(pvalue)$retries = retries
if(saveData) {
	attributes(pvalue)$data = datamat
}
	

return(pvalue)
}



#### floor(approx( cumsum(c(0, 0.2, 0.5, 0.3)), 1:4, runif(1))$y)
  ##which(as.logical(rmultinom(1, 1, c(0.2, 0.5, 0.3))))


##plot function starts here

seqPower = function(parameterName,  parSeq, 
   allcoefs, CommonCancer,SimulationTime, agegroup, lambdaF, 
  lambdaM, ageCancer, Enrollmentprobs,Followup, YesNoRate, NoYesRate) {

  if(! parameterName %in% names(allcoefs))
    warning(c(parameterName, names(allcoefs)))



# get data together in a format patrick's rewrite wants
# parameters
parameters = list(
	   geno=allcoefs["geno"], env=allcoefs["env"], "geno:env"=allcoefs["geno:env"],
		 size=allcoefs["size"], probGeno=allcoefs["probGeno"], 
		 probEnv=allcoefs["probEnv"], genoMiss= allcoefs["GenoMiss"],
		 envMiss = allcoefs["EnvMiss"], oneEnvPerCommunity =oneEnvPerCommunity,
     cancerMissRate=YesNoRate, falseCancerRate=NoYesRate, sdGroup=allcoefs["sdGroup"], 
		 sdPerson=allcoefs["sdPerson"], Enrollmentprobs = Enrollmentprobs, agerange = c(35, 69), 
		 Followup=Followup
)

populationData = getPopData()            

lambda = list("F"=lambdaF, M=lambdaM)    

parameters[[parameterName]] = parSeq[1]
result= Power(parameters, populationData, lambda, CommonCancer,  SimulationTime)

 
AllResult = array(NA, c(length(parSeq), dim(result)))
dimnames(AllResult) = c(list(paste(parameterName, "=" ,parSeq, sep="")), dimnames(result)) 
AllResult[1,,,,]=result

		        
for (i in 2:length(parSeq)){
parameters[[parameterName]] = parSeq[i]
result= Power(parameters, populationData, lambda, CommonCancer,  SimulationTime)
AllResult[i,,,,]=result

}

return(AllResult)

}

seqPowerList = function(varying,
	parameters, CommonCancer, SimulationTime, 
	populationData,deathFile, cancerRatesByType, verbose=F) {
# varying is a list, with elements having names being parameters in the model
# and the elements themselves being vectors of parameter values
	
if(!all(names(varying) %in% names(parameters)) )
    warning(c(names(varying), names(parameters)))

for(Dvarying in names(varying))
	parameters[[Dvarying]] = varying[[Dvarying]][1]


	
result= Power(parameters, populationData,  CommonCancer,  SimulationTime, deathFile, cancerRatesByType, verbose)


AllResult = array(NA, c(length(varying[[1]]), dim(result)))

fordimnames = as.data.frame(varying)
names(fordimnames) = names(varying)

for(Dvarying in names(varying))
	fordimnames[,Dvarying] = paste(Dvarying, fordimnames[,Dvarying], sep="=")
fordimnames = apply(fordimnames, 1, toString)
fordimnames = gsub(" ", "", fordimnames)

dimnames(AllResult) = c(
	list(fordimnames), 
	dimnames(result)) 

AllResult[1,,,,]=result


 
		        
for (i in seq(from=2, len=length(varying[[1]])-1, by=1) ){
	for(Dvarying in names(varying))
		parameters[[Dvarying]] = varying[[Dvarying]][i]
	
	result= Power(parameters, populationData, CommonCancer,   SimulationTime, deathFile, cancerRatesByType, verbose)

	AllResult[i,,,,]=result
}

return(AllResult)

}

getPopData = function(file="../data/On129Communities.csv") {



# read in and format data of population by CSD
thedata = read.csv(file, sep=",", header=T,
	as.is=T)
	
	
populationData=list(
	F=cbind(thedata$Female__35_to_39_years,thedata$Female__40_to_44_years,
			thedata$Female__45_to_49_years,
			thedata$Female__50_to_54_years, thedata$Female__55_to_59_years,
			thedata$Female__60_to_64_years,thedata$Female__65_to_69_years),

	M=cbind( 
		thedata$Male__35_to_39_years,thedata$Male__40_to_44_years,
		thedata$Male__45_to_49_years,
		thedata$Male__50_to_54_years,thedata$Male__55_to_59_years,
		thedata$Male__60_to_64_years,thedata$Male__65_to_69_years) 
		
	
	)
	
colnames(populationData$F) = 	colnames(populationData$M) = 
	as.character(seq(35, 65, by=5))
	
#simulate gender
populationData$all =cbind(
	propFemale=thedata$Female__35_69_years/thedata$Total__35_69_years,
	total=thedata$Total__35_69_years)  
rownames(populationData$all) = rownames(populationData$M) =
	rownames(populationData$F) = thedata$GeoName	

populationData

}
