
#if(F) {
#CommonCancer=c("Breast", "Prostate","Colo","Lung", "Stomach")
#CommonCancer=c("Colo", "Stomach")

#Nsim=10


#}


CaseCI=function(parameters, lambda, populationData, deathFile, cancerRatesByType, CommonCancer,verbose=F,
theGender=c("F","M"), Nsim, Percentile=c(0.025, 0.975)) {

if(is.null(cancerRatesByType)){
CommonCancer = "All"
}else{
if(is.null(CommonCancer)){
CommonCancer=	c(colnames(cancerRatesByType[[1]]) [ !colnames(cancerRatesByType[[1]])== "Total" ],"[[:alpha:]]")
}
}
carray=array(NA, c(length(CommonCancer),length(parameters$Followup),length(theGender),Nsim))
dimnames(carray)= list(c(CommonCancer), as.character(parameters$Followup),theGender,NULL)


for (Dsim in 1:Nsim){
carray[,,,Dsim]=CaseNum(parameters, lambda, populationData, deathFile, cancerRatesByType, CommonCancer, theGender=c("F", "M"), verbose=F)
                            }

dimnames(carray) = list(c(CommonCancer), as.character(parameters$Followup),theGender)
EstiCase= apply(carray, c(1,2,3), function(qq) mean(qq))
######### For Confidence Interval
EstiCI = array(NA, c(length(CommonCancer), length(parameters$Followup), length(theGender), length(Percentile)))
dimnames(EstiCI) = list(c(CommonCancer), as.character(parameters$Followup),theGender, as.character(Percentile))

EstiCI = apply(carray, c(1,2,3), function(ci) (quantile(ci,probs=Percentile) ))

#save(EstiCase, file="EstiCase.Rdata")
#save(EstiCI, file="EstiCI.Rdata")
 return(list(CI = EstiCI, cases = EstiCase)
)
}




CaseNum = function(parameters, lambda, populationData, deathFile, cancerRatesByType, CommonCancer, theGender=c("F", "M"), verbose=F){


datamat = simCohort(parameters, lambda, populationData, deathFile=deathFile, cancerRatesByType)

size=parameters$size
Followup=as.character(parameters$Followup)

if(is.null(cancerRatesByType)){
datamat$cancerRatesByType=datamat$Cancer
}

NumOfCases=array(NA, c(length(CommonCancer),length(Followup),length(theGender)))
dimnames(NumOfCases)= list(c(CommonCancer), Followup,theGender)


for(Dcancer in CommonCancer) {
  
	if(verbose) cat(Dcancer)
	thisCancer = rep(F, parameters$size)
# find the cancer we're interested in
      if(is.null(cancerRatesByType)){
      thisCancer=datamat$Cancer
      } else{
	    thisCancer[grep(Dcancer, datamat$cancerRatesByType)] = T
      }
                            
	for(Dfollowup in Followup) {
	if(verbose) cat(Dfollowup)

	 	# set thisEventTime equal to age at end of followoup period
		ageEndFollowup = datamat$Age + as.numeric(Dfollowup) - datamat$Enrollment
		hadEvent = datamat$Event < ageEndFollowup

		# eventType=T if the cancer of interest occurred during the followup period
		datamat$thisEventType = hadEvent &  thisCancer
    NumOfCases[Dcancer,Dfollowup,1]=sum(as.numeric(datamat$thisEventType& datamat$Gender=="F"))
    NumOfCases[Dcancer,Dfollowup,2]=sum(as.numeric(datamat$thisEventType&datamat$Gender=="M"))

		}
   }
NumOfCases
}