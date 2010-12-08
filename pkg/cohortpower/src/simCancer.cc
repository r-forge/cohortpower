#include <Rmath.h>

#define MIN(a,b) ((a)>(b)?(b):(a))

extern "C" void simCancer(double *Ri, double *start, double *death, 
  double *result, int *Nindiv, 
  double *ageSeq, double *lambdaSeq, double *diffLambda, 
  double *thisAgeSeq, double *thisLambda, int *Nlambda) {

  // called from the R function simCancer. 

  register int Dindiv, Dage, Dlambda;
  double mui, Ncancers, thisCancer, tryCancer;

   
  for(Dindiv=0;Dindiv<*Nindiv;++Dindiv) {
    //get mui
    //find first age bracket greater than cohort joining age
    for(Dage=0;start[Dindiv]>ageSeq[Dage];++Dage);

    //found the age bracket which contains the start age
    thisAgeSeq[0] = start[Dindiv];
    thisLambda[0] = 0 ;
    thisAgeSeq[1] = ageSeq[Dage];

    //if the person dies before the next age bracket starts, 
    //compute mu up to the death age
    if(death[Dindiv] < ageSeq[Dage]) {
      mui = thisLambda[1] = (death[Dindiv] - start[Dindiv]) * lambdaSeq[Dage-1];
    //otherwise add mu up to the next age bracket
    } else {
      mui = thisLambda[1] = (ageSeq[Dage] - start[Dindiv]) * lambdaSeq[Dage-1];
      //loop through all age brackets contained within the individual's life time
      Dage++;
      for(Dlambda = 2;death[Dindiv]>ageSeq[Dage];++Dlambda,++Dage) {
        thisAgeSeq[Dlambda] = ageSeq[Dage];
        mui += diffLambda[Dage-1] * lambdaSeq[Dage-1];
        thisLambda[Dlambda] = mui;
      } 
      thisAgeSeq[Dlambda] = death[Dindiv];
      //add the relevant proportion of the final age group
      mui = thisLambda[Dlambda] = 
        mui + lambdaSeq[Dage-1] * (death[Dindiv] - ageSeq[Dage-1]);
    }

   Ncancers = rpois( (mui * Ri[Dindiv] ) );

   if(Ncancers) {
    //loop through cancer incidents (if any) and keep the earliest one
    thisCancer = runif(0,1);
    for(Dage=1;Dage < Ncancers;++Dage) {
      tryCancer = runif(0,1);
      thisCancer = MIN(thisCancer, tryCancer) ;
    }
    //rather than scaling the lambda by mu, multiply the uniform 0 1 by mu
    thisCancer = thisCancer * mui;
    // find the cdf range containing thisCancer
    for(Dage=0;thisLambda[Dage] < thisCancer;++Dage);
    // assign cancer as the quantile of the cdf of thisCancer
    result[Dindiv] = thisAgeSeq[Dage-1] + (death[Dindiv]-thisAgeSeq[Dage-1]) *
         (thisCancer - thisLambda[Dage-1]) / (thisLambda[Dage] - thisLambda[Dage-1]);  
    
  } else {
    //assign cancer age 200
    result[Dindiv] = 200 ;   
  }

  }


}
