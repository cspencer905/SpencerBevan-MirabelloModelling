
#modified from sptools (Bevan 2011)
#Spencer - October 2016

condlogRMCppm <- function(datarea, dumarea, npts, areadepend=TRUE, nsim, trend, covariates, verbose=TRUE,testmode=FALSE,...){
 
  #Run loop
  cat('Sample run ')
  for (i in 1:nsim){

    # Message
    if (verbose){
      if (i==nsim){
        cat(paste(i,".\n",sep=""))
      } else {
        if (isWholeNumber(i/10)){
          cat(paste(i,'...\n', sep=""))
        } else {
          cat(paste(i,"...", sep=""))           
        }
      }
    }

# spsample() picks presence points from inside area specified by the datarea arg (in this case either within site polygons from other periods, or in the wider landscape)
# if areadepend=TRUE, then npoints arg must also corrspond to the density of points wanting to be sampled	
   if (areadepend){
	for (j in 1:length(datarea)){
		rpt <- spsample(datarea[j,],n=1,type="random", iter=20)
			if (j==1){ rdat <- rpt }
			else { rdat <- rbind(rdat,rpt) }
	}
	if (length(rdat) < npts){
		rpts <- spsample(datarea,n=(npts-length(rdat)), type="random", iter=20)
		rdat <- rbind(rdat,rpts)
	}
   }
   else { rdat <- spsample(datarea, n=npts, type="random", iter=20) }		
# if npts=number of site points, make areadepend=FALSE

# spsample() picks absence points from dumarea arg
   rdum <- spsample(dumarea,n=npts, type="random",iter=20) #increase iter necessary - likely due to n dummy points choice against sample area size

# pass a set of pre-cooked and randomised data and dummy points to be used as presence/absence data in logistic regression, so very close to standard glm
   myq <- quadscheme.logi(as.ppp(coordinates(rdat),W=as.owin(study2)),as.ppp(coordinates(rdum),W=as.owin(study2)))

# with print out distributions of each data/dummy point sample for verification
   if (testmode){
		cat(paste("n data: ",npoints(myq$data),'\n',sep=""))
		cat(paste("n dummy: ",npoints(myq$dummy),'\n',sep=""))
		dev.new()
		plot(myq$data, main=c("sim no. ",i),pch=19,cex=0.3)
		plot(myq$dummy, pch=19, cex=0.3, col="red", add=TRUE) }

# logistic regression on specified presence/absence points from quadrature scheme, using stepAIC method to select covariates to build first order trend
   cov.ppm <-step(ppm(myq,trend=trend, interaction=NULL, covariates=covariates, method="logi"))
# produce dataframe with coeffs and intercepts from each simulation
   coef <- data.frame(t(coef(cov.ppm)))
   if (i==1){ simcoef <- coef} else { simcoef <- merge(simcoef,coef,all=TRUE)}}

return(simcoef)

}


##

isWholeNumber <- function(x, tol = .Machine$double.eps^0.5) {
  #Checks whether a value is a whole number (used in nnhistMC reporting)
  abs(x - round(x)) < tol
}

reworkMCr <- function(modres, modf){
  
  covs <- data.frame(Covariates=all.vars(modf)[2:length(all.vars(modf))], Letter=LETTERS[1:length(all.vars(modf))-1])
  
  for (i in 1:length(covs$Covariates)){
    if (!covs$Covariates[i] %in% names(modres)){
      modres$tmp <- NA
      names(modres)[names(modres)=="tmp"] <- as.character(covs$Covariates[i])
    }
  }
  nm <- data.frame(Covariates=names(modres))
  mm <- merge(nm, covs, by="Covariates", sort=FALSE)
  names(modres) <- mm$Letter
  modres <- modres[ ,c(sort(names(modres)))]
  
  return(modres)
}

##