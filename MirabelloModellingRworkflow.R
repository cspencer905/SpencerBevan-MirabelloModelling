# Mirabello PPM #
# Spencer and Bevan 2016 #

library(spatstat)
library(raster)
library(rgdal)
library(maptools)
library(arm)
#library(dismo)

# import data:
setwd("X://MirabelloModellingPaper")

#image and raster layer of dem
#dem <- as.im(as(raster("data/dem_surveyed.asc"), 'SpatialGridDataFrame')) 	#dem as image
#dem_surv <- raster("data/dem_surveyed.asc")					#dem as raster

# polygon layers for survey boundaries and coast
coastc <- readOGR(dsn="data/shapefiles/mirab_polyc.shp", layer="mirab_polyc")	#shape polygon
study2 <-readOGR(dsn="data/shapefiles/surv_total.shp", layer="surv_total")		# study area polygon (dem coastline)
study2.proj <-proj4string(study2)

# All period site data
sitepts <- read.csv("data/csvs/MirabSiteData.csv", header=T)						#site points
sites_buff <- readOGR(dsn="data/shapefiles/sites_buff/sites_buff.shp", layer="sites_buff")	#site polygons

# create ppp with site-size marks
#i <- "EMIII-MMIA" 			# choose appropriate period
#i <- "proto" 
i <- "neo"
#i <- "postpal"
sites <- sitepts[sitepts$palatial==i,] 
#sites <- sitepts[sitepts$phase==i,] 	# choose for EMIII-MMIA only

pppw <- ppp(sites$centroid_x, sites$centroid_y, marks=sites$cat_area, window=as.owin(study))	#ppp of selected sites with site size marks
ppps <- unmark(pppw)										#ppp without marks

# import environmental variables as images (SCALED [0,1])
#sqrtelev <- as.im(as(raster("data/covs/for_ppm/sqrtelev.asc"), 'SpatialGridDataFrame'))
#logslope <- as.im(as(raster("data/covs/for_ppm/logslope.asc"), 'SpatialGridDataFrame'))
#asp180 <-  as.im(as(raster("data/covs/for_ppm/asp180.asc"), 'SpatialGridDataFrame'))
#flat <-  as.im(as(raster("data/covs/for_ppm/flat.asc"), 'SpatialGridDataFrame'))
#sqtwi <- as.im(as(raster("data/covs/for_ppm/sqtwi.asc"), 'SpatialGridDataFrame'))
#distcst <-  as.im(as(raster("data/covs/for_ppm/distcst.asc"), 'SpatialGridDataFrame'))
#distlime <-  as.im(as(raster("data/covs/for_ppm/distlime.asc"), 'SpatialGridDataFrame'))
#ridge <-  as.im(as(raster("data/covs/for_ppm/ridge.asc"), 'SpatialGridDataFrame'))
#channel <-  as.im(as(raster("data/covs/for_ppm/channel.asc"), 'SpatialGridDataFrame'))
#nearview <-  as.im(as(raster("data/covs/for_ppm/nearview.asc"), 'SpatialGridDataFrame'))
#farview <-  as.im(as(raster("data/covs/for_ppm/farview.asc"), 'SpatialGridDataFrame'))

#covlist <- list(sqrtelev, logslope,asp180, flat, sqtwi, distcst, distlime, ridge, channel, nearview, farview) # updated list
#names(covlist) <- c("sqrtelev", "logslope", "asp180", "flat", "sqtwi", "distcst", "distlime", "ridge", "channel", "nearview", "farview")

#####

source("r/ppmfunc-July2017.R")

# create objects for input parameters for modelling

# choose from one of the following:
study_abs1 <- sites_buff[sites_buff$palatial!=i,] 	# subtract sites of selected period from all sites polygon shapefile
#study_abs1 <- sites_buff[sites_buff$phase!=i,] 	# for EMIII-MMIA only 

study_abs2 <- study2 - sites_buff 				# no evidence of activity in any period in study area

study_pre <- sites_buff[sites_buff$palatial==i,] 	# site polygons for selected period
#study_pre <- sites_buff[sites_buff$phase==i,]  	# for EMIII-MMIA only 

# specify covariates in logistic regression model (environmental raster files as im objects)
modf <- PPP ~ sqrtelev + logslope + asp180 + flat + sqtwi + distcst + distlime + ridge + channel + nearview + farview
modt <- ~ sqrtelev + logslope + asp180 + flat + sqtwi + distcst + distlime + ridge + channel + nearview + farview

# specify number of presence and absence points to be sampled
npts <- ceiling((sum(sites$cat_area)/1000)) 		# WHEN: area dependency sampling
#npts <- npoints(ppps) 						# WHEN: n sample points = n sites

nsim <- 1000								# select number of simulations to be iterated

### run models
## function (adapted from logRMCr) to produce dataframe with coeffs and intercepts from each spsample data and dummy point ppm

#site locations for selected period compared against site locations in other periods
#modabs1 <- condlogRMCppm(study_pre,study_abs1,npts=npts, areadepend=TRUE, nsim=nsim,trend=modt,covariates=covlist, verbose=T, testmode=FALSE)

#site locations for selected period compared against empty landscape within surveyed area
#modabs2 <- condlogRMCppm(study_pre,study_abs2,npts=npts, areadepend=TRUE, nsim=nsim,trend=modt,covariates=covlist, verbose=T, testmode=FALSE)

#save(modabs1, file="r/mods/neomodabs1.RData") # change file name depending on period
#save(modabs2, file="r/mods/neomodabs2.RData")
#write.csv(modabs1, file="r/plots/neomod1.csv")	  #save output as backup csv file
#write.csv(modabs2, file="r/plots/neomod2.csv")

#### sort dataframe, change column names to letters
modabs1[1] <- NULL
modabs1 <- reworkMCr(modabs1,modf)
modabs2[1] <- NULL
modabs2 <- reworkMCr(modabs2,modf)

# frequency of chosen cov in fitted model
abs1probs <- data.frame(Covariates=names(modabs1), Prob=(nsim - as.numeric(colSums(is.na(modabs1)))) / nsim)
abs2probs <- data.frame(Covariates=names(modabs2), Prob=(nsim - as.numeric(colSums(is.na(modabs2)))) / nsim)

# boxplot of results
y <- 15
boxplot(modabs2, pch=19, ylim=c(-(y),y), cex=0.4, whisklty="solid",main="sampled from general survey landscape", cex.main=0.8)
abline(h=0,lty="dotted")
text(x=c(1:length(abs2probs$Prob)), y=y, labels=as.character(round(abs2probs$Prob,2)), cex=0.75, font=3)
boxplot(modabs1, pch=19, ylim=c(-(y),y), cex=0.4, whisklty="solid",main="sampled from site locations in other periods", cex.main=0.8)
abline(h=0,lty="dotted")
text(x=c(1:length(abs1probs$Prob)), y=y, labels=as.character(round(abs1probs$Prob,2)), cex=0.75, font=3)

###

# Second-order modelling using inut from first-order PPM Models
# Area-dependent sampling example

# load covs
#sqrtelev <- raster("data/covs/for_ppm/sqrtelev.asc")
#logslope <- raster("data/covs/for_ppm/logslope.asc")
#asp180 <-  raster("data/covs/for_ppm/asp180.asc")
#flat <-  raster("data/covs/for_ppm/flat.asc")
#sqtwi <- raster("data/covs/for_ppm/sqtwi.asc")
#distcst <-  raster("data/covs/for_ppm/distcst.asc")
#distlime <-  raster("data/covs/for_ppm/distlime.asc")
#ridge <-  raster("data/covs/for_ppm/ridge.asc")
#channel <-  raster("data/covs/for_ppm/channel.asc")
#nearview <-  raster("data/covs/for_ppm/nearview.asc")
#farview <-  raster("data/covs/for_ppm/farview.asc")

# mask to study area
s <- stack(sqrtelev, logslope,asp180, flat, sqtwi, distcst, distlime, ridge, channel, nearview, farview) 
s <- mask(s, study2)

# Generate first-order trend surface from PPM results
#load and clean mod results
load("r/mods/neomodabs2_s1000_study2.RData")
intrcpt <- modabs2[,1]
modabs2[1] <- NULL
modabs2 <- reworkMCr(modabs2,modf)

# calc mean for cov coeff
for (i in 1:ncol(modabs2)){
	if (i==1){
	cov_avg <- mean(modabs2[,i],na.rm=T) }
	else {
	avg <- mean(modabs2[,i],na.rm=T)
	cov_avg <- c(cov_avg,avg)}}

# calc mean for intercept
int_avg <- mean(intrcpt)

# turn insignif covs to zeros
insig <- c(3,4,5,6,9)			# input manually
for (j in 1:length(insig)){
	cov_avg[insig[j]] <- 0}

# log-odds surface based on the average coefficient values from simulations
logodds <- int_avg+(s[[1]]*cov_avg[1])+(s[[2]]*cov_avg[2])+(s[[3]]*cov_avg[3])+(s[[4]]*cov_avg[4])+(s[[5]]*cov_avg[5])+(s[[6]]*cov_avg[6])+(s[[7]]*cov_avg[7])+(s[[8]]*cov_avg[8])+(s[[9]]*cov_avg[9])+(s[[10]]*cov_avg[10])+(s[[11]]*cov_avg[11])
relprob <- (exp(logodds))/(1+(exp(logodds)))

# raster surface as im object for input in envelope()
lambda <- as.im(relprob)

# make ppp of simulated households using area dependency sampling method
datarea <- sites_buff[sites_buff$palatial==i,]
npts <- ceiling((sum(sites$cat_area)/1000))

for (j in 1:length(datarea)){
	rpt <- spsample(datarea[j,],n=1,type="random", iter=20)
		if (j==1){ rdat <- rpt }
		else { rdat <- rbind(rdat,rpt) }}
if (length(rdat) < npts){
	rpts <- spsample(datarea,n=(npts-length(rdat)), type="random", iter=20)
	rdat <- rbind(rdat,rpts) }


# creation of simulated point patterns from intensity surface
rlmbd <- rpoispp(lambda=lambda,nsim=999)

# overwrite oberservation window to match ppp window
for (i in 1:length(rlmbd)){
	rlmbd[[i]]$window <- ppps$window }

# pair correlation function using cov-simulated ppp of absence points based on trend surface
#hldsPCFcov <- envelope(as.ppp(coordinates(rdat),W=as.owin(study2)),pcf,simulate=rlmbd,nsim=999,nrank=25)


#plotting of pcf results
load(file="r/pcf/neo_pcfsamp_data.RData")			# load pcfdat
load(file="r/pcf/neo_hldsPCF999cov.RData") 		# load pcf results
hldsPCFcov <- neo_hldsPCF
x2 <- data.frame(hldsPCFcov)
r <- x2[,1]

#is clustered when: min>x2[1:513,5]
#is dispersed when: max<x2[1:513,4]

plot(hldsPCFcov,ylim=c(0,4),legend=FALSE,xlab="",ylab="",bty="l",main="")
title(ylab=expression(italic("g(r)")), line=2,cex.lab=1)
title(xlab=expression(italic("inter-site distance")), line=2,cex.lab=1)
rect(r[2],x2[2:29,3],r[28],4.15,col="#FCBDB6", border=NA)
rect(r[81],x2[81:84,3],r[84],4.15,col="#FCBDB6", border=NA)
rect(r[91],x2[91:95,3],r[95],4.15,col="#FCBDB6", border=NA)
rect(r[259],-0.15,r[261],x2[259:261,3],col="#BCEEFB", border=NA)
rect(r[268],-0.15,r[277],x2[268:277,3],col="#BCEEFB", border=NA)
rect(r[346],-0.15,r[361],x2[346:361,3],col="#BCEEFB", border=NA)
rect(r[379],-0.15,r[394],x2[379:394,3],col="#BCEEFB", border=NA)
rect(r[428],-0.15,r[429],x2[428:429,3],col="#BCEEFB", border=NA)
plot(hldsPCFcov,ylim=c(0,4),legend=FALSE,add=TRUE,bty="l",main="")
polygon(c(r,rev(r)), c(max,rev(min)),col="#666666C8",border=NA)
lines(r,x2[,2],col="black")
legend(2000,4, legend=c("Observed statistic","Expected random statistic"),lty=c(1,2),col=c("black","red"),bty="n",cex=0.9)
legend(1990,3.74,legend=c("Simulation envelope","Household sampling variation","Significant segregation","Significant clustering"),pch=0,fill=c("grey","gray40","#BCEEFB","#FCBDB6"),col=NA,border=NA,bty="n",cex=0.9)

##



