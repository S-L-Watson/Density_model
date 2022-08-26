         
#-----What:--------- BRT Analysis, Presence Pseudo Absence model for benthic invertebrates  -----------###
# Step 1  - Run the Presence Absence Model
# Step 2  - Run the Abundance model 
# Step 3  - Combine the outputs from each model, scale or normalize

#--- Supporting script for publication -- Watson, Stephanie Louise and Stephenson, Fabrice and Conrad, Pilditch A. and Lundquist, Carolyn, 
# Improving Predictions of Coastal Benthic Invertebrate Occurrence and Density Using a Multi-Scalar Approach.
# Available at SSRN: https://ssrn.com/abstract=4142118 or http://dx.doi.org/10.2139/ssrn.4142118

                                       ################ ---- START ----######################

# load packages
require(raster);require(dismo);require(gbm); require(sp); require(pROC);require(FinCal);require(ggplot2);require(dplyr)
require(rgdal); require(corrplot); require(devEMF);require(rcompanion); require(DescTools); require(BBmisc); require(tools)
                    
#----LOADING OF DATA STARTS----
#---1---LOAD SAMPLE DATA (species data, clipped to mask area) ####
setwd("C:/Users/sw277/Desktop/R/PhD/Invert_data/Oct_20/abundance")
Ama <- read.csv("AMA_TGB_ABU.csv")
Ama$pa <- 1 # add new column for presence absence 
Ama <- data.frame(Ama$POINT_X, Ama$POINT_Y, Ama$Abu_1km, Ama$pa) # bind relevant columns
colnames(Ama) <- c("X", "Y", "Abu", "pres") # rename - count data now Abundance "Abu"
head(Ama) 
                
#---2---LOAD ENVIRO PREDS ##########
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB/Env_preds") 
f250m <- list.files(getwd()) # 14 env preds
imp.var <-c("X","Y","FID",c(file_path_sans_ext(f250m)))
ras250m <- lapply(f250m, raster)
f250m <- file_path_sans_ext(f250m)
# plot(ras250m[[2]]) # check through layers make sure same extent and projection
env250m_rs <- raster::stack(ras250m) # calling the function directly within the package
# #( if using multiple packages, call them as needed, as different packages can block commands if they have the same name)

# #--3--- CREATE SPDF: Turn raster stack and spatial points data into a dataframe for model prediction#####
env250m_spdf <- rasterToPoints(env250m_rs, spatial = T) # turns data into a large spatial points dataframe
# Extract the coordinates, for the spatial points df, and bind it with the env data from the raster stack,
# and save as a new dataframe
env250m_df <- data.frame(env250m_spdf@data, X =env250m_spdf@coords[,1],Y = env250m_spdf@coords[,2])
head(env250m_df)
# 
# #---4--SET UNIQUE FID, EXTRACT X & Y( row and column numbers, not lat and long) data from raster stack ####
#env250m_spdf@coords # a command to see the x and y from the raster data
env250m_df$FID <- cellFromXY(env250m_rs, env250m_df[,c("X", "Y")])
env250m_df <- env250m_df[,c(35:36,1:34,37)] # re-shuffle the order of columns, to put x and y at front
str(env250m_df) # check everything is named correctly
# save for use next time to speed things up
TGB.all.env.vars.train <- env250m_df
setwd("C:/Users/sw277/Desktop/R/PhD/Environ_preds/Environ_preds/250m_standardised")
save(TGB.all.env.vars.train, file =  "TGB.all.env.vars.train.Rdata")

#---4b---LOAD ENV_VAR_DF ONCE CREATED - from steps 1-4 above----
setwd("C:/Users/sw277/Desktop/R/PhD/Environ_preds/Environ_preds/250m_standardised")
load("TGB.250m.all.env.vars.train.Rdata") # load all env variables DF

#5--- EXTRACT ENV.VAR INFO where presences occur ----
Ama$FID <- cellFromXY(env250m_rs, Ama[,c("X", "Y")])# from the object (raster stack), return the raster cells (row +cols), in which
#the coordinates from Ama (in the form of X and Y) are found within the raster, and place them into newly created  FID column 
Ama_F <- cbind(Ama ,raster::extract(env250m_rs, Ama[,c("X", "Y")])) # extract all env data where Ama occurs in env var raster stack

#-6---AGGREGATE ABUNDANCE to 1km --
#/ group by FID 
#we bring in the abundance column from our Raw data which is in count format - sometimes calculated abundance, and then we create new 
# column called mean abundance (mean_abu) - multiple records are aggregated per cell 
Ama_F <- Ama_F %>% # 
group_by(FID) %>%
summarize(mean_abu = (sum(Abu)/sum(pres))) %>% # make sure that multiple points are calculated as a proportion for that raster cell /geographic sample area.
select(mean_abu = mean_abu, FID) %>%
left_join(Ama_F,mean_abu ,by= "FID") # combine the newly created column with existing dataset like cbind
# shuffle order so that mean_abu,pres, wg, and FID are last columns (remove abundance as its no longer needed)
Ama_F.M2 <- as.data.frame(Ama_F[,c(3:4,1,7:13,6,2)]) # HERE - change to select 7 preds  # omitting the 5th column as dont need Abu

#6b SAVE/LOAD FINAL SPECIES DATA---- 
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
save(Ama_F.M2, file = "Ama_F.M2.Rdata")
rm(Ama_F.M2)

#LOAD FINAL SPECIES DATA 
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
#species data from steps 5 & 6 above
load("Ama_F.M2.Rdata")
# Important Env vars, extracted from step 4 above (already refined from gbm simplify)
# define the important vars per species and model type then subset
imp.var <-c("X","Y","FID",c(f250m))  # create character index from env var files
TGB.250m.all.env.vars.train <- TGB.250m.all.env.vars.train[, imp.var]

# INSPECT DATA DISTRIBUTION to see if transformation needed # (# 3 Options)----
#library(rcompanion)
#1 - DO NOTHING
# sometimes the natural way in which the data comes in, is more evenly distributed than it 
# would be if you log transformed it. Therefore its better to apply no transformation
# and simply use the Gaussian model option. to view the distribution use 
plotNormalHistogram(Ama_F.M2$mean_abu)

#--2 - Round abundance data and use Poisson family in gbm 
# if heavily skewed (used Poisson family in gbm) -and skip step #3 (the log transform) below 
#plotNormalHistogram(Ama_F.M2$mean_abu) 

#3 IF using Guassian -( more normally distributed) Log Transform abundance data and use Guassian family in gbm
# https://rcompanion.org/handbook/I_12.html
#Ama_F.M2$mean_abu <- log1p(Ama_F.M2$mean_abu +1) #  a
#plotNormalHistogram(Ama_F.M2$mean_abu) 

#--MODEL OBJECTS(-setup objects to save information from model runs to compare metrics with diff ENV Preds)####
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB/1st_run_all_env_preds")
n.boot <- 100
length.map <- nrow(TGB.250m.all.env.vars.train)
boot_mat <- array(0, c(length.map,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)
deviance_mat <- array(0, c(n.boot,6)) # for saving model fit metrics (length of bootstaps with 6 columns)
influence_mat <- array(0,c(7,n.boot)) # for saving information on predictors (15 rows, 100 columns) # ths number of rows here, 
#has to match the number of env preds used, or will throw out an error.

# test 1 loop of your model first, check metrics, before bootstrapping  

#Setup for partial dependance plots
PredMins <- apply(na.omit(Ama_F.M2)[,f250m],2,min)
PredMaxs <- apply(na.omit(Ama_F.M2)[,f250m],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(f250m)))) # store a DF with a sequence of env ranges from min to max, from all raw values in f250ms
names(EnvRanges) <- f250m
for (i in c(1:length(f250m))) {
   EnvRanges[,i] <- seq(PredMins[f250m[i]],PredMaxs[f250m[i]],length = 100)
}
# for each model run (bootstrap) store the results of environmental ranges, through the range
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot)) # #3D vector ( an array) with columns given same names as f250ms
dimnames(PD)[[2]] <- f250m
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="") # the 3rd dimension is names the rep number of the bootstrap

#7---CREATE TRAINING AND EVAL DATA ----####
# * NOTES - if dataset has <100 samples, use all data and DONT split into training and eval
# * Only log data if normally distributed, otherwise leave as is and try and run Guassian - as Possion only works with whole numbers
# * if model only runs sometimes, will need to leave while statement on to complete through bootstrap

counter = 0
for (i in 1:n.boot){ # bootstrap loop for M2 (Abu model) 
#(1) create randomised sample of presence data 
# ind_pres <- sample(seq_len(nrow(Ama_F.M2)), size = floor(0.85 * nrow(Ama_F.M2))) # create an index vector of rows for 75% of presence data
# # using the above index, subset from the orgional dataset to create training and evaluation presences 
# pres_train <-Ama_F.M2[ind_pres, ] # 75% - IF dataset is less than 100 points change to 100% and dont have eval portion
# pres_eval <-Ama_F.M2[-ind_pres, ] # 25%

#7A--- M2 TRAINING MODEL for Abundance / export to rast ##### (use model 8A if small sample size <100 otherwise use 8b) ----
# #gbm.step - Fits a gbm model to one or more response variables, using cross-validation to estimate the optimal number of trees.
NoTrees <- 0
lr <- 0.01
           
while(NoTrees < 3000) {   # comment out if model runs fine - only needed for small sample size / difficult model that stops

   M2<-try(gbm.step(data= Ama_F.M2, gbm.x = c("BPI_fine","Carbonate","Ebed","Gravel","Mud","POCFlux","Slope","VGPM"), #	 (or column number of env preds to be used in analysis)
                        gbm.y = 3, # column number of abundance data (aka the response variable)
                        family = "gaussian", tree.complexity = 2, #(aka Interaction Depth) - play around with settings to obtain 1500-3000 trees
                        learning.rate = lr, bag.fraction = 0.5, n.folds=5,  # bag fraction between 0.5-0.75 is good, reduce n.folds down to 5 if model not converging
                        max.trees = 10000, step.size = 25, plot.main = T, tolerance.method = "fixed", # change 'plot main' to T when doing initial runs
                        tolerance = 0.01, verbose = T)) # switch this to T to visualize AUC /tree number for first run before looping ( then switch to F)
   NoTrees <- M2$n.trees
   if(is.null(NoTrees)){
      NoTrees <- 1
      lr <- lr *0.8
   }else
   {next}
}
   
   #partial dependance plots
for (j in 1:length(f250m)){
   grid <- gbm::plot.gbm(M2, i.var = c(paste(f250m[[j]])), return.grid = T)
   PD[,f250m[[j]],i] <- loess(grid$y~EnvRanges[,f250m[[j]]])$y # LOESS - local Polynomial Regression Fitting 
}  #allows for the extraction of best model fit, within the whole range of env preds used

#--internal model fit metrics
int.null.deviance <- M2$self.statistics$mean.null # trying to explain the deviance.If it is 0 it cannot explain any of the deviance
int.residual.deviance <- M2$cv.statistics$deviance.mean 
deviance_mat[i,1] <- (int.null.deviance-int.residual.deviance)/int.null.deviance # internal Deviance explained

#-- RMSE  & P.CORR -calculate the modeled and observed values for RMSE 
m = M2$fit # mean model fit
o = Ama_F.M2$mean_abu # mean abundance observed from data frame
#plot(o,m) # plot observed vs modeled obs  - blank out for bootstrap
deviance_mat[i,5] <- RMSE(m,o)   #True
deviance_mat[i,6] <- cor(o,m) # pearsons correlation between model fit and mean abundance 0-1. 
#Correlation coefficient between X and Y.if it is between >0.5 then it is a strong correlation

# #model fit comparison with withheld evaluation data (only if split between training and eval >100)
pred <- predict.gbm(M2, pres_eval, n.trees = M2$gbm.call$best.trees, type = "response")
ext.residual.deviance <- calc.deviance(pres_eval$pres, pred, family = "gaussian" ,calc.mean=T)
ext.null.deviance <- calc.deviance(pres_eval$pres,family = "gaussian", rep(mean(pres_eval$pres),nrow(pres_eval)), calc.mean=T)
deviance_mat[i,2]<-(ext.null.deviance - ext.residual.deviance)/ext.null.deviance  # deviance explained for eval data
cor(pred, eval.df)  # pearsons correlation (only if >`100 data points`)

#--Internal AUC
deviance_mat[i,3] <- M2$cv.statistics$correlation.mean

#--Env pred contribution
M2_contrib <- as.data.frame(M2$contributions) 
env_var_ord <- M2_contrib[order(M2_contrib$var),]
influence_mat[,i]<-env_var_ord[,2]
# env_var_ord

# predict spatially to map/study area (change out to smaller prediction area, than size of training tiffs if needed)
boot_mat[,i] <- round(predict.gbm(M2, TGB.250m.all.env.vars.train, n.trees = M2$gbm.call$best.trees, family = "bernoulli", type = "response"), digits = 2)
 
counter = counter +1
print(counter)
}
#------------------------------------END OF MODEL LOOP #-------------------------------------------####
summary(M2)

#-8---CREATE OUTPUTS & save data-####
# save created data
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
save(boot_mat, file="M2_boot_mat.Rdata")
save(deviance_mat, file="M2_dev_mat.Rdata")
save(influence_mat, file="M2_inf_mat.Rdata")

#-- PARTIAL DEPENDENCE PLOTS #----
#require(BBmisc)  #= normalize

PD1 <- normalize(PD, method = "range", range = c(0, 1), margin = 1) # normalize the Y axis on PD plots to be on same scale 
# plot PDs + 95 PI and save to file
emf(file = "PD_Ama.mean.TGB.emf", emfPlus = FALSE) 
par(mar=c(4, 2, 1, 1)) # PAR = margins / change layout of plot grids
par(mfrow=c(3,3))
#explanation: Go through the imp vars, and for each value row, find the corresponding value in the environmental ranges, and plot 
# the mean predictor response for that value. 
for (i in c(1:length(f250m))) { # or f250m
        plot(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,mean), col = 1,type= 'l', # this is a 3D vector [(rows, columns, depth /the rest)]
        xlab = paste(f250m[i]), 
        ylab = '',
        ylim = c(min(PD1[,,]), max(PD1[,,]))) 
  
# SOLUTION fOR RANGE IN SHADED POLYGON# in a DF
 maxs <- data.frame(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,max))
 mins <- data.frame(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,min))
 means <- data.frame(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,mean))
 DF3 <- data.frame(mins[,1], mins[,2],maxs[,2], means[,2])
colnames(DF3) <- c("EnvR", "PDmin", "PDmax", "PDmean")
head(DF3)
# plot min-max ranges
polygon(x=c(DF3$EnvR,rev(DF3$EnvR)),y=c(DF3$PDmax ,rev(DF3$PDmean)),col="grey",border=NA)
polygon(x=c(DF3$EnvR,rev(DF3$EnvR)),y=c(DF3$PDmin,rev(DF3$PDmean)),col="grey",border=NA)
# add average line (black)
meancolor <- "black"
lines(x=DF3$EnvR,y=DF3$PDmean,col=meancolor, lwd = 1, ljoin = 0)
 
   rug(quantile(Ama_F.M2[,f250m[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)
}

dev.off()


#CORRELATION MATRIX #----
##### CO-LINEARITY OFF PREDICTOR VARIABLES #
cordf <- Ama_F.M2[,c(4:11)]
# CORRELATION OF PREDICTOR VARIABLES  for correlation matrix below

cordf <- na.omit(cordf)
cordf<-cor(cordf)
head(round(cordf,2))

cor.mtest <- function(mat, ...) {
   mat <- as.matrix(mat)
   n <- ncol(mat)
   p.mat<- matrix(NA, n, n)
   diag(p.mat) <- 0
   for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
         tmp <- cor.test(mat[, i], mat[, j], ...)
         p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
   }
   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
   p.mat
}


# CORRELATION MATRIX # 
# matrix of the p-value of the correlation (for each species, can look at the corr between regional and national)
#require(corrplot)
emf(file = "Ama.corr.matrix.ABU.TGB.emf", emfPlus = FALSE, custom.lty = F) 
par(mar=c(4, 2, 1, 1)) # PAR = margins / change layout of plot grids
#par(mfrow=c(3,2))
p.mat <- cor.mtest(cordf)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corr.plt <- corrplot(cordf, method="color", col=col(200),    
                     type="upper", order="hclust", 
                     addCoef.col = "black", # Add coefficient of correlation
                     tl.col="black", tl.srt=45, #Text label color and rotation
                     # Combine with significance
                     p.mat = p.mat, sig.level = 0.05, insig = "blank", 
                     # hide correlation coefficient on the principal diagonal
                     diag=F)
dev.off()

# Importance of Environmental predictors
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
row.names(influence_mat)<-as.character(env_var_ord[,1])
preds_influences<-apply(influence_mat, 1, function(x) c(mean = mean(x), sd = sd(x)))   #Calculate mean and standard error of the relative influecne of each preds
write.csv(preds_influences, file="M2.preds.influences_AMA_ABU_TGB.csv")

# MODEL FIT METRICS (best to take RMSE over sd)----
# for Presence /rel abs mode use - Internal AUC as metric
#for Abundanced model, use Pearsons correlation (between  modelled data/model fit [m] and observed mean abundance [o]) = Internally Validated 
colnames(deviance_mat) <- c("Dev.Exp.int","Dev.Exp.eval","AUC.int","AUC.eval", "RMSE", "P.corr") # add RMSE column ( instead of deviacne explained for Abu, best do both and remove not needed as will be needed for PA but not Abu)
mean.model.fit <-apply(deviance_mat, 2, function(x) c(mean = mean(x), sd = sd(x)))
write.csv(mean.model.fit, file="mean.model.fit.metrics.M2.AMA.TGB.csv")

# calculate mean suitability and CV spatially
Ama.boot.mean <-apply(boot_mat,1,mean) # apply mean function, over the bootstrap matrix's rows. (1 = row, 2 = columns)
#uncertainty ( is the standard deviation between model runs)
Ama.boot.sd <-apply(boot_mat,1,sd)
Ama.boot.cv<-Ama.boot.sd/Ama.boot.mean # calculation of Coefficient of variation if needed

# export the map into X,Y,Z format
Ama.map.mean <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Ama.boot.mean) # Env Vars X & Y (column) and Z (mean prediction from model bootstrap matrix)
Ama.map.UC <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Ama.boot.sd) # uncertainty 
# and others
Ama.map.mean.pi.5 <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Ama.boot.pi.5)
Ama.map.mean.pi.95 <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Ama.boot.pi.95)
# convert to raster
BRT_Ama.mean <- rasterFromXYZ(data.frame(x =Ama.map.mean[,1], 
                                        y =Ama.map.mean[,2], 
                                        z =Ama.map.mean[,3]),
                             crs = crs("+init=epsg:9191")) 

BRT_Ama.UC <- rasterFromXYZ(data.frame(x =Ama.map.UC[,1], 
                                      y =Ama.map.UC[,2], 
                                      z =Ama.map.UC[,3]),
                           crs = crs("+init=epsg:9191")) 


plot(BRT_Ama.mean)
plot(BRT_Ama.UC)


# EXPORT AS TIFF - file for viewing in GIS----
writeRaster(BRT_Ama.mean,filename = "AMA_BRT_boots_mean_ABU_TGB.tif", overwrite=T, progress = "window")
writeRaster(BRT_Ama.UC,filename = "AMA_BRT_boots_mean_ABU_UC_TGB.tif", overwrite=T, progress = "window")

#-9---USE SIMPLIFY.GBM TO REFINE ENV PREDS (AND FOLLOW STEPS 1-8 ABOVE AGAIN BEFORE BOOSTRAPAMAG-------------
# we were able to run the GBM model, with just the presences and abscences.
# we split into 75%/25% training and eval data, which we will later bootstrap. 

#1st run
simp.sub2.run1 <- gbm.simplify(M2)
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
write.csv(simp.sub2.run1$drop.count, file="pred_simplify_output_run1.csv")

# will now re-run the above model, subsetting the removal list generated
# by simplify.
#the latter can be used as an argument to gbm.step e.g., gbm.step(data = data, gbm.x = simplify.object$pred.list[[4]]... 
#would implement a new analysis with the original predictor set, minus its four lowest contributing predictors

#2nd run
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
simp.sub2.run2 <- gbm.simplify(M2)
write.csv(simp.sub2.run2$drop.count , file="pred_simplify_output_run2.csv")
# as per output, only removing 2 lowest varaibles this time

#3rd run
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
simp.sub2.run3 <- gbm.simplify(M2)
write.csv(simp.sub2.run3$drop.count , file="pred_simplify_output.csv")

#-10---SCALING AND COMBINING OF ALL MODELS ####

#---HURDLE 1  = Pres / relative absence (mean from bootstrap) -----------------------------------------
#P/Ra
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/AMA/TGB")
H1_rast <- raster("AMA_BRT_boots_mean_Pres.Ra_TGB.tif") # load in pres/Abs model you created 
# Abundance 
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/ABU/AMA/TGB")
H1.2_rast <- raster("AMA_BRT_boots_mean_ABU_TGB.tif") # abundance model you just created from this script
# Hurdle creation 
H1 <- H1_rast * H1.2_rast # combine together to make hurdle model 
raster::scale(H1) # normalizing raster to be on a scale of 0-1
plot(H1)
# save as  
setwd("C:/Users/sw277/Desktop/R/PhD/R_working/ModelOutput_Oct_20/HUR")
writeRaster(H1,filename = "H1_AMA_HUR_TGB.tif", overwrite=T, progress = "window")


##----------------------------END---------------------####
