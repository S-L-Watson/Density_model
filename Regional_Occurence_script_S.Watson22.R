         
 #-----What:--------- BRT Analysis, Presence Pseudo Absence model for benthic invertebrates  -----------###
# Step 1  - Run the Presence Absence Model
# Step 2  - Run the Abundance model 
# Step 3  - Combine the outputs from each model, scale or normalize
 
#--- Supporting script for publication -- Watson, Stephanie Louise and Stephenson, Fabrice and Conrad, Pilditch A. and Lundquist, Carolyn, 
# Improving Predictions of Coastal Benthic Invertebrate Occurrence and Density Using a Multi-Scalar Approach.
# Available at SSRN: https://ssrn.com/abstract=4142118 or http://dx.doi.org/10.2139/ssrn.4142118
                

 
                                          ################ ---- START ----######################

#  packages
require(raster);require(dismo);require(gbm); require(sp); require(pROC);require(FinPec);require(ggplot2);require(dplyr)
require(rgdal); require(corrplot); require(devEMF);require(rcompanion); require(DescTools); require(tools);require(BBmisc)             
#----LOADING OF DATA STARTS----
#---1---LOAD species data clipped to mask ####
setwd("D:/R/PhD/Invert_data/Oct_20/presence")
Pec <- read.csv("PEC_TGB_pres.csv")
Pec$pa <- 1 # add new column for presence absence 
Pec <- data.frame(Pec$POINT_X, Pec$POINT_Y, Pec$pa) # bind relevant columns
colnames(Pec) <- c("X", "Y", "pres") 
head(Pec) # a data frame with presence points/rows, with lat, long 
                
#---2---LOAD ENVIRO PREDS ########## 
#TGB Wide Env Preds
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/PEC/TGB/Env_preds")
f250m <- list.files(getwd()) 
imp.var <-c("X","Y","FID",c(file_path_sans_ext(f250m)))
ras250m <- lapply(f250m, raster)
f250m <- file_path_sans_ext(f250m)
plot(ras250m[[1]]) # check through layers make sure same extent and projection
env250m_rs <- raster::stack(ras250m) 

# #--3--- CREATE SPDF: Turn original raster stack and spatial points data into a data frame for model prediction#####
env250m_spdf <- rasterToPoints(env250m_rs, spatial = T) # turns data into a large spatial points dataframe
# Extract the coordinates, for the spatial points df, and bind it with the env data from the raster stack,
# and save as a new dataframe
env250m_df <- data.frame(env250m_spdf@data, X =env250m_spdf@coords[,1],Y = env250m_spdf@coords[,2])
head(env250m_df)

# #---4--SET UNIQUE FID, EXTRACT X & Y( row and column numbers, not lat and long) data from raster stack ####
env250m_spdf@coords # a command to see the x and y from the raster data
env250m_df$FID <- cellFromXY(env250m_rs, env250m_df[,c("X", "Y")])
env250m_df <- env250m_df[,c(35:36,37,1:34)] # re-shuffle the order of columns, to put x and y at front
str(env250m_df) # check everything is named corectly
# --save for use next time to speed things up
EEZ.250m.data_extent.env.vars.train <- env250m_df
save(EEZ.250m.data_extent.env.vars.train, file = "EEZ.250m.data_extent.env.vars.train.Rdata")


#---4b---LOAD ENV_VAR_DF ONCE CREATED - from steps 1-4 above ----
setwd("D:/R/PhD/Environ_preds/Environ_preds/250m_standardised")
load("TGB.250m.all.env.vars.train.Rdata") # load all env variables DF

#5--- EXTRACT ENV.VAR INFO where presences occur ----
Pec$FID <- cellFromXY(env250m_rs, Pec[,c("X", "Y")])# from the object (raster stack), return the raster cells (row +cols), in which
#the coordinates from Pec (in the form of X and Y) are found within the raster, and place them into newly created  FID column 
Pec_F <- cbind(Pec ,raster::extract(env250m_rs, Pec[,c("X", "Y")])) # extract all env data where Pec occurs in env var raster stack
 
#6 AGGREGATE OCCURENCES BY FID TO 250M / ADD EQUAL WEIGHTING
#so multiple samples within 1 cell are a proportion of 1. This will provide EQUAL WEIGHTING for when we generate absences.
Pec_F <- Pec_F %>% # using the pipe operator as part of the dplyr package converts data to a tible, so use as.data.frame to convert back below
group_by(FID) %>%
summarize(wg = (1/sum(pres))) %>% # make sure that multiple points are calculated as a proportion for that raster cell /geographic sample area.
select(wg = wg, FID) %>%
left_join(Pec_F,wg,by= "FID")
# shuffle order so that pres, wg, and FID are last columns
Pec_F <- Pec_F[,c(3:4,6:12,5,1:2)]
# remove abundance as not needed in this model
Pec_F.M1 <- as.data.frame(Pec_F) # final, for model 1 (presence only)

#6b SAVE/LOAD FINAL SPECIES DATA---- 
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/PEC/TGB")
save(Pec_F.M1, file = "Pec_F.M1.Rdata")
rm(Pec_F.M1)

# Important Env vars, extracted from step 4 above (already refined from gbm simplify)
# define the important vars per species and model type then subset
imp.var <-c("X","Y","FID",c(f250m))  # create character index from env var files, and remove file ext eg .tif
TGB.250m.all.env.vars.train <- TGB.250m.all.env.vars.train[, imp.var]

#--7----TARGET GROUP BACKGROUND ABSENCES SELECTION----
#the UD (Utilization Distribution for the species was calculated, and absences were clipped to those boundaries, without the target species. 
#bottom 5% of lowest density records removed leaving the top 95% home range 
#load home range kernel UD background absences
setwd("D:/R/PhD/Invert_data/Oct_20/absence") # set directory 
abs.95 <- read.csv("PEC_abs_TGB.csv", header = T) # load 95% home range absences
abs.95$pa <- 0 # add new column for presence absence 
backg.abs.95 <- data.frame(abs.95$X, abs.95$Y, abs.95$pa) # bind relevant columns
colnames(backg.abs.95) <- c("X", "Y", "pres") # 
head(backg.abs.95) 
#Extract corresponding FID into this dataset, as we did for presences, so we can have common value for comparing
backg.abs.95$FID <- cellFromXY(env250m_rs, backg.abs.95[,c("X", "Y")])

# Extract all env data where absences occur (#dataframe with 11 variables)
backg.abs.95 <- cbind(backg.abs.95 ,raster::extract(env250m_rs, backg.abs.95[,c("X", "Y")])) 
# group by FID
backg.abs.95 <- as.data.frame(backg.abs.95 %>% group_by(FID))# using the pipe operator as part of the dplyr package converts data to a tible, so use as.data.frame to convert back below
# shuffle order so that abundance,pres and FID are last columns
backg.abs.95 <- backg.abs.95[,c(1:2,5:11,3:4)]
backg.abs.95 <- na.omit(backg.abs.95)
rm(env250m_rs)

#--MODEL OBJECTS(-setup objects to save information from model runs to compare metrics with diff ENV Preds)####
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/PEC/TGB")
n.boot <- 100
length.map <- nrow(TGB.250m.all.env.vars.train)  # previously env250m_df (THIS CAN BE YOUR SMALLER STUDY AREA)
boot_mat <- array(0, c(length.map,n.boot)) # matrix where rows will store predictions for each bootstrap (columns)
deviance_mat <- array(0, c(n.boot,6)) # for saving model fit metrics ( length of bootstaps with 4 columns)
influence_mat <- array(0,c(7,n.boot)) #for saving information on predictors (15 rows, 100 columns)

#Setup for partial dependance plots
PredMins <- apply(na.omit(Pec_F.M1)[,f250m],2,min)
PredMaxs <- apply(na.omit(Pec_F.M1)[,f250m],2,max)

EnvRanges <- as.data.frame(array(0,c(100,length(f250m)))) # store a DF with a sequence of env ranges from min to max, from all raw values in imp.vars
names(EnvRanges) <- f250m
for (i in c(1:length(f250m))) {
   EnvRanges[,i] <- seq(PredMins[f250m[i]],PredMaxs[f250m[i]],length = 100)
}
# for each model run (bootstrap) store the results of environmental ranges, through the range
PD <- array(0,c(length(EnvRanges[[1]]),length(EnvRanges),n.boot)) # #3D vector ( an array) with columns given same names as imp.vars
dimnames(PD)[[2]] <- f250m
dimnames(PD)[[3]] <- paste('Rep_',seq(1,n.boot),sep="") # the 3rd dimension is names the rep number of the bootstrap

#---8-CREATE RANDOMLY SAMPLED PRESENCE /ABSCENSE DATA POINTS ---####
counter = 0
for (i in 1:n.boot){ # bootstrap loop for M1 ( presence / Absence) 
               
#(1) create randomised sample of presence data  (#presence selection, before absence)
ind_pres <- sample(seq_len(nrow(Pec_F.M1)), size = floor(0.75 * nrow(Pec_F.M1))) # create an index vector - 75% of data
# using the above index, subset from the original dataset to create training and evaluation presences 
#PRESENCES
pres_train <-Pec_F.M1[ind_pres, ] # 75%
pres_eval <-Pec_F.M1[-ind_pres, ] # 25%

#(2) Generate absence pool

absence.pool <- backg.abs.95[!backg.abs.95$FID %in% Pec_F.M1$FID, ] # creates pool of, of our known absences minus where a presence FID is found. ( all the FID's, that are not in Pec)

#(3) create randomised sample of absence data - 
ind_abs_train <- sort(sample(seq_len(nrow(absence.pool)),size = nrow(pres_train)*2, replace = FALSE)) # create random training sample from the absences pool 
ind_abs_eval <- sort(sample(seq_len(nrow(absence.pool)),size = nrow(pres_eval)*2, replace = FALSE))# create random eval sample from the absences pool 
#ABSENCES
abs_train <-absence.pool[ind_abs_train, ] #just the randomly sampled abscenes (same length as presences) #498
abs_eval <-absence.pool[ind_abs_eval, ]  # just the randomly sample eval rows, accross all columns of data  #166  

#(4)
#ABSCENCES TRAINING SUSBSET
abs.train.wg <- abs_train
abs.train.wg$wg <- 1
#ABSCENCES EVALUATION SUSBSET 
abs.eval.wg <- abs_eval
abs.eval.wg$wg <- 1

#(5) LOOP FOR EQUAL WEIGHTING - apply loop so that weighting of abscences from sub-sample will equal that of presences
#Temp = abs data (a randomised, then weighted, large pool of random absences selected from within kernel density range of presences);
# thresh = your sum of weightings from your presences
thresh.tr <- sum(pres_train$wg)
temp.tr <- abs.train.wg[1,] # first row of abscenes data 

# TRAINING DATA 
for (j in 2:nrow(abs.train.wg)){
 temp.tr <- rbind(temp.tr, abs.train.wg[j,]) # new samples each time through each iteration 
if(sum(temp.tr$wg) >= thresh.tr){break}
}
Pec.PresAbs.tr <- as.data.frame(rbind(pres_train,temp.tr)) # final presence/absence data (training group 75%) 
 
# EVALUATION DATA  
thresh.ev <- sum(pres_eval$wg)
temp.ev <- abs.eval.wg[1,]

for (j in 2:nrow(abs.eval.wg)){
temp.ev <- rbind(temp.ev, abs.eval.wg[j,])
if(sum(temp.ev$wg) >= thresh.ev){break}
}

Pec.PresAbs.ev <- as.data.frame(rbind(pres_eval,temp.ev)) # final pres / abs data (evaluation group 25%) 

#---- 9-- M1 TRAINING MODEL for presence and PSEUDO absence / export to raster #####
# once you run through the model once, use gbm simplify to reduce the env preds, may need to run more than once. 

# #gbm.step - Fits a gbm model to one or more response variables, using cross-validation to estimate the optimal number of trees.
NoTrees <- 0
lr <- 0.001
while(NoTrees < 500) {
M1<-gbm.step(data=Pec.PresAbs.tr,
           gbm.x =f250m,      
           gbm.y = 10, # column number of presence/ relative absence data (aka the response variable)
           site.weights = Pec.PresAbs.tr$wg, # weighting column
           family = "bernoulli", tree.complexity = 2, #(aka Interaction Depth) - play around with settings to obtain 1500-3000 trees
           learning.rate = lr, bag.fraction = 0.7, n.folds=5,  # bag fraction between 0.5-0.75 is good   #n.folds = upto 10, take it down to 5 if not converging
           max.trees = 10000, step.size = 25, plot.main = T, tolerance.method = "fixed", # change 'plot main' to T when doing initial runs
           tolerance = 0.01, verbose = T) # switch this to T to visualize AUC /tree number for first run before looping ( then switch to F)
   NoTrees <- M1$n.trees
   if(is.null(NoTrees)){
      NoTrees <- 1
      lr <- lr *0.8
   }else
   {next}
}


#In gbm, stochasticity is controlled through a 'bag fraction' that specifies the proportion
#of data to be selected at each step. The default bag fraction is 0·5, meaning that, at each iteration, 
# 50% of the data are drawn at random, without replacement, from the full training set
#partial dependence plots
for (j in 1:length(f250m)){
   grid <- gbm::plot.gbm(M1, i.var = c(paste(f250m[[j]])), return.grid = T)
   PD[,f250m[[j]],i] <- loess(grid$y~EnvRanges[,f250m[[j]]])$y # LOESS - local Polynomial Regression Fitting 
}  #allows for the extraction of best model fit, within the whole range of env preds used

#--internal model fit metrics
int.null.deviance <- M1$self.statistics$mean.null # trying to explain the deviance. If it is 0 it cannot explain any of the deviance
int.residual.deviance <- M1$cv.statistics$deviance.mean 
deviance_mat[i,1] <- (int.null.deviance-int.residual.deviance)/int.null.deviance # internal Deviance explained

#model fit comparison with withheld evaluation data
pred <- predict.gbm(M1, Pec.PresAbs.ev, n.trees = M1$gbm.call$best.trees, type = "response") 
ext.residual.deviance <- calc.deviance(Pec.PresAbs.ev$pres, pred, family = "bernoulli" ,calc.mean=T) 
ext.null.deviance <- calc.deviance(Pec.PresAbs.ev$pres,family = "bernoulli", rep(mean(Pec.PresAbs.ev$pres),nrow(Pec.PresAbs.ev)), calc.mean=T) 
deviance_mat[i,2]<-(ext.null.deviance - ext.residual.deviance)/ext.null.deviance

#-- RMSE  & P.CORR -calculate the modelled and observed values for RMSE 
m = M1$fit # mean model fit (predicted from model)
o = Pec.PresAbs.tr$pres # observed from input data frame
plot(o,m) # plot observed vs modeled obs  - blank out for bootstrap
deviance_mat[i,5] <- RMSE(m,o)   #True
deviance_mat[i,6] <- cor(o,m) # pearsons correlation between model fit and presence 
#Correlation coefficient between X and Y.if it is between >0.5 then it is a strong correlation

# external AUC
deviance_mat[i,4] <- roc(Pec.PresAbs.ev$pres, pred)$auc
#internal AUC
deviance_mat[i,3] <- M1$cv.statistics$discrimination.mean

# Env pred contribution
M1_contrib <- as.data.frame(M1$contributions) # M1 = your BRT model object
env_var_ord <- M1_contrib[order(M1_contrib$var),]
influence_mat[,i]<-env_var_ord[,2]
#env_var_ord

# predict spatially to map/study area (change out to smaller prediction area, than size of training tiffs if needed)
boot_mat[,i] <- round(predict.gbm(M1, TGB.250m.all.env.vars.train, n.trees = M1$gbm.call$best.trees, family = "bernoulli", type = "response"), digits = 2)
                              
counter = counter + 1
print(counter)
}
#------------------------------------END OF MODEL LOOP #-------------------------------------------####
summary(M1)

#-10---CREATE OUTPUTS & save data-####
# save created data
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/PEC/TGB")
save(boot_mat, file="M1_boot_mat.Rdata")
save(deviance_mat, file="M1_dev_mat.Rdata")
save(influence_mat, file="M1_inf_mat.Rdata")

#-- PARTIAL DEPENDENCE PLOTS #----
#require(BBmisc)  #= normalize
PD1 <- PD
PD1 <- normalize(PD, method = "range", range = c(0, 1), margin = 1) # normalize the Y axis on PD plots to be on same scale 
# plot PDs + 95 PI and save to file
emf(file = "PD_Pec.mean.TGB.emf", emfPlus = FALSE) 
par(mar=c(4, 2, 1, 1)) 
par(mfrow=c(3,3))
#explanation: Go through the imp vars, and for each value row, find the corresponding value in the environmental ranges, and plot 
# the mean predictor response for that value. 
for (i in c(1:length(f250m))) { # or imp.var
        plot(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,mean), col = 1,type= 'l', 
        xlab = paste(f250m[i]), 
        ylab = '',
        ylim = c(min(PD1[,,]), max(PD1[,,]))) 
 
# SOLUTION fOR RANGE IN SHADED POLYGON# in a DF
 maxs <- data.frame(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,max))
 mins <- data.frame(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,min))
 means <- data.frame(EnvRanges[,i],apply(PD1[,i, ,drop=F],1,mean))
 DF3 <- data.frame(mins[,1], mins[,2],maxs[,2], means[,2])
colnames(DF3) <- c("EnvR", "PDmin", "PDmax", "PDmean")

# plot min-max ranges
polygon(x=c(DF3$EnvR,rev(DF3$EnvR)),y=c(DF3$PDmax ,rev(DF3$PDmean)),col="grey",border=NA)
polygon(x=c(DF3$EnvR,rev(DF3$EnvR)),y=c(DF3$PDmin,rev(DF3$PDmean)),col="grey",border=NA)
# add average line (black)
meancolor <- "black"
lines(x=DF3$EnvR,y=DF3$PDmean,col=meancolor, lwd = 1, ljoin = 0)
   rug(quantile(Pec_F.M1[,f250m[i]],seq(0,1,0.1), na.rm = T), ticksize = 0.05, side = 1, lwd = 0.75)

}

dev.off()

#CORRELATION MATRIX #----
##### CO-LINEARITY OFF PREDICTOR VARIABLES #
cordf <- Pec_F.M1[,f250m] 
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
emf(file = "Pec.corr.matrix.TGB.emf", emfPlus = FALSE, custom.lty = F) 
par(mar=c(4, 2, 1, 1))
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
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/PEC/TGB")
row.names(influence_mat)<-as.character(env_var_ord[,1])
preds_influences<-apply(influence_mat, 1, function(x) c(mean = mean(x), sd = sd(x)))   #calculate mean and standard error of the relative influence of each preds
write.csv(preds_influences, file="M1.preds.influences_Pec_Pres.Ra_TGB.csv")

# MODEL FIT METRICS 
# for Presence /rel abs mode use - Internal AUC as metric
#for Abundance model, use Pearson correlation (between  modelled data/model fit [m] and observed mean abundance [o]) = Internally Validated 
colnames(deviance_mat) <- c("Dev.Exp.int","Dev.Exp.eval","AUC.int","AUC.eval", "RMSE", "P.corr") 
mean.model.fit <-apply(deviance_mat, 2, function(x) c(mean = mean(x), sd = sd(x)))
write.csv(mean.model.fit, file="mean.model.fit.metrics.M1.Pec.TGB.csv")

# calculate mean suitability and CV spatially
# want to create a map for each of these ( p.5%, mean, and p.95% PREDICTION INTERVALS)
Pec.boot.mean <-apply(boot_mat,1,mean) # apply mean function, over the bootstrap matrix's rows. (1 = row, 2 = columns)
#uncertainty ( is the standard deviation between model runs)
Pec.boot.sd <-apply(boot_mat,1,sd)
Pec.boot.cv<-Pec.boot.sd/Pec.boot.mean # calculation of Coefficient of variation if needed

# export the map into X,Y,Z format
Pec.map.mean <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Pec.boot.mean) # Env Vars X & Y (column) and Z (mean prediction from model bootstrap matrix)
Pec.map.UC <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Pec.boot.sd) # uncertainty 
# and others
Pec.map.mean.pi.5 <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Pec.boot.pi.5)
Pec.map.mean.pi.95 <- cbind(TGB.250m.all.env.vars.train[, c("X", "Y")],Pec.boot.pi.95)
# convert to raster
BRT_Pec.mean <- rasterFromXYZ(data.frame(x =Pec.map.mean[,1], 
                                        y =Pec.map.mean[,2], 
                                        z =Pec.map.mean[,3]),
                             crs = crs("+init=epsg:9191")) 

BRT_Pec.UC <- rasterFromXYZ(data.frame(x =Pec.map.UC[,1], 
                                      y =Pec.map.UC[,2], 
                                      z =Pec.map.UC[,3]),
                           crs = crs("+init=epsg:9191")) 


plot(BRT_Pec.mean)
plot(BRT_Pec.UC)


# EXPORT AS TIFF - file for viewing in GIS----
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/PRES_RA/PEC/TGB")
writeRaster(BRT_Pec.mean,filename = "PEC_BRT_boots_mean_Pres.Ra_TGB.tif", overwrite=T, progress = "window")
writeRaster(BRT_Pec.UC,filename = "PEC_BRT_boots_mean_UC_Pres.Ra_TGB.tif", overwrite=T, progress = "window")

#-11---USE SIMPLIFY.GBM TO REFINE ENV PREDS (AND FOLLOW STEPS 1-10 ABOVE AGAIN BEFORE BOOSTRAPPING-

#1st run
simp.sub2.run1 <- gbm.simplify(M1)
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/ABU/Pec/TGB/1st_run_all_Env_preds")
write.csv(simp.sub2.run1$drop.count, file="pred_simplify_output_run1.csv")

# will now re-run the above model, subsetting the removal list generated by simplify.
#the latter can be used as an argument to gbm.step e.g., gbm.step(data = data, gbm.x = simplify.object$pred.list[[4]]... 
#would implement a new analysis with the original predictor set, minus its four lowest contributing predictors

#2nd run
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/ABU/Pec/TGB/1st_run_all_Env_preds")
simp.sub2.run2 <- gbm.simplify(M1)
write.csv(simp.sub2.run2$drop.count , file="pred_simplify_output_run2.csv")
# as per output, only removing 2 lowest varaibles this time

#3rd run
setwd("D:/R/PhD/R_working/ModelOutput_Oct_20/ABU/Pec/TGB/1st_run_all_Env_preds")
simp.sub2.run3 <- gbm.simplify(M1)
write.csv(simp.sub2.run3$drop.count , file="pred_simplify_output_run3.csv")


##----------------------------END---------------------####

# for Hurdle model creation, now run the "Regional Abundance " script
# and combine final outputs from each to make hurdle tiff. 
