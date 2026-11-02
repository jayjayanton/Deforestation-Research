library(caret)
library(corrplot)
library(car)
library(randomForest)
library(PresenceAbsence)
library(stargazer)
library(randomForest)
library(rpart)
library(rpart.plot)
library(ROCR)
library(ecospat)
library(dismo)
library(dplyr)
library(spdep)
library(sandwich)
library(factoextra)
library(readxl) # Reads excel files
library(spdep)

#####
#DATA
#####

setwd('C:/Users/jayja/Documents/Statistical analysis')
# read in data: note this already has redundant variables and NAs removed, and distance units converted to km
# this dataset also has an 'id cell' column, that is equivalent to an objectID
env_dat <- read_excel("Amazonia_2004_2017_03_24_sel_coord.xlsx")
attach (env_dat)

env_dat <- replace(env_dat,grepl( "_dist" , names( env_dat )),
                          env_dat[ , grepl( "_dist" , names( env_dat ))]*.001)

str(env_dat)
head(env_dat)
summary(env_dat)
names(env_dat)
env_dat$DEFORESTATION<- as.factor(env_dat$DEFORESTATION)

is.na(env_dat) #checking for NaNs
sum(is.na(env_dat))
### Replace NA values with mean value in the columns
for(i in 1:ncol(env_dat)){
  env_dat[is.na(env_dat[,i]), i] <- mean(env_dat[,i], na.rm = TRUE)
}


env_dat<-na.omit(env_dat)
sum(is.na(env_dat))
#####
#GLMS 
#####

which( colnames(env_dat)=="DEFORESTATION" ) # Identify Column 60
which( colnames(env_dat)=="RIVERS_SUM_METERS" )

# FULL MODEL
f <- paste(names(env_dat)[60], "~", paste(names(env_dat)[c(-1:-5, -9,-10, -12, -14:-17, -23, -24, -26, -29, -30:-35,-38, 
                                                           -40:-42, -44:-49, -50, -53:-59, -67:-68, -71,-72)], collapse="+")) 

defor.m <- glm(f, data=env_dat, family="binomial")
summary(defor.m)
defor.m.cov <- vcovHC(defor.m, type = "HC3");defor.m.rob.std.err <- sqrt(diag(defor.m.cov)) # calculate Huber White clustered robust st errors

# STEPWISE MODEL
defor.step <- step(defor.m, direction = c("both")) # stepwise drops 3 vars: peat_lands_dist, roads_2nd_dist, tp_rubber_sqkm 
defor.step.cov <- vcovHC(defor.step, type = "HC3"); defor.step.rob.std.err <- sqrt(diag(defor.step.cov)) 

# Variables dropped from model: 
which( colnames(env_dat)=="FB_POP_SUM_GRIDCODE") #[6]
which( colnames(env_dat)=="water_lines_dist" ) #[13]
which( colnames(env_dat)=="agrar_reform_settlements_dist" )#[22]
which( colnames(env_dat)=="oil_gas_concessions_dist") #[25]

# NAIVE MODEL
defor.naive <- glm(DEFORESTATION~ ELEVATION_SRTM90_AMAZONIA, data=env_dat, family="binomial")
defor.naive.cov <- vcovHC(defor.naive, type = "HC3");defor.naive.rob.std.err <- sqrt(diag(defor.naive.cov)) 

#compare these side-by-side:
stargazer(defor.m, defor.step, defor.naive, type = "html", no.space= TRUE, title = "Alternative Models",
          column.labels = c('Full', 'Stepwise', 'Naive', 'RAC'), 
          se=list(defor.m.rob.std.err, defor.step.rob.std.err, defor.naive.rob.std.err),
          digits = 3, model.numbers=FALSE, df=FALSE, out="Primary_models_Brazil.htm")

# RAC MODEL (following Crase et al. 2012)
xy <- cbind(env_dat$CENTROID_X, env_dat$CENTROID_Y) #extract xy co-ordinates from data
ac<-autocov_dist(defor.m$residuals ,xy,type="inverse",style="B",nbs = 10,longlat=TRUE,zero.policy=TRUE) #derive the residual autocovariate term 
env_dat_ac <-cbind(env_dat,ac)#add residual autocovariate term to dataset
g <- paste(names(env_dat_ac)[60], "~", paste(names(env_dat_ac)[c(-1:-6, -9,-10, -12:-17, -22:-26, -29, -30:-35,-38, 
                                                                 -40:-42, -44:-49, -50, -53:-59, -67:-68, -71,-72)], collapse="+")) #uses same covariates as in the stepwise model
defor.RAC <- glm(g, data=env_dat_ac, family="binomial")
summary(defor.RAC)
defor.RAC.cov <- vcovHC(defor.RAC, type = "HC3");defor.RAC.rob.std.err <- sqrt(diag(defor.RAC.cov)) 

#compare these side-by-side:
stargazer(defor.m, defor.step, defor.naive, defor.RAC, type = "html", no.space= TRUE, title = "Alternative Models",
          column.labels = c('Full', 'Stepwise', 'Naive', 'RAC'), 
          se=list(defor.m.rob.std.err, defor.step.rob.std.err, defor.naive.rob.std.err, defor.RAC.rob.std.err),
          digits = 3, model.numbers=FALSE, df=FALSE, out="Alt_models_Brazil.htm")


######
# EVAL
######

# generate a test data set of structure: ID, observation, prediction (where ID is the row number of the dataset) 
test.defor <- data.frame(cbind(as.numeric(rownames(env_dat)), env_dat$DEFORESTATION, defor.m$fitted.values, defor.step$fitted.values, defor.naive$fitted.values, defor.RAC$fitted.values))
names(test.defor)<-c("ID","Observed","Predicted.Full", "Predicted.Step", "Predicted.Naive", "Predicted.RAC")
test.defor$Observed[test.defor$Observed==1] <- 0; test.defor$Observed[test.defor$Observed==2] <- 1 # not sure why, but these are pulled from factors to numbers and assigned values of 1 and 2 rather than 0 and 1
glm.cut <- as.numeric(optimal.thresholds(test.defor, opt.methods=10, req.sens = 0.90)[2]) ## Select a threshold where sensitivity is at least 90%
defor.acc <- presence.absence.accuracy(test.defor, threshold=glm.cut, st.dev=F)
defor.acc$TSS <- defor.acc$sensitivity+defor.acc$specificity-1  # add true skills statistic to our accuracy metrics
defor.acc
auc.roc.plot(test.defor, main = "Brazil Model ROC Plots", color = c('Red', 'Blue', 'Green', 'Yellow')) # plot ROC curve

###########
# RANDOM CV 
###########

#create training and testing sets:
Train <- createDataPartition(env_dat$DEFORESTATION, p=0.7, list=FALSE)
training <- env_dat[ Train, ]; testing <- env_dat[ -Train, ]

# create confusion matrix, estimate Auc and plot Roc
ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, savePredictions = TRUE, classProbs = TRUE, summaryFunction = twoClassSummary)
levels(training$DEFORESTATION) <- c("zero", "one"); levels(testing$DEFORESTATION) <- c("zero", "one")
mod_fit <- train(DEFORESTATION ~ ELEVATION_SRTM90_AMAZONIA + SLOPE_SRTM90_AMAZONIA
                 + water_areas_dist + Comun_Indigenas_dist + Territorios_Indigenas_dist
                 + fb_population_counts_2018_dist + afrobrazilian_settlements_dist + reg_rural_prop_dist 
                 + tree_plant_dist + AFROBRA_SETTLEM_SUM_HA + AFROBRA_SETTLEM_SUM_familia
                 + MINING_SUM_HA + TERRIT_INDIGENAS_SUM_HA +TERRIT_INDIGEN_HABITAN_COUNT
                 + c2004_mining_dist + c2004_cropland_dist + c2004_pasture_dist
                 + c2004_urbaninfra_dist + roads_primary_secondary_dist
                 + roads_tertiary_local_dist + roads_prim_sec_km + roads_tert_loc_km
                 , data=training, method="glm", family="binomial", trControl = ctrl, metric="ROC")
summary(mod_fit) 
mean(mod_fit$resample$ROC); sd(mod_fit$resample$ROC) # mean and st dev of AUC across all folds
mean(mod_fit$resample$Sens); sd(mod_fit$resample$Sens) # mean and st dev of sensitivity across all folds
mean(mod_fit$resample$Spec); sd(mod_fit$resample$Spec) # mean and st dev of specificity across all folds

pred <- predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$DEFORESTATION)


############
#spatial CV 
############

#  Generate spatial clusters to be used in spatial K-fold Cross-Validation (Spatial K-fold CV)
mdist <- dist(env_dat[c("CENTROID_X", "CENTROID_Y")])
hc <- hclust(mdist, method='complete')
d=150*1000  # the maximum distance between pixels within clusters (150 m * 1000 m = 150 km)                   
env_dat$Clust_150km = cutree(hc, k=10, h=d) 
plot(env_dat$CENTROID_X, env_dat$CENTROID_Y, col=as.factor(env_dat$Clust_150km), xlab="Longitude", ylab="Latitude")
legend("topright", legend=paste('group', 1:10), pch=16, col=(1:10), cex=0.8) 
#SPATIAL K-fold CV
h <- paste(names(env_dat)[60], "~", paste(names(env_dat)[c(-1:-6, -9,-10, -12:-17, -22:-26, -29, -30:-35,-38, 
                                                           -40:-42, -44:-49, -50, -53:-59, -67:-68, -71,-72)], collapse="+")) #uses same covariates as in the stepwise model
out.glm <- data.frame()
for (j in 1:10)
{
  id = which(env_dat$Clust_150km == j)
  training_data= env_dat[-id,]
  validation_data = env_dat[id,]  
  
  model <- glm(h, data=training_data, family="binomial")

  # prediction
  pred_obj = predict(model, newdata = validation_data, type="response")
  PRED.TEST = data.frame(id_cell = validation_data$OBJECTID,
                         deforestation = validation_data$DEFORESTATION,
                         pred = pred_obj)
  PRED.TEST$deforestation<- as.numeric(PRED.TEST$deforestation)
  PRED.TEST$deforestation[PRED.TEST$deforestation==1] <- 0; PRED.TEST$deforestation[PRED.TEST$deforestation==2] <- 1
  #performance
  thresh <- (optimal.thresholds(PRED.TEST, opt.methods="MaxSens+Spec")[2]) # Throws NA Value
  acc <- presence.absence.accuracy(PRED.TEST, threshold = 0.69, find.auc = TRUE, st.dev=TRUE)
  acc$TSS <- acc$sensitivity+acc$specificity-1 ## $ operator is invalid for atomic vectors
  acc$model <- j
  out.glm <- rbind(out.glm, acc)
  
  # store results
  if(j == 1) { PRED.TEST_tot = PRED.TEST }
  if(j > 1)  { PRED.TEST_tot = rbind(PRED.TEST_tot, PRED.TEST) }
  
} # end of spatial K-fold CV loop
out.glm
#PRED.TEST_tot  






