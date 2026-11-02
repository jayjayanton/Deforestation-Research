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

library(ggpubr)
library(factoextra)

#####
#DATA
#####

setwd('C:/Users/jayja/Documents/Statistical analysis')
# read in data: note this already has redundant variables and NAs removed, and distance units converted to km
# this dataset also has an 'id cell' column, that is equivalent to an objectID
env_dat_Col <- read.csv("Colombia_Analysis_03_07_2021_No_IDs.csv")
attach (env_dat_Col)
# 
env_dat_Col <- replace(env_dat_Col,grepl( "_Dist" , names( env_dat_Col )),
                   env_dat_Col[ , grepl( "_Dist" , names( env_dat_Col ))]*.001)
# 
# str(env_dat_Col)
# head(env_dat_Col)
# summary(env_dat_Col)
# names(env_dat_Col)
env_dat_Col$Deforestation<- as.factor(env_dat_Col$Deforestation)
#env_dat_Col$Deforestation[env_dat_Col$Deforestation==1] <- 0; env_dat_Col$Deforestation[env_dat_Col$Deforestation==2] <- 1

env_dat_Col[is.na(env_dat_Col)] <- 0
sum(is.na(env_dat_Col))

is.na(env_dat_Col) #checking for NaNs
sum(is.na(env_dat_Col))
### Replace NA values with mean value in the columns
# 
# env_dat_Col<-na.omit(env_dat_Col)
# sum(is.na(env_dat_Col))

#####
#GLMS 
#####

 which( colnames(env_dat_Col)=="Deforestation" ) # Identify Column [7]
# 
 env_dat_Col <- env_dat_Col[c(-2:-6, -8)]
 names(env_dat_Col)


# FULL MODEL
f_Col <- paste(names(env_dat_Col)[2], "~", paste(names(env_dat_Col)[-c(1,2,5,6:15,17,19,20,26,28:32,
                                                                       35:38,40,41,43:45, 
                                                                       55,57:59,64:71,73:75)], collapse="+")) 
#Irrigated Cropland Dist Also Taken Out
defor_Col.m <- glm(f_Col, data=env_dat_Col, family="binomial")
summary(defor_Col.m)
vif(defor_Col.m)
defor_Col.m.cov <- vcovHC(defor_Col.m, type = "HC3");defor_Col.m.rob.std.err <- sqrt(diag(defor_Col.m.cov)) # calculate Huber White clustered robust st errors

# STEPWISE MODEL
Col_defor.step <- step(defor_Col.m, direction = c("both")) #  
Col_defor.step.cov <- vcovHC(Col_defor.step, type = "HC3"); Col_defor.step.rob.std.err <- sqrt(diag(Col_defor.step.cov)) 

# Variables dropped from model: 

which( colnames(env_dat_Col)=="Elevation" ) #[4]
which( colnames(env_dat_Col)=="Universal_Mills_Dist" )#[16]
which( colnames(env_dat_Col)=="PAs_Dist" ) #[34]
which( colnames(env_dat_Col)=="Ld_Grassland_Dist" ) #[42]
which( colnames(env_dat_Col)=="Ld_Cropland_Grassland_Mosaic_Dist") #46
which( colnames(env_dat_Col)=="Rd_Prim_Sec_Km") #[48]
which( colnames(env_dat_Col)=="Riv_Km") #[50]
which( colnames(env_dat_Col)=="Tr_Oil_Palm_SqKm") #[56]
which( colnames(env_dat_Col)=="Ld_Cropland_Grassland_Mosaic_SqKm") #[61]
which( colnames(env_dat_Col)=="Ld_Cropland_Woodland_Mosaic_SqKm") #[62]
which( colnames(env_dat_Col)=="Ld_Grassland_SqKm") #[63]
which( colnames(env_dat_Col)=="PAs_Dist") #[25]

# NAIVE MODEL
Col_defor.naive <- glm(Deforestation~ Elevation, data=env_dat_Col, family="binomial")
Col_defor.naive.cov <- vcovHC(Col_defor.naive, type = "HC3");Col_defor.naive.rob.std.err <- sqrt(diag(Col_defor.naive.cov)) 

#compare these side-by-side:
stargazer(defor_Col.m, Col_defor.step, Col_defor.naive, type = "html", no.space= TRUE, title = "Alternative Models",
          column.labels = c('Full', 'Stepwise', 'Naive', 'RAC'), 
          se=list(defor_Col.m.rob.std.err, Col_defor.step.rob.std.err, Col_defor.naive.rob.std.err),
          digits = 5, model.numbers=FALSE, df=FALSE, out="Primary_models_Colombia_V5.html")

# RAC MODEL (following Crase et al. 2012)
Col_xy <- cbind(env_dat_Col$Lattitude, env_dat_Col$Longitude) #extract Col_xy co-ordinates from data
Col_ac<-autocov_dist(defor_Col.m$residuals ,Col_xy,type="inverse",style="B",nbs = 10,longlat=TRUE,zero.policy=TRUE) #derive the residual autocovariate term 
env_dat_Col_ac <-cbind(env_dat_Col,Col_ac)#add residual autocovariate term to dataset
g_col <-paste(names(env_dat_Col)[2], "~", paste(names(env_dat_Col_ac)[-c(1,2,4,5,6:15,16,17,19,20,25,26,28:32,
                                                                         34,35:38,40,41,43:45, 46,48,50,
                                                                         55,56,57:59,61:63,64:71,73:75)], collapse="+"))  #uses same covariates as in the stepwise model
Col_defor.RAC <- glm(g_col, data=env_dat_Col_ac, family="binomial")
summary(Col_defor.RAC)
Col_defor.RAC.cov <- vcovHC(Col_defor.RAC, type = "HC3");Col_defor.RAC.rob.std.err <- sqrt(diag(Col_defor.RAC.cov)) 

#compare these side-by-side:
stargazer(defor_Col.m, Col_defor.step, Col_defor.naive, Col_defor.RAC, type = "html", no.space= TRUE, title = "Alternative Models",
          column.labels = c('Full', 'Stepwise', 'Naive', 'RAC'), 
          se=list(defor_Col.m.rob.std.err, Col_defor.step.rob.std.err, Col_defor.naive.rob.std.err, Col_defor.RAC.rob.std.err),
          digits = 3, model.numbers=FALSE, df=FALSE, out="Alt_models_Colombia_v5.html")


######
# EVAL
######

# generate a test data set of structure: ID, observation, prediction (where ID is the row number of the dataset) 
test.defor <- data.frame(cbind(as.numeric(rownames(env_dat_Col)), env_dat_Col$Deforestation, defor_Col.m$fitted.values, Col_defor.step$fitted.values, Col_defor.naive$fitted.values, Col_defor.RAC$fitted.values))
names(test.defor)<-c("ID","Observed","Predicted.Full", "Predicted.Step", "Predicted.Naive", "Predicted.RAC")
test.defor$Observed[test.defor$Observed==1] <- 0; test.defor$Observed[test.defor$Observed==2] <- 1 # not sure why, but these are pulled from factors to numbers and assigned values of 1 and 2 rather than 0 and 1
glm.cut <- as.numeric(optimal.thresholds(test.defor, opt.methods=10, req.sens = 0.90)[2]) ## Select a threshold where sensitivity is at least 90%
defor.acc <- presence.absence.accuracy(test.defor, threshold=glm.cut, st.dev=F)
defor.acc$TSS <- defor.acc$sensitivity+defor.acc$specificity-1  # add true skills statistic to our accuracy metrics
defor.acc
auc.roc.plot(test.defor, main = "Colombia Model ROC Plots", color = c('Red', 'Blue', 'Green', 'Yellow')) # plot ROC curve

###########
# RANDOM CV 
###########

#create training and testing sets:
Train <- createDataPartition(env_dat_Col$Deforestation, p=0.7, list=FALSE)
training <- env_dat_Col[ Train, ]; testing <- env_dat_Col[ -Train, ]

# create confusion matrix, estimate Auc and plot Roc
ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, savePredictions = TRUE, classProbs = TRUE, summaryFunction = twoClassSummary)
levels(training$Deforestation) <- c("zero", "one"); levels(testing$Deforestation) <- c("zero", "one")

mod_fit <- train(Deforestation ~Slope+Main_Cities_Dist
                 + Rivers_Dist + Rd_Prim_Sec_Dist +Rd_Ter_Loc_Dist + Population_Dist + Tr_Oil_Palm_Dist
                 + Ld_Woodland_Mosaic_Dist+Ld_Savanna_Dist+Ld_Dryland_Cropland_Pasture_Dist
                 + Rd_Ter_Loc_Km + Operated_Contracts_SqKm + PAs_SqKm + Tr_Banana_SqKm
                 + Tr_Gr_Palms_SqKm + Ld_Irrigated_Cropland_and_Pasture_SqKm
                 + Ld_Savanna_SqKm
                 , data=training, method="glm", family="binomial", trControl = ctrl, metric="ROC")
summary(mod_fit) 
mean(mod_fit$resample$ROC); sd(mod_fit$resample$ROC) # mean and st dev of AUC across all folds
mean(mod_fit$resample$Sens); sd(mod_fit$resample$Sens) # mean and st dev of sensitivity across all folds
mean(mod_fit$resample$Spec); sd(mod_fit$resample$Spec) # mean and st dev of specificity across all folds

pred <- predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$Deforestation)


############
#spatial CV 
############

#  Generate spatial clusters to be used in spatial K-fold Cross-Validation (Spatial K-fold CV)
mdist <- dist(env_dat_Col[c("Lattitude", "Longitude")])
hc <- hclust(mdist, method='complete')
d=150*1000  # the maximum distance between pixels within clusters (150 m * 1000 m = 150 km)                   
env_dat_Col$Clust_150km = cutree(hc, k=10, h=d) 
plot(env_dat_Col$Lattitude, env_dat_Col$Longitude, col=as.factor(env_dat_Col$Clust_150km),  xlab="Longitude", ylab="Latitude")
#legend("topright", legend=paste('group', 1:10), pch=16, col=(1:10), cex=0.8) 
### Added Line of Code; please check if it is ok


#SPATIAL K-fold CV
h <- paste(names(env_dat_Col)[2], "~", paste(names(env_dat_Col)[-c(1,2,4,5,6:15,16,17,19,20,25,26,28:32,
                                                                   34,35:38,40,41,43:45, 46,48,50,
                                                                   55,56,57:59,61:63,64:71,73:75)], collapse="+")) #uses same covariates as in the stepwise model

out.glm <- data.frame()
for (j in 1:10)
{
  id = which(env_dat_Col$Clust_150km == j)
  training_data= env_dat_Col[-id,]
  validation_data = env_dat_Col[id,]  
  
  model <- glm(h, data=training_data, family="binomial")
  
  # prediction
  pred_obj = predict(model, newdata = validation_data, type="response")
  PRED.TEST = data.frame(id_cell = validation_data$ID,
                         deforestation = validation_data$Deforestation,
                         pred = pred_obj)
  PRED.TEST$deforestation<- as.numeric(PRED.TEST$deforestation)
  PRED.TEST$deforestation[PRED.TEST$deforestation==1] <- 0; PRED.TEST$deforestation[PRED.TEST$deforestation==2] <- 1
  #performance
  # thresh <- (optimal.thresholds(PRED.TEST, opt.methods="MaxSens+Spec")[2])
  acc <- presence.absence.accuracy(PRED.TEST, threshold = 0.31, find.auc = TRUE, st.dev=TRUE)
  acc$TSS <- acc$sensitivity+acc$specificity-1
  acc$model <- j
  out.glm <- rbind(out.glm, acc)
  
  # store results
  if(j == 1) { PRED.TEST_tot = PRED.TEST }
  if(j > 1)  { PRED.TEST_tot = rbind(PRED.TEST_tot, PRED.TEST) }
  
} # end of spatial K-fold CV loop
out.glm
#PRED.TEST_tot  






