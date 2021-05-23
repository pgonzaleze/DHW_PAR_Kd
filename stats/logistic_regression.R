#############################################################################
#############################################################################
##########  LOGISTIC REGRESSION MODEL TO PREDICT CORAL BLEACHING  ###########
#############################################################################
#############################################################################

'
Author: Pedro C Gonzalez-Espinosa
University of British Columbia
Departent of Geography
Date: Aug/08/18
'

############################################################################
########################## loading libraries ###############################
############################################################################

library(caTools) # library to split the database in test and training
library(ggplot2)  # For plotting 
library(ROCR) # library to plot cutoffs of logistic regression
library(rcompanion) # to get the pseudoR2
library(caret) # For CrossValidation using "repeatedcv" method. Other option include "cv" and "LOOCV"
library(corrplot) # For correlation plot
library(plotly) # for 3D plot 

#### loading .csv file #### 
CCB = read.csv("DCW_PAR_Kd.csv") # import dataset
CCB$site_year <- paste(CCB$site_name,CCB$year) #combine columns site and year and create a new column
CCB$bleaching_status <- ifelse(CCB$bleaching == 0, "No", "Yes") 
CCB[,22] <- as.factor(CCB[,22])  # the "response" should be binary and as factor 
CCB <-subset(CCB, DCW != 0)   ## clean all DCW = 0 since there are no chances to have
                              ## a bleaching with no thermal stress

###### dimensions and structure ####
dim(CCB)  # dimensions 
str(CCB)  # structure of the dataframe

###########################################################################
####################### split data in train and test ######################
###########################################################################

## First glimpse of logistic regression model
## Two options:
set.seed(11) # this is applied for replicatibility

#1)
row.number = sample(1:nrow(CCB), 0.8*nrow(CCB)) # split the data to fit 5-fold
train = CCB[row.number,]                        
test = CCB[-row.number,]
dim(train)
dim(test)

#2) slpliting using library caTools (alternative method)
split <- sample.split(CCB, SplitRatio = 0.8)
train = subset(CCB, split == TRUE)
test = subset(CCB, split == FALSE)

attach(train)
model1 = glm(bleaching ~ CS,
             data=train, family=binomial(link = 'logit'))
model2 = glm(bleaching ~ DCW,
             data=train, family=binomial(link = 'logit'))
model3 = glm(bleaching ~ CS + LS,
             data=train, family=binomial(link = 'logit'))
model4 = glm(bleaching ~ DCW + LS,
             data=train, family=binomial(link = 'logit'))# quit a variable is only an option if AIC decreases
model5 = glm(bleaching ~ CS + DLW,
             data=train, family=binomial(link = 'logit'))
model6 = glm(bleaching ~ DCW + DLW,
             data=train, family=binomial(link = 'logit'))
model7 = glm(bleaching ~ CS + dDLW,
             data=train, family=binomial(link = 'logit'))
model8 = glm(bleaching ~ DCW + dDLW,
             data=train, family=binomial(link = 'logit'))
model9 = glm(bleaching ~ CS + LS + Kd490a,
             data=train, family=binomial(link = 'logit'))
model10 = glm(bleaching ~ DCW + LS + Kd490a,
             data=train, family=binomial(link = 'logit'))
#### confirm and use the one with better performance (AIC and R2)####
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
summary(model6)
summary(model7)
summary(model8)
summary(model9)
summary(model10)
nagelkerke(model1) # gives the pseudoR2
nagelkerke(model2)
nagelkerke(model3)
nagelkerke(model4)
nagelkerke(model5)
nagelkerke(model6)
nagelkerke(model7)
nagelkerke(model8) # *
nagelkerke(model9)
nagelkerke(model10)
#### Defining the threshold ####
res <- predict(model8, train, type = "response")
ROCRPred <- prediction(res, train$bleaching)
ROCRPerf <- performance(ROCRPred, "tpr", "fpr") #tpr=true positive and fpr=false positive 
plot(ROCRPerf, colorize=TRUE, print.cutoffs.at=seq(0.1, by=0.1))

#### Predict for training data and find training accuracy ####
pred.prob = predict(model2, type="response")
pred.prob = ifelse(pred.prob > 0.5, 1, 0)
table(pred.prob, train$bleaching) ## Iterate for each model

#### Predict for test Data and find the test accuracy ####
attach(test)
pred.prob = predict(model2, newdata= test, type="response")
pred.prob = ifelse(pred.prob > 0.5, 1, 0)
table(pred.prob, test$bleaching)

#### Predict for training data and find training accuracy ####
y_pred <- factor(pred.prob, levels=c(0, 1))
y_act <- bleaching
mean(y_pred == y_act)

####verify which cases were correctly classified ####
attach(test)
pred.probM1 = predict(model1, newdata= test, type="response")
pred.probM1 = ifelse(pred.probM1 > 0.5, 1, 0)
pred.probM1  ###do the same for all models iterating "M2..M3..M4"
pred.probM2 = predict(model2, newdata= test, type="response")
pred.probM2 = ifelse(pred.probM2 > 0.5, 1, 0)
pred.probM2
pred.probM6 = predict(model6, newdata= test, type="response")
pred.probM6 = ifelse(pred.probM6 > 0.5, 1, 0)
pred.probM6
pred.probM8 = predict(model8, newdata= test, type="response")
pred.probM8 = ifelse(pred.probM8 > 0.5, 1, 0)
pred.probM8

###########################################################################
##############################   DCW ONLY #################################
###########################################################################

### Predict for training data and find training accuracy
pred.prob1 = predict(model1, type="response")
pred.prob1 = ifelse(pred.prob1 > 0.4, 1, 0)
table(pred.prob1, train$bleaching)

## Predict for test Data and find the test accuracy.
attach(test)
pred.prob1 = predict(model1, newdata= test, type="response")
pred.prob1 = ifelse(pred.prob1 > 0.5, 1, 0)
table(pred.prob1, test$bleaching)

###Predict for training data and find training accuracy
yy_pred <- factor(pred.prob1, levels=c(0, 1))
yy_act <- test$bleaching
mean(yy_pred == yy_act)

##### comparing models using ANOVA 
anova(model8, model2, test="LRT")

###########################################################################
####################### Cross-validation analysis #########################
###########################################################################

# 3 options 
#### 1) K-fold Cross validation
# define training control
train_control <- trainControl(method="cv", number=5,
                              savePredictions = TRUE)
#### 2) Repeated K-fold cross validation 
# define training control
train_control <- trainControl(method="repeatedcv", number=5, repeats = 200,
                              savePredictions = TRUE)
#### 3) LOOCV 
# define training control
train_control <- trainControl(method="LOOCV", savePredictions = TRUE)

# train the model   # regardless the selected option 
model1 <- train(bleaching ~ CS, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model2 <- train(bleaching ~ DCW, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model3 <- train(bleaching ~ DCW1_5, data=CCB,
              trControl=train_control, method="glm", family="binomial")
model4 <- train(bleaching ~ DCW2, data=CCB,
              trControl=train_control, method="glm", family="binomial")
model5 <- train(bleaching ~ CSM.CS, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model6 <- train(bleaching ~ DCW.MCS0, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model7 <- train(bleaching ~ DCW + DLW, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model8 <- train(bleaching ~ DCW + dDLW, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model9 <- train(bleaching ~ DCW + dDLW2m, data=CCB,
                trControl=train_control, method="glm", family="binomial")    
model10 <- train(bleaching ~ DCW + dDLW5m, data=CCB,
                trControl=train_control, method="glm", family="binomial")   
model11 <- train(bleaching ~ DCW + dDLWsite, data=CCB,
                trControl=train_control, method="glm", family="binomial")                                        
model12 <- train(bleaching ~ DCW + LS + Kd490a , data=CCB,
                trControl=train_control, method="glm", family="binomial")
model13 <- train(bleaching ~ CountDays , data=CCB,
                trControl=train_control, method="glm", family="binomial")
model14 <- train(bleaching ~ consecutive_days , data=CCB,
                 trControl=train_control, method="glm", family="binomial")
model15 <- train(bleaching ~ CountDays + CPAR, data=CCB,
                 trControl=train_control, method="glm", family="binomial")
model16 <- train(bleaching ~ CountDays + CPARprevm, data=CCB,
                 trControl=train_control, method="glm", family="binomial")
model17 <- train(bleaching ~ LS, data=CCB,
                trControl=train_control, method="glm", family="binomial")
model18 <- train(bleaching ~ C.PAR , data=CCB,
                trControl=train_control, method="glm", family="binomial")

##### print cross-validation scores #####
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)
summary(model6)
summary(model7)
summary(model8)
summary(model9)
summary(model10)
summary(model11)
summary(model12)
summary(model13)
summary(model14)
summary(model15)
summary(model16)

##### gives the accuracy and kappa values #####
model1$results
model2$results
model3$results 
model4$results
model5$results
model6$results
model7$results
model8$results
model9$results
model10$results
model11$results
model12$results
model13$results
model14$results
model15$results
model16$results

######### confusion matrix of cross-validation #######
confusionMatrix(model1)
confusionMatrix(model2)
confusionMatrix(model3)
confusionMatrix(model4)
confusionMatrix(model5)
confusionMatrix(model6)
confusionMatrix(model7)
confusionMatrix(model8)
confusionMatrix(model9)
confusionMatrix(model10)
confusionMatrix(model11)
confusionMatrix(model12)
confusionMatrix(model13)
confusionMatrix(model14)
confusionMatrix(model15)
confusionMatrix(model16)

#### getting the probabilities of each "event" of the test dataset ####
# Model 6 is used as example
pred.probM6 = predict(model5, type="prob")
pred.probM6 = ifelse(pred.probM6 > 0.5, 1, 0)
pred.probM6

### Create new columns to be used in the 3D plot (this should be done separately from the previous step)
pred.probM2 = predict(model2, type="prob")
pred.probM2$no_bleaching <- pred.probM2[,1]
pred.probM2$bleaching <- pred.probM2[,2]
pred.probM2

pred.probM3 = predict(model3, type="prob")
pred.probM3$no_bleaching <- pred.probM3[,1]
pred.probM3$bleaching <- pred.probM3[,2]
pred.probM3

pred.probM4 = predict(model4, type="prob")
pred.probM4$no_bleaching <- pred.probM4[,1]
pred.probM4$bleaching <- pred.probM4[,2]
pred.probM4

pred.probM5 = predict(model5, type="prob")
pred.probM5$no_bleaching <- pred.probM5[,1]
pred.probM5$bleaching <- pred.probM5[,2]
pred.probM5

############################################################################
############################### ploting LOGIT ##############################
############################################################################

#### ####
CCBplot = read.csv("DCW_PAR_Kd.csv") # import dataset
CCBplot[,15] <- as.numeric(CCBplot[,15])
CCBplot$site_year <- paste(CCBplot$site_name, CCBplot$year)
CCBplot <-subset(CCBplot, DCW != 0)
CCBplot <- subset(CCBplot, !site_name  %in% c('ESC','GAV', 'CAN',
                                              'CAR', 'LNV')) 
str(CCBplot)              

#### ploting usin plot() ####
DCW <- plot(CCB$DCW,pred.probM2$bleaching)
g=glm(bleaching~DCW,family=binomial,CCBplot) 
curve(predict(g,data.frame(DCW=x),type="resp"),add=TRUE)

### plotting using ggplot2
(DCW <- ggplot(CCBplot, aes(x=DCW, y=bleaching, colour=dDLWln, size=dDLWln)) + 
    geom_point(position = position_jitter(w=0, h=0)) + 
    scale_color_gradient(low="red4", high="yellow2") +
    theme_classic() + 
    stat_smooth(method="glm", method.args = list(family="binomial"),se=FALSE,
                fullrange=T, colour = "gray25"))+
  scale_x_continuous(breaks=seq(-24,0,4))


########################################################################
############### Scatterplot bleaching points filled ####################
########################################################################

# Step 1 create a "temporary data frame" for plotting
probM6 = predict(model6, type="prob")
probM6 = ifelse(probM6 > 0.5, 1, 0)
probM6 <- as.data.frame(probM6)
probM6$no_bleaching <- probM6[,1]
probM6$bleaching <- probM6[,2]
probM6
CCB$probM6 <- paste(probM6$bleaching)
CCB$bleaching_probM6 <- ifelse(CCB$probM6== 0, "No", "Yes") 

#### Plotting the probability of bleaching and its associated errors
# with label of "bleaching" on top and horizontal

spDCW_dDLWln <- ggplot(CCB, aes(DCW, dDLWln, label=site_year)) +
  geom_point(aes(shape=bleaching_status, colour=bleaching_status), size=4) +
  labs(shape="Bleaching", colour="Bleaching") +
  scale_shape_manual(values=c(1,19)) +
  scale_color_manual(values=c("gray1", "midnightblue")) +
  geom_text(data=subset(CCB, bleaching_status == "Yes"), size= 5, 
            nudge_x = 0.06, nudge_y = 0.06, check_overlap = F) +
  scale_y_continuous(name="(ln)dDLW", position="right", breaks = c(2,3,4), 
                     limits = c(1.5, 4.6)) +
  scale_x_continuous(breaks = c(-4,-8,-12,-16,-20), limits = c(-24,0)) +
  theme_classic() +
  theme(legend.position = "top", legend.box = "horizontal", 
        legend.text=element_text(size=14), 
        legend.title = element_text(size = 14))+
  theme(axis.text.x = element_text(color="gray15", size=12),
        axis.text.y = element_text(color="gray15", size=12))
spDCW_dDLWln 

spDCW_dDLW <- ggplot(CCB, aes(DCW, dDLW, label=site_year)) +
  geom_point(aes(shape=bleaching_status, colour=bleaching_status), size=4) +
  labs(shape="Bleaching", colour="Bleaching") +
  scale_shape_manual(values=c(1,19)) +
  scale_color_manual(values=c("gray1", "midnightblue")) +
  geom_text(data=subset(CCB, bleaching_status == "Yes"), size= 5, 
            nudge_x = 0.06, nudge_y = 0.06, check_overlap = F) +
  scale_y_continuous(name="dDLW", position="right", breaks = c(10,20,30,40,50), 
                     limits = c(1, 50)) +
  scale_x_continuous(breaks = c(-4,-8,-12,-16,-20), limits = c(-24,0)) +
  theme_classic() +
  theme(legend.position = "top", legend.box = "horizontal", 
        legend.text=element_text(size=14), 
        legend.title = element_text(size = 14))+
  theme(axis.text.x = element_text(color="gray15", size=12),
        axis.text.y = element_text(color="gray15", size=12))
spDCW_dDLW

#######################################################################
################## 2D & 3D plot (plotly package) ######################
#######################################################################

##Graphing your 2D & 3d scatterplot using plotly's scatter3d type:
# To have plain background first set 

# example using pred.probM2
xaxis <- list(title= "DCW", showgrid = FALSE, 
              autotick = FALSE, ticks = "inside",  tick0 = 0,  dtick = 2,
              ticklen = 5,  tickwidth = 2,  tickcolor = toRGB("black")) 
yaxis <- list(title= "Probability of bleaching", showgrid = FALSE,
              autotick = FALSE, ticks = "outside",  tick0 = 0,  dtick = 0.25,
              ticklen = 5,  tickwidth = 2,  tickcolor = toRGB("black"))

(M2_2D <- plot_ly(data=CCB, x=CCB$DCW, y=pred.probM2$bleaching,
                 color=CCB$bleaching_status, colors = "Greys",
                 type = "scatter", mode = "markers",
                 marker = list(size = 8,
                               line = list(color = 'black',
                                           width = 1))) %>%
  layout(title = " ",
         xaxis = xaxis,
         yaxis = yaxis))

# example using pred.probM4
xaxis <- list(title= "(ln)dDLW", showgrid = TRUE)
yaxis <- list(title= "DCW", showgrid = TRUE)
zaxis <- list(title= "Pblx", showgrid = TRUE)
(M4_3D <- plot_ly(x=CCB$DLWln, y=CCB$DCW, z=pred.probM4$bleaching,
                 type="scatter3d", mode="markers", color=CCB$DLWln) %>%
  layout(
    title = "Cold coral bleaching DCW + (ln)DLW",
    scene = list(
      xaxis = xaxis, #list(title = "(ln)DLW"),
      yaxis = yaxis, #list(title = "DCW"),
      zaxis = zaxis)#list(title = "Pblx")
    ))

M2_2D
M4_3D
