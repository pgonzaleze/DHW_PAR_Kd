#############################################################################
#############################################################################
################  DETERMINE CLUSTERS USING K-MEDOIDS   ######################
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
library(tidyverse)
library(cluster)
library(factoextra)
#library(devtools)
#install_github("vqv/ggbiplot")

CCB = read.csv("DCW_PAR_Kd.csv") # import dataset
CCB$site_year <- paste(CCB$site_name,CCB$year) #combine collumns site and year and create a new collumn
CCB$bleaching_status <- ifelse(CCB$bleaching== 0, "No", "Yes") 
CCB[,14] <- as.factor(CCB[,14])  # the "response" should be binary and as factor 
CCB <-subset(CCB, DCW != 0)   ## clean al DCW = 0 since there are no chances to have bleaching with no thermal stress

df <- CCB
df <- df %>%
  select(DCW, dDLWln) 
df <- na.omit(df)

# Determine best number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# K-means cluster
k2 <- kmeans(df, centers = 2, nstart = 25)
#str(k2)
fviz_cluster(k2, data = df)

# PAM 
pam <- pam((scale(df)),stand = TRUE, 2) 
(pam_plot <- fviz_cluster(pam, ellipse.type = 'convex', geom = 'point') +
    scale_color_manual(values=c("gray1", "black")) +
    scale_shape_manual(values=c(1,19)) +
    scale_fill_manual(values=c("blue", "white")) +
    theme_classic())
pam_plot + theme(axis.text.x = element_text(size = 18),
                 axis.text.y = element_text(size = 18))

#### PCA ####
CCB.pca <- prcomp(CCB[,c(5,9)], center = TRUE,scale. = TRUE)
summary(CCB.pca)
fviz_eig(CCB.pca)
fviz_pca_biplot(CCB.pca)
