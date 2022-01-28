library(MASS)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)
library(reshape2)
library(plyr)
library(factoextra)
library(randomForest)
library(permimp)

raw.data<-read.csv(url("https://github.com/Dansieg91/AMSA-KUL-DCT/raw/main/amsa_data_dct_p.csv"))
raw.data$TREAT<-as.factor(raw.data$TREAT)
str(raw.data)


var_list<-c('INSTN',"TREAT","x40","x13","x3","x25","x14","x17","x31","Pwev","x7",
            "x5","x22","x16","x15","x47","x48","x51","x52","x43","x41","x53",
            "x45","x46","x63","x64","x60","x26","x19","x30")

df <- raw.data[,names(raw.data) %in% var_list]
                                         
# creating average by genotype within each treatment
df_ni <- df[df$TREAT == "NI", ]
df_rec <- df[df$TREAT == "REC", ]
df_td <- df[df$TREAT == "TD", ]
dim(df_ni)
dim(df_rec)
dim(df_td)

df_ni <- aggregate(df_ni[, 3:27], list(df_ni$INSTN), FUN=mean)
df_ni$treat<-"NI"
df_rec <- aggregate(df_rec[, 3:27], list(df_rec$INSTN), FUN=mean)
df_rec$treat<-"REC"
df_td <- aggregate(df_td[, 3:27], list(df_td$INSTN), FUN=mean)
df_td$treat<-"TD"



# Merging all averages by treatment into a single df
df_avg<-rbind(df_ni,df_rec,df_td)
names(df_avg)[names(df_avg) == 'Group.1'] <- 'INSTN'

# Checking missing values in all variables
apply(is.na(df_avg),2,sum)


#fixing the levels
levels(as.factor(df_avg$Pwev))
df_avg$Pwev<-ifelse(df_avg$Pwev >= 4.4 & df_avg$Pwev <=5.8,5,df_avg$Pwev)
df_avg$Pwev<-ifelse(df_avg$Pwev >= 6 & df_avg$Pwev <=7.8,7,df_avg$Pwev)
df_avg$Pwev<-ifelse(df_avg$Pwev >= 8 & df_avg$Pwev <=8.8,9,df_avg$Pwev)
levels(as.factor(df_avg$Pwev))


############################################################
###########     Fitting LDA - first approach    ############
############################################################
df_lda0 <-na.omit(df_avg)
df_lda1 <- as.data.frame(scale(df_lda0[,2:26],center = T, scale = T))
df_lda1<-as.data.frame(cbind(df_lda0[,5],df_lda1))
df_lda1 <- df_lda1[,!names(df_lda1) %in% c('Pwev',"x26","x19","x30")]
names(df_lda1)[names(df_lda1) == 'df_lda0[, 5]'] <- 'Pwev'
str(df_lda1)


set.seed(1)
intrain.1 <- sample(nrow(df_lda1), round(0.50*nrow(df_lda1)))
train.1 <- df_lda1[intrain.1, ]
test.1 <- df_lda1[-intrain.1, ]


#### train ####
(lda1.tr <- lda(as.factor(Pwev)~ ., data=df_lda1, CV = F))
lda1.ptr <- predict(lda1.tr, newdata =train.1 )
sum(diag(prop.table(table(train.1$Pwev, lda1.ptr$class))))
table(train.1$Pwev, lda1.ptr$class)
#plot(lda1.tr,col=as.numeric(train.1$Pwev))

#------- test ---------#
lda1.pte <- predict(lda1.tr, newdata=test.1)
# Results table
tr.lda1<- table(test.1$Pwev,lda1.pte$class)
tr.lda1<-as.data.frame.matrix(tr.lda1)
tr.lda1[c("x1","Orig")]<-count(as.numeric(test.1$Pwev))
tr.lda1[c("x2","Pred")]<-count(lda1.pte$class)
tr.lda1<-tr.lda1[,-c(6,8)]
tr.lda1[c("h_1","h_3","h_5","h_7","h_9")]<-(round((as.data.frame.matrix(prop.table(table(test.1$Pwev, lda1.pte$class)))),4))
tr.lda1$Hit_rate<-c("","","","",round((sum(diag(prop.table(table(test.1$Pwev, lda1.pte$class))))),4))
tr.lda1
#plot(lda1.tr,col=as.numeric(test.1$Pwev))


#--------- LOCV ---------#
lda1.cv <- lda(as.factor(Pwev)~ ., data=df_lda1, CV = T)
sum(diag(prop.table(table(df_lda1$Pwev, lda1.cv$class))))
table(df_lda1$Pwev, lda1.cv$class)

# Results table
cv.lda1<- table(df_lda1$Pwev, lda1.cv$class)
cv.lda1<-as.data.frame.matrix(cv.lda1)
cv.lda1[c("x1","Orig")]<-count(as.numeric(df_lda1$Pwev))
cv.lda1[c("x2","Pred")]<-count(lda1.cv$class)
cv.lda1<-cv.lda1[,-c(6,8)]
cv.lda1[c("h_1","h_3","h_5","h_7","h_9")]<-(round((as.data.frame.matrix(prop.table(table(df_lda1$Pwev, lda1.cv$class)))),4))
cv.lda1$Hit_rate<-c("","","","",round((sum(diag(prop.table(table(df_lda1$Pwev, lda1.cv$class))))),4))
cv.lda1

str(df_lda1)
### Making plots
# Color item
#library(RColorBrewer)
myColors <- colors()[c(258, 496, 87, 144, 75)]
df_lda1$Pwev <- factor(test.1$Pwev, levels=c("No Wilting", "Slight", "Medium", "Severe", "Very Severe"))
(names(myColors) <- levels(as.factor(test.1$Pwev)))
colScale <- scale_colour_manual(name = "Plant Wilting",values = myColors)

# Train PLOT
prop.lda =  lda1.tr$svd^2/sum(lda1.tr$svd^2)
plda1 <- predict(object = lda1.tr, newdata = test.1)
dataset1 = data.frame(P.Wilt = as.factor(test.1[,"Pwev"]), lda = plda1$x)
p1 <- ggplot(dataset1) + geom_point(aes(lda.LD1, lda.LD2, colour = P.Wilt), size = 2) + xlim(-6, 6) + ylim(-6, 6) +
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) +
  ggtitle(' ') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
grid.arrange(p1 + colScale)





#################################################
########### Fitting LDA Second approach #########
#################################################

df_lda0 <-na.omit(df_avg)
df_lda2 <- as.data.frame(scale(df_lda0[,2:26],center = T, scale = T))
df_lda2<-as.data.frame(cbind(df_lda0[,27],df_lda2))
df_lda2 <- df_lda2[,!names(df_lda2) %in% c('treat','Pwev',"x26","x19","x30")]
names(df_lda2)[names(df_lda2) == 'df_lda0[, 27]'] <- 'treat'

set.seed(1)
intrain.2 <- sample(nrow(df_lda2), round(0.50*nrow(df_lda2)))
train.2 <- df_lda2[intrain.2, ]
test.2 <- df_lda2[-intrain.2, ]
str(test.2)

#### train ####
(lda2.tr <- lda(as.factor(treat)~ ., data=df_lda2, CV = F))
lda2.ptr <- predict(lda2.tr, newdata =train.2 )
sum(diag(prop.table(table(train.2$treat, lda2.ptr$class))))
table(train.2$treat, lda2.ptr$class)
#plot(lda2.tr,col=as.numeric(train.2$treat))

#------- test ---------#
lda2.pte <- predict(lda2.tr, newdata=test.2)
# Results table
tr.lda2<- table(test.2$treat,lda2.pte$class)
tr.lda2<-as.data.frame.matrix(tr.lda2)
tr.lda2[c("x1","Orig")]<-count(test.2$treat)
tr.lda2[c("x2","Pred")]<-count(lda2.pte$class)
tr.lda2<-tr.lda2[,-c(4,6)]
tr.lda2[c("h_NI","h_REC","h_TD")]<-(round((as.data.frame.matrix(prop.table(table(test.2$treat, lda2.pte$class)))),4))
tr.lda2$Hit_rate<-c("","",round((sum(diag(prop.table(table(test.2$treat, lda2.pte$class))))),4))
tr.lda2

plot(lda2.tr,col=as.numeric(test.2$treat))


#--------- LOCV ---------#
lda2.cv <- lda(as.factor(treat)~ ., data=df_lda2, CV = T)
sum(diag(prop.table(table(df_lda2$treat, lda2.cv$class))))
table(df_lda2$treat, lda2.cv$class)
str(df_lda2)


# Results table
cv.lda2<- table(df_lda2$treat, lda2.cv$class)
cv.lda2<-as.data.frame.matrix(cv.lda2)
cv.lda2[c("x1","Orig")]<-count(test.2$treat)
cv.lda2[c("x2","Pred")]<-count(lda2.cv$class)
cv.lda2<-cv.lda2[,-c(4,6)]
cv.lda2[c("h_NI","h_REC","h_TD")]<-(round((as.data.frame.matrix(prop.table(table(df_lda2$treat, lda2.cv$class)))),4))
cv.lda2$Hit_rate<-c("","",round((sum(diag(prop.table(table(df_lda2$treat, lda2.cv$class))))),4))
cv.lda2


### Making plots
# Color item
#library(RColorBrewer)
myColors <- colors()[c(258, 144, 75)]
df_lda2$treat <- factor(test.2$treat, levels=c("Control", "Medium", "Severe"))
(names(myColors) <- levels(as.factor(test.2$treat)))
colScale <- scale_colour_manual(name = "treat",values = myColors)

# Train PLOT
prop.lda =  lda2.tr$svd^2/sum(lda2.tr$svd^2)
plda2 <- predict(object = lda2.tr, newdata = test.2)
dataset2 = data.frame(treat = as.factor(test.2[,"treat"]), lda = plda2$x)
p2 <- ggplot(dataset2) + geom_point(aes(lda.LD1, lda.LD2, colour = treat), size = 2) + xlim(-6, 6) + ylim(-6, 6) +
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) +
  ggtitle(' ') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
grid.arrange(p2 + colScale)



##### Assesing importance ####

# with LDA
lda2.tr

# With Random forest (RF)
r_forest <- randomForest(ordered(treat) ~ ., data = train.2)
pred3 <- predict(r_forest, newdata = test.2)
(ct <- table(test.2$treat, pred3))
sum(diag(prop.table(ct)))
gini <- data.frame(r_forest$importance)
gini$trait <- rownames(gini)
gini[sort.int(gini$MeanDecreaseGini, decreasing = T, index.return = T)$ix, ]

# optional RF with permutation (package mentions results are similar as in regular RF varimp)
perm<-permimp(r_forest,conditional = T)



###################################
#####    LDA third approach    ####
###################################

df_lda0 <-na.omit(df_avg)
df_lda3 <- as.data.frame(scale(df_lda0[,2:26],center = T, scale = T))
df_lda3<-as.data.frame(cbind(df_lda0[,27],df_lda3))
df_lda3 <- df_lda3[,!names(df_lda3) %in% c('treat','PwevDAP')]
names(df_lda3)[names(df_lda3) == 'df_lda0[, 27]'] <- 'treat'

set.seed(1)
intrain.3 <- sample(nrow(df_lda3), round(0.50*nrow(df_lda3)))
train.3 <- df_lda3[intrain.3, ]
test.3 <- df_lda3[-intrain.3, ]
str(test.3)

#### train ####
(lda3.tr <- lda(as.factor(treat)~ x26 + x43 + x19 + x30, data=df_lda3, CV = F))
lda3.ptr <- predict(lda3.tr, newdata =train.3 )
sum(diag(prop.table(table(train.3$treat, lda3.ptr$class))))
table(train.3$treat, lda3.ptr$class)
plot(lda3.tr,col=as.numeric(train.3$treat))


#------- test ---------#
lda3.pte <- predict(lda3.tr, newdata=test.3)
# Results table
tr.lda3<- table(test.3$treat,lda3.pte$class)
tr.lda3<-as.data.frame.matrix(tr.lda3)
tr.lda3[c("x1","Orig")]<-count(test.3$treat)
tr.lda3[c("x2","Pred")]<-count(lda3.pte$class)
tr.lda3<-tr.lda3[,-c(4,6)]
tr.lda3[c("h_NI","h_REC","h_TD")]<-(round((as.data.frame.matrix(prop.table(table(test.3$treat, lda3.pte$class)))),4))
tr.lda3$Hit_rate<-c("","",round((sum(diag(prop.table(table(test.3$treat, lda3.pte$class))))),4))
tr.lda3

#plot(lda3.tr,col=as.numeric(test.3$treat))

#--------- LOCV ---------#
lda3.cv <- lda(as.factor(treat)~ x26 + x43 + x19 + x30, data=df_lda3, CV = T)
sum(diag(prop.table(table(df_lda3$treat, lda3.cv$class))))
table(df_lda3$treat, lda3.cv$class)
str(df_lda3)


# Results table
cv.lda3<- table(df_lda3$treat, lda3.cv$class)
cv.lda3<-as.data.frame.matrix(cv.lda3)
cv.lda3[c("x1","Orig")]<-count(test.3$treat)
cv.lda3[c("x2","Pred")]<-count(lda3.cv$class)
cv.lda3<-cv.lda3[,-c(4,6)]
cv.lda3[c("h_NI","h_REC","h_TD")]<-(round((as.data.frame.matrix(prop.table(table(df_lda3$treat, lda3.cv$class)))),4))
cv.lda3$Hit_rate<-c("","",round((sum(diag(prop.table(table(df_lda3$treat, lda3.cv$class))))),4))
cv.lda3


### Making plots
# Color item
#library(RColorBrewer)
myColors <- colors()[c(258, 144, 75)]
df_lda3$treat <- factor(test.3$treat, levels=c("Control", "Medium", "Severe"))
(names(myColors) <- levels(as.factor(test.3$treat)))
colScale <- scale_colour_manual(name = "Treat",values = myColors)

# Train PLOT
prop.lda =  lda3.tr$svd^2/sum(lda3.tr$svd^2)
plda3 <- predict(object = lda3.tr, newdata = test.3)
dataset3 = data.frame(treat = as.factor(test.3[,"treat"]), lda = plda3$x)
p3 <- ggplot(dataset3) + geom_point(aes(lda.LD1, lda.LD2, colour = treat), size = 2) + xlim(-6, 6) + ylim(-6, 6) +
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) +
  ggtitle(' ') +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
grid.arrange(p3 + colScale)


















































