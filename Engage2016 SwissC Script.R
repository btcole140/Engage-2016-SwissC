setwd("/Users/karensamis/Google Drive/NSERC Engage/Methods and Data/SwissChard")

EngSC <- read.csv("Engage2016_SwissChard.csv")

str(EngSC)

EngSC$Rep <- as.factor(EngSC$Rep)
EngSC$Trtmt <- as.factor(EngSC$Trtmt)
EngSC$Tray <- as.factor(EngSC$Tray)
EngSC$Row <- as.factor(EngSC$Row)
EngSC$Col <- as.factor(EngSC$Col)
EngSC$Survival <- as.factor(EngSC$Survival)
EngSC$HvstDate <- as.Date(EngSC$HvstDate, "%d-%b-%y")

#Subset data
EngSC1 <- subset(EngSC, Rep == "1")
EngSC2 <- subset(EngSC, Rep == "2")
EngSC3 <- subset(EngSC, Rep == "3")

#incase you need to remove NAs... >EngSC$WetWt[!is.na(EngSC$WetWt)]

#**************************Rep 1 only************************************
##histograms and transformations
#Lvs
par(mfrow = c(1,1))
hist(EngSC1$Lvs) #has a good overall shape, but gaps on either side of tallest column **
EngSC1$LogLvs <- log10(EngSC1$Lvs)
hist(EngSC1$LogLvs) #not much better than raw
EngSC1$sqrtLvs <- sqrt(EngSC1$Lvs+0.5)
hist(EngSC1$sqrtLvs) #more centred than raw, but does not deal with gaps
EngSC1$rankLvs <- rank(EngSC1$Lvs)

#WetWt
par(mfrow = c(1,1))
hist(EngSC1$WetWt) #great shape **
EngSC1$LogWetWt <- log10(EngSC1$WetWt+1)
hist(EngSC1$LogWetWt) #not better than raw
EngSC1$sqrtWetWt <- sqrt(EngSC1$WetWt+0.5)
hist(EngSC1$sqrtWetWt) #not better than raw
EngSC1$rankWetWt <- rank(EngSC1$WetWt)

#DryWt
par(mfrow = c(1,1))
hist(EngSC1$DryWt) #great shape - not many columns
EngSC1$LogDryWt <- log10(EngSC1$DryWt+1)
hist(EngSC1$LogDryWt) #not better than raw
EngSC1$sqrtDryWt <- sqrt(EngSC1$DryWt+0.5)
hist(EngSC1$sqrtDryWt) #better than raw because has more columns **
EngSC1$rankDryWt <- rank(EngSC1$DryWt)

#WaterWt
par(mfrow = c(1,1))
hist(EngSC1$WaterWt) #great shape **
EngSC1$LogWaterWt <- log10(EngSC1$WaterWt+1)
hist(EngSC1$LogWaterWt) #not better than raw
EngSC1$sqrtWaterWt <- sqrt(EngSC1$WaterWt+0.5)
hist(EngSC1$sqrtWaterWt) #not better than raw
EngSC1$rankWaterWt <- rank(EngSC1$WaterWt)


##AOV
#Lvs
aovSC1Lvs  <- aov(Lvs~Trtmt, data=EngSC1)
summary(aovSC1Lvs)
plot(aovSC1Lvs) # plots look ok
  # no sig diff in number of leaves between treatments (p=0.324 f=1.162)
boxplot(Lvs~Trtmt, data=EngSC1)
TukeyHSD(aovSC1Lvs) 
#             diff         lwr       upr     p adj
#0.5-0 -0.10000000 -0.39014148 0.1901415 0.8108086
#2-0    0.10833333 -0.18180814 0.3984748 0.7707473
#5-0    0.02500000 -0.26514148 0.3151415 0.9961290
#2-0.5  0.20833333 -0.08180814 0.4984748 0.2508029
#5-0.5  0.12500000 -0.16514148 0.4151415 0.6831048
#5-2   -0.08333333 -0.37347481 0.2068081 0.8806811

#WetWt
aovSC1WetWt  <- aov(WetWt~Trtmt, data=EngSC1)
summary(aovSC1WetWt)
plot(aovSC1WetWt) # plots look ok
  # sig diff in wet weights between treatments (p=<0.0001 f=17.8)
boxplot(WetWt~Trtmt, data=EngSC1)
TukeyHSD(aovSC1WetWt) 
#            diff        lwr         upr     p adj
#0.5-0 -0.6700833 -1.2226583 -0.11750838 0.0101113 *
#2-0    0.7480000  0.1954250  1.30057496 0.0029590 **
#5-0   -0.5154167 -1.0679916  0.03715829 0.0775328
#2-0.5  1.4180833  0.8655084  1.97065829 0.0000000 ***
#5-0.5  0.1546667 -0.3979083  0.70724162 0.8884633
#5-2   -1.2634167 -1.8159916 -0.71084171 0.0000000 ***

#DryWt (sqrt)
aovSC1DryWt  <- aov(sqrtDryWt~Trtmt, data=EngSC1)
summary(aovSC1DryWt)
plot(aovSC1DryWt) # plots look ok
  # sig diff in wet weights between treatments (p=<0.0001 f=8.376)
boxplot(sqrtDryWt~Trtmt, data=EngSC1)
TukeyHSD(aovSC1DryWt) 
#              diff         lwr           upr     p adj
#0.5-0 -0.024358358 -0.05427886  0.0055621423 0.1549300
#2-0    0.022370811 -0.00754969  0.0522913108 0.2178329
#5-0   -0.028942781 -0.05886328  0.0009777197 0.0621350
#2-0.5  0.046729169  0.01680867  0.0766496688 0.0003824 ***
#5-0.5 -0.004584423 -0.03450492  0.0253360776 0.9790858
#5-2   -0.051313591 -0.08123409 -0.0213930909 0.0000716 ***

#WaterWt
aovSC1WaterWt  <- aov(WaterWt~Trtmt, data=EngSC1)
summary(aovSC1WaterWt)
plot(aovSC1WaterWt) # plots look ok
  # sig diff in wet weights between treatments (p=<0.0001 f=17.86)
boxplot(WaterWt~Trtmt, data=EngSC1)
TukeyHSD(aovSC1WaterWt) 
#            diff        lwr         upr     p adj
#0.5-0 -0.6074417 -1.1089405 -0.10594280 0.0102289 **
#2-0    0.6928917  0.1913928  1.19439054 0.0022833 **
#5-0   -0.4436083 -0.9451072  0.05789054 0.1040223
#2-0.5  1.3003333  0.7988345  1.80183220 0.0000000 ***
#5-0.5  0.1638333 -0.3376655  0.66533220 0.8342487
#5-2   -1.1365000 -1.6379989 -0.63500113 0.0000001 ***
