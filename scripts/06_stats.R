#################################
#Load data and libraries
#################################
source("scripts/04_combine-MIC-flow.R")
library(lmerTest)
#################################
#Fitness in FLC1
#################################
# The majority of evolved replicates grew to a higher optical density (OD) in 1ug/mL than the ancestral strains, when measured at either 24 h or 72 h (Figure 1, paired t-test results in Table S1).
comp1 <- list()
for(i in 1:20){
  sub0 <- subset(all0.1, line==i)
  sub10 <- subset(all10.1, line==i)
  test <- t.test(sub0$all0.1, sub10$all10.1)
  comp1[[i]] <- c(test$statistic, test$parameter, test$p.value)
}
FLC1.ttest <- do.call("rbind", comp1)
FLC1.ttest <- data.frame(line =1:20, FLC1.ttest)
names(FLC1.ttest)[4] <- "p.value"
FLC1.ttest$sig <- ifelse(FLC1.ttest$p.value < 0.05, "*", "")
FLC1.ttest$t <- round(FLC1.ttest$t, 2)
FLC1.ttest$df <- round(FLC1.ttest$df, 2)
FLC1.ttest$p.value <- round(FLC1.ttest$p.value, 4)
write.csv(FLC1.ttest, "manuscript/tables/TablsS1a-FLC1ttest_24h.csv", row.names=FALSE)

comp1.72 <- list()
for(i in 1:20){
  sub0 <- subset(all0.1.72, line==i)
  sub10 <- subset(all10.1.72, line==i)
  test <- t.test(sub0$ all0.1.72, sub10$all10.1.72)
  comp1.72[[i]] <- c(test$statistic, test$parameter, test$p.value)
}
FLC1.ttest.72 <- do.call("rbind", comp1.72)
FLC1.ttest.72 <- data.frame(line =1:20, FLC1.ttest.72)
names(FLC1.ttest.72)[4] <- "p.value"
FLC1.ttest.72$sig <- ifelse(FLC1.ttest.72$p.value < 0.05, "*", "")
print(FLC1.ttest.72)
FLC1.ttest.72$t <- round(FLC1.ttest.72$t, 2)
FLC1.ttest.72$df <- round(FLC1.ttest.72$df, 2)
FLC1.ttest.72$p.value <- round(FLC1.ttest.72$p.value, 4)
write.csv(FLC1.ttest.72, "manuscript/tables/TableS1b-FLC1ttest_72h.csv", row.names=FALSE)

##########################
#https://paperpile.com/view/adc154f8-62c6-08e2-9d84-eb0bbc7c7222
#The initial ratio was computed as log(wi_final/wmax_start) whereas the ratio of fitness gains was calculated as log(wi_final/wi_start) – log(wmax_final/wmax_start).

wmax_start <- max(fitFlow.ag$all0.1)
wmax_final <- max(fitFlow.ag$all10.1)
init <- 1+log10(fitFlow.ag$all0.1/wmax_start)/100
fitnessGain <- 1+(log10(fitFlow.ag$all10.1/fitFlow.ag$all0.1) - log10(wmax_final/wmax_start))/100

plot(init, fitnessGain)
summary(lm(fitnessGain~init))
#r-squared 0.7191

wmax_start72 <- max(fitFlow.ag$all0.1.72)
wmax_final72 <- max(fitFlow.ag$all10.1.72)
init72 <- 1+log(fitFlow.ag$all0.1.72/wmax_start72)/100
fitnessGain72 <- 1+(log(fitFlow.ag$all10.1.72/fitFlow.ag$all0.1.72) - log(wmax_final72/wmax_start72))/100
summary(lm(fitnessGain72~init72))
#r-squared 0.955
plot(init72, fitnessGain72)

##################
#whole strain set
##################
#The degree of fitness improvement was strongly influenced by initial strain fitness, while neither mating type nor clade had a significant effect (Figure 2, ANOVA tests with change in fitness values transformed using the log-modululs transformation to meet parametric assumptions; 24 h---initial fitness: F1, 12 = 15.55, p = 0.002, clade: F4, 12 = 1.48, p = 0.27, MAT zygosity: F1, 12 = 1.70, p = 0.22; 72 h---initial fitness: F1, 12 = 7.20, p = 0.020, clade: F4, 12 = 12.37, p = 0.16, MAT zygosity: F1, 11 = 0.04, p = 0.84).
#https://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg.html

#on the strain average
delta1.aov <- aov(delta1.24~ all0.1 +zygosity + as.factor(clade) , data=fitFlow.ag)
# Df Sum Sq Mean Sq F value Pr(>F)
# all0.ag.1         1 0.5410  0.5410   7.994 0.0143 *
# zygosity          1 0.0184  0.0184   0.271 0.6112
# as.factor(clade)  4 0.1734  0.0433   0.640 0.6430
# Residuals        13 0.8798  0.0677

#mixed-effect treating strain as a random effect - exact same result (woo!)
delta1.lmer <- lmer(delta1.24 ~ all0.1+ as.factor(clade) +zygosity + (1|strain), data = fitFlow)
anova(delta1.lmer, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)
# all0.1           0.74992 0.74992     1 25.173 11.2827 0.002494 **
# as.factor(clade) 0.17111 0.04278     4 13.646  0.6436 0.640561
# zygosity         0.03481 0.03481     1 13.672  0.5237 0.481477


delta1.72.aov <- aov(delta1.72~ all0.ag.1.72 +zygosity + clade , data=fitFlow.ag)
# Df Sum Sq Mean Sq F value       Pr(>F)
# all0.ag.1.72  1 1.3723  1.3723 158.730 0.0000000116 ***
# clade         4 0.0681  0.0170   1.970        0.159
# zygosity      1 0.0021  0.0021   0.248        0.626
# Residuals    13 0.1124  0.0086

#same result whether on average or looking at variability
delta1.72.lmer <- lmer(delta1.72 ~ all0.1.72+ as.factor(clade) +zygosity + (1|strain), data = fitFlow)
anova(delta1.72.lmer, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)
# all0.1.72        2.24616 2.24616     1 20.912 157.5214 3.355e-11 ***
# as.factor(clade) 0.11471 0.02868     4 13.523   2.0112    0.1500
# zygosity         0.00000 0.00000     1 13.466   0.0001    0.9941

###Look at variability
#can only do sd looking among replicates
delta1.variance24.aov <- aov(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1) ~ fitFlow.ag$all0.1 + fitFlow.sd$zygosity + as.factor(fitFlow.sd$clade))
# Df  Sum Sq Mean Sq F value    Pr(>F)
# fitFlow.ag$all0.1            1 0.24434 0.24434  40.941 0.0000235 ***
#   fitFlow.sd$zygosity          1 0.01325 0.01325   2.220     0.160
# as.factor(fitFlow.sd$clade)  4 0.01950 0.00488   0.817     0.537
# Residuals                   13 0.07759 0.00597

delta1.variance72.aov <- aov(c(fitFlow.sd$all10.1.72-fitFlow.sd$all0.1.72) ~ fitFlow.ag$all0.1.72 +zygosity + clade , data=fitFlow.ag)
# Df  Sum Sq Mean Sq F value Pr(>F)
# fitFlow.ag$all0.1.72  1 0.03437 0.03437   6.352 0.0256 *
#   zygosity              1 0.00614 0.00614   1.136 0.3060
# clade                 4 0.02195 0.00549   1.014 0.4356
# Residuals            13 0.07033 0.00541


# There is also a significant negative correlation between the degree of fitness change in low drug and the variance in evolved fitness among replicates – replicate lines were more variable in strain backgrounds that increased in fitness the most.
cor.test(fitFlow.ag$delta1.24, fitFlow.sd$all10.ag.1, paired=TRUE, method="spearman") #S = 694, p-value = 0.03451, rho = 0.48
cor.test(fitFlow.ag$delta1.72, fitFlow.sd$all10.ag.1.72, paired=TRUE, method="spearman") #S = 476, p-value = 0.00286, rho = 0.64

#################################
#Higher drug
#################################
# After ~100 generations of evolution at a low level of fluconazole, only 18 replicates from seven strain backgrounds increased in MIC beyond 1ug/mL measured through a broth microdilution assay (BMD). In most of these replicates, MIC increased only slightly, to 4 ug/mL (five replicates from strain A3, four replicates from A8, and one replicate from A4 and A7). Clinical resistance (MIC >= 8 at 24 hours, Fothergill 2014) was acquired in only eight replicates, six from strain A5 (which had an initial MIC of 4), and one replicate each from strains A17 and A19, two backgrounds with a very low initial MIC (Figure 1A).

MICto1 <- subset(fitFlow, MIC24 < 1 & MIC24.10 == 1) #97 of 132 reps increased to 1, 12 increased above 1
MICinc1 <- subset(fitFlow, MIC24 <= 1 & MIC24.10 > 1)
MICinc8 <- subset(fitFlow, MIC24 <= 8 & MIC24.10 >= 8)

nrow(subset(fitFlow, MIC24 == 1 & MIC24.10 < 1)) #4
nrow(subset(fitFlow, MIC24 > 1 & MIC24.10 == 1)) #40

 # Interestingly, the MIC of some replicates from initially-high MIC strains actually decreased to 1ug/mL  (two replicates from A5, eight replicates from A12, all 12 replicates from A20).

MICdec <- subset(fitFlow, change24 < 0)

#SMG all replicates - this is using the range of ancestral values as the limits for up/down in evolved relicates
SMG24up <- subset(fitFlow, SMG24.10.2 > rep(SMG24.up.2, each=12)) #110
SMG24down <- subset(fitFlow, SMG24.10.2 < rep(SMG24.down.2, each=12)) #52
SMG48up <- subset(fitFlow, SMG48.10.2 > rep(SMG48.up.2, each=12)) #68
SMG48down <- subset(fitFlow, SMG48.10.2 < rep(SMG48.down.2, each=12)) #91
SMG72up <- subset(fitFlow, SMG72.10.2 > rep(SMG72.up.2, each=12)) #74
SMG72down <- subset(fitFlow, SMG72.10.2 < rep(SMG72.down.2, each=12)) #88

SMG72up.MIC1 <- subset(fitFlow, SMG72.10.2 > rep(SMG72.up.2, each=12) & MIC24.10 == 1) #42 in  14 strains (3, 12 unique)
SMG72down.MIC1 <- subset(fitFlow, SMG72.10.2 < rep(SMG72.down.2, each=12) & MIC24.10 == 1) #70 in 15  strains (14, 16, 8 unique)
# A1 & A18 had no replicates MIC=1, A5 started with very high tolerance and did not change

table(SMG72up.MIC1$strain)
# 10 11 12 13 15 17 19  2 20  3  4  6  7  9
# 3  2  6  3  1  1  2  2  5  5  3  3  5  1

table(SMG72down.MIC1$strain)
# 10 11 13 14 15 16 17 19  2 20  4  6  7  8  9
# 3  1  6  9  8  9  6  1  4  3  7  3  1  5  4

cor.test(fitFlow.ag$SMG72, fitFlow.ag$MIC24, method="spearman")
S = 1068.8, p-value = 0.4066

fitFlow_MIC1 <- subset(fitFlow, MIC24.10 == 1)
fitFlow_MIC1$clade <- as.factor(fitFlow_MIC1$clade)
fitFlow_MIC1.ag <- aggregate(fitFlow_MIC1[c("change72.SMG.2", "delta1.72", "all10.1.72", "all0.1.72")], fitFlow_MIC1[c("line", "clade", "zygosity", "col72")], mean)

fitFlow_MIC1.sd <- aggregate(fitFlow_MIC1[c("change72.SMG.2", "delta1.72", "all10.1.72", "all0.1.72")], fitFlow_MIC1[c("line", "clade", "zygosity", "col72")], sd)


fitFlow_MIC1.aov <- aov(change72.SMG.2 ~ delta1.72 + as.factor(clade) +as.factor(zygosity), data = fitFlow_MIC1.ag)
summary(fitFlow_MIC1.aov)
# Df  Sum Sq Mean Sq F value Pr(>F)
# delta1.72            1 0.00199 0.00199   0.074 0.7910
# as.factor(clade)     4 0.21955 0.05489   2.031 0.1593
# as.factor(zygosity)  1 0.15925 0.15925   5.894 0.0335 *
#   Residuals           11 0.29722 0.02702
beeswarm(fitFlow_MIC1.ag$change72.SMG.2~as.factor(fitFlow_MIC1.ag$zygosity))
fitFlow_MIC1.ag
#really just driven by two strains (3 and 12) being high and homozygous and one strain (14) being low and heterozygous

plot(fitFlow_MIC1.ag$delta1.72, fitFlow_MIC1.ag$change72.SMG.2)
cor.test(fitFlow_MIC1.ag$delta1.72, fitFlow_MIC1.ag$change72.SMG.2)

fitFlow_MIC1.aov.sd <- aov(change72.SMG.2 ~ delta1.72 + as.factor(clade) +as.factor(zygosity), data = fitFlow_MIC1.sd)
summary(fitFlow_MIC1.aov.sd)
# Df  Sum Sq  Mean Sq F value Pr(>F)
# delta1.72            1 0.01726 0.017264   2.365  0.152
# as.factor(clade)     4 0.04268 0.010671   1.462  0.279
# as.factor(zygosity)  1 0.00404 0.004041   0.554  0.472
# Residuals           11 0.08029 0.007299

#none of this is significant.
cor.test(fitFlow_MIC1.sd$delta1.72, fitFlow_MIC1.sd$change72.SMG.2)
cor.test(fitFlow_MIC1.ag$delta1.72, fitFlow_MIC1.sd$change72.SMG.2)
cor.test(fitFlow_MIC1.sd$all10.1.72, fitFlow_MIC1.sd$change72.SMG.2)
cor.test(fitFlow_MIC1.ag$all10.1.72, fitFlow_MIC1.sd$change72.SMG.2)



#HERE
hist(fitFlow_MIC1$change72.SMG.2)
quantile(fitFlow_MIC1$change72.SMG.2)
# 0%         25%         50%         75%        100%
# -0.66702455 -0.11654942 -0.01970274  0.06392017  0.63506899
highSMG.10 <- subset(fitFlow_MIC1, change72.SMG.2 > quantile(fitFlow_MIC1$change72.SMG.2)[4])
table(highSMG.10$strain)
# 11 12 13 15 17 19 20  3  4  6  7  8  9
# 2  6  3  1  1  3  6  5  3  3  5  1  1
lowSMG.10 <- subset(fitFlow_MIC1, change72.SMG.2 < quantile(fitFlow_MIC1$change72.SMG.2)[2])
table(lowSMG.10$strain)
# 10 13 14 15 16 17 19 20  4  7  8  9
# 3  5  8  7  8  1  1  2  1  1  2  1

# Strain 13 has 3 of the highest SMG replicates and 5 of the lowest
# Still not enough variation in fitness to be significant (although tempting)
fitfFlow_MIC1_13 <- subset(fitFlow_MIC1, strain == 13)
plot(fitfFlow_MIC1_13$all10.1.72, fitfFlow_MIC1_13$change72.SMG.2)
abline(lm(fitfFlow_MIC1_13$change72.SMG.2~fitfFlow_MIC1_13$all10.1.72))
cor.test(fitfFlow_MIC1_13$change72.SMG.2,fitfFlow_MIC1_13$all10.1.72)

# Strain 20 hass 6 of the highest SMG replicates and 2 of the lowest
# Really not enough variation in fitness to say anything (and can easily see they are well spread out)
fitfFlow_MIC1_20 <- subset(fitFlow_MIC1, strain == 20)
plot(fitfFlow_MIC1_20$all10.1.72, fitfFlow_MIC1_20$change72.SMG.2)
abline(lm(fitfFlow_MIC1_20$change72.SMG.2~fitfFlow_MIC1_20$all10.1.72))
cor.test(fitfFlow_MIC1_20$change72.SMG.2,fitfFlow_MIC1_20$all10.1.72)

#Full model
SMGevol.lmer <- lmer(change72.SMG.2~ delta1.72 + zygosity+as.factor(clade)+ (1|strain), data=fitFlow_MIC1, REML = FALSE)
anova(SMGevol.lmer, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)
# delta1.72        0.14029 0.14030     1 61.900  5.4427 0.022916 *
# zygosity         0.21865 0.21865     1 15.746  8.4823 0.010306 *
# as.factor(clade) 0.53585 0.13396     4 16.018  5.1970 0.007055 **
summary(SMGevol.lmer)

qqnorm(resid(SMGevol.lmer))
qqline(resid(SMGevol.lmer))

#Not influenced by a single strain
for(i in unique(fitFlow_MIC1$strain)){
  print(i)
  test_lmer_drop <- lmer(change72.SMG.2~ delta1.72 + zygosity+as.factor(clade) + (1|strain), data=subset(fitFlow_MIC1, strain!=i, REML = FALSE))
  print(anova(test_lmer_drop, type = 3))
}

#fitness: A10, A13, A14
#zygosity: A3, A5, A6, A7, A8, A9, A12, A13, A14, A15, A16, A17, A20
#clade: A4, A5, A6, A7, A8, A9, A10, A12, A13, A14, A15, A17, A20

cor.test(fitFlow.sd$change72.SMG.2, fitFlow.ag$all0.1, method="spearman") #S = 2076, p-value = 0.01131, rho = -0.55
cor.test(fitFlow.sd$change72.SMG.2, fitFlow.ag$all0.1.72, method="spearman") #S = 1838, p-value = 0.09739

cor.test(fitFlow.sd$change72.SMG.2, fitFlow.sd$delta1.24, method="spearman") #S = 586, p-value = 0.01157
cor.test(fitFlow.sd$change72.SMG.2, fitFlow.sd$delta1.72, method="spearman") #S = 556, p-value = 0.008166



#is there really a point at looking at the aggreagate among strain data?
fit0_flow.variable <- cor.test(fitFlow.ag.variable$all0.1, fitFlow.ag.variable$t10.G1.1) #t = -0.20031, df = 13, p-value = 0.8443
fit10.72_flow.variable <- cor.test(fitFlow.ag.variable$all10.1.72, fitFlow.ag.variable$t10.G1.1) #t = -0.21898, df = 13, p-value = 0.8301

fit_flow.variable.aov <- lmer(t10.G1.1 ~ delta1.72 + zygosity + as.factor(clade) + (1 | strain), data=fitFlow.variable, na.action=na.omit)
anova(fit_flow.variable.aov)
# Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# delta1.72           21.3    21.3     1 46.408  0.0048 0.9450
# zygosity             7.3     7.3     1  8.688  0.0016 0.9686
# as.factor(clade) 16659.1  4164.8     4  9.206  0.9402 0.4825


fitD_flow.variable24 <- cor.test(fitFlow.ag.variable$change24.SMG.2, fitFlow.ag.variable$t10.G1.1, method="spearman") #S = 700.13, p-value = 0.3684
fitD_flow.variable72 <- cor.test(fitFlow.ag.variable$change72.SMG.2, fitFlow.ag.variable$t10.G1.1, method="spearman") #S = 1027.3, p-value = 0.04325, rho = -0.51
plot(fitFlow.ag.variable$change72.SMG.2, fitFlow.ag.variable$t10.G1.1)


# Interestingly, there was a significant correlation between the variance of evolved low drug fitness (measured at 24 h) and variance of evolved genome size (Pearson's correlation, t14 = 5.65, p-value = 0.00002, cor = 0.80, this remains significant when the four non-variable lines are removed: t14 = 2.31, p-value = 0.037, cor = 0.51).

shapiro.test(fitFlow.ag$t10.G1.1) #W = 0.95807, p-value = 0.5061
qqnorm(fitFlow.ag$t10.G1.1)
shapiro.test(fitFlow.sd$all10.ag.1) #W = 0.93818, p-value = 0.2214
shapiro.test(fitFlow.sd.variable$all10.ag.1) #W = 0.97123, p-value = 0.8581
shapiro.test(fitFlow.sd.sub$FoG20.10)
qqnorm(fitFlow.sd.sub$FoG20.10)

fit10_flowCV <- cor.test(fitFlow.sd$all10.1, fitFlow.sd$t10.G1.1, method="spearman") #S = 384, p-value = 0.0006306, rho = 0.71
fit10_flowCV <- cor.test(fitFlow.sd$all10.1.72, fitFlow.sd$t10.G1.1, method="spearman") #S = 554, p-value = 0.007973, rho = 0.0080
fit10_flowCV.variable <- cor.test(fitFlow.sd.variable$all10.ag.1, fitFlow.sd.variable$t10.G1.1) #t = 2.3115, df = 14, p-value = 0.03655

cor.test(fitFlow.sd$t10.G1.1, fitFlow.ag$all0.1, method="spearman") #S = 2076, p-value = 0.01131, rho = -0.56
cor.test(fitFlow.sd$t10.G1.1, fitFlow.ag$all0.1.72, method="spearman") #S = 1908, p-value = 0.05692, rho = -0.43

cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$delta1.24, method="spearman") #S = 362, p-value = 0.0004083, rho = 0.73
cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$delta1.72, method="spearman") #S = 440, p-value = 0.00166, rho = 0.67

cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$change72.SMG.2, method="spearman") #S = 708, p-value = 0.039, rho = 0.47

 #################################
 #Fitness in YPD
 #################################
 # When fitness was assessed in a permissive environment (standard lab YPD, the base medium for the evolutionary environment), there was a significant reduction in fitness for evolved strains compared to ancestral most strain backgrounds at 24 h

 comp0 <- list()
 stat <- c()
 for(i in 1:20){
   sub0 <- subset(all0.0, line==i)
   s0 <- shapiro.test(sub0$data)
   sub10 <- subset(all10.0, line==i)
   s10 <- shapiro.test(sub10$data)
   if(s0$p.value > 0.05 & s10$p.value > 0.05){
     v <- var.test(sub0$data, sub10$data)
     if(v$p.value > 0.05){
       test <- t.test(sub0$data, sub10$data, var.equal=TRUE)
       stat <- append(stat, "t-test")
     }
     else{
       test <- t.test(sub0$data, sub10$data, var.equal=FALSE)
       stat <- append(stat, "t-test-welch")
     }
   }
   else{
     test <- wilcox.test(sub0$data, sub10$data)
     test$parameter <- NA
     stat <- append(stat, "wilcox")
   }
   comp0[[i]] <- c(statistic = test$statistic, df = test$parameter, "p.value"= test$p.value, "imp" = ifelse(mean(sub10$data) > mean(sub0$data), "evol", "anc"))
 }
 FLC0.ttest <- do.call("rbind", comp0)
 FLC0.ttest <- data.frame(line =1:20, "test" = stat, FLC0.ttest)
 names(FLC0.ttest) <- c("line", "test", "statistic", "df", "p.value", "imp")
 FLC0.ttest$sig <- ifelse(FLC0.ttest$p.value < 0.05, "*", "")
 FLC0.ttest$statistic<- round(as.numeric(FLC0.ttest$statistic), 2)
 FLC0.ttest$df <- round(as.numeric(FLC0.ttest$df), 2)
 FLC0.ttest$p.value <- round(as.numeric(FLC0.ttest$p.value), 4)
 write.csv(FLC0.ttest, "/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/tables/170802FLC0ttest.csv", row.names=FALSE)


 comp0.72 <- list()
 stat <- c()
 for(i in 1:20){
   sub0 <- subset(all0.0.72, line==i)
   s0 <- shapiro.test(sub0$data)
   sub10 <- subset(all10.0.72, line==i)
   s10 <- shapiro.test(sub10$data)
   if(s0$p.value > 0.05 & s10$p.value > 0.05){
     v <- var.test(sub0$data, sub10$data)
     if(v$p.value > 0.05){
       test <- t.test(sub0$data, sub10$data, var.equal=TRUE)
       stat <- append(stat, "t-test")
     }
     else{
       test <- t.test(sub0$data, sub10$data, var.equal=FALSE)
       stat <- append(stat, "t-test-welch")
     }
   }
   else{
     test <- wilcox.test(sub0$data, sub10$data)
     test$parameter <- NA
     stat <- append(stat, "wilcox")
   }
   comp0.72[[i]] <- c(statistic = test$statistic, df = test$parameter, "p.value"= test$p.value, "imp" = ifelse(mean(sub10$data) > mean(sub0$data), "evol", "anc"))
 }
 FLC0.ttest.72 <- do.call("rbind", comp0.72)
 FLC0.ttest.72 <- data.frame(line =1:20, "test" = stat, FLC0.ttest.72)
 names(FLC0.ttest.72) <- c("line", "test", "statistic", "df", "p.value", "imp")
 FLC0.ttest.72$sig <- ifelse(FLC0.ttest.72$p.value < 0.05, "*", "")
 FLC0.ttest.72$statistic<- round(as.numeric(FLC0.ttest.72$statistic), 2)
 FLC0.ttest.72$df <- round(as.numeric(FLC0.ttest.72$df), 2)
 FLC0.ttest.72$p.value <- round(as.numeric(FLC0.ttest.72$p.value), 4)
 write.csv(FLC0.ttest.72, "/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/tables/170802FLC0ttest-72h.csv", row.names=FALSE)

 test0 <-data.frame(line =1:20, YPD=FLC0.ttest$sig, YPD.imp = FLC0.ttest$imp, YPD.72=FLC0.ttest.72$sig, YPD.imp.72 = FLC0.ttest.72$imp)
 test0

 cor.test(fitFlow.ag$delta1.24, fitFlow.sd$all10.ag.1, paired=TRUE, method="spearman") #S = 694, p-value = 0.03451, rho = 0.48

 # For our low drug evolved strains, however, we recovered a significant positive relationship among strain backgrounds when we compare mean fitness assessed at 24 h at low drug and the permissive environment (spearman’s rank correlation test; S = 396, p = 0.0008, rho = 0.702; there is not enough variation in the permissive environment at 72 h for correlative purposes)
 cor.test(all10.ag.0$data, fitFlow.ag$all10.ag.1, method="spearman")

 parameter <- c()
 p.value <- c()
 rho <- c()
 for(i in 1:20){
   test <- cor.test(subset(all10.0, line==i)$data, subset(all10.1, line==i)$data, method="spearman")
   parameter[i] <- test$statistic
   p.value[i] <- test$p.value
   rho[i] <- test$estimate
 }

 ypdCor <- data.frame(line = 1:20, s =parameter, p.value, rho)
 write.csv(ypdCor, "/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/tables/171030YPD-FLC1_cor-24h.csv", row.names=FALSE)

