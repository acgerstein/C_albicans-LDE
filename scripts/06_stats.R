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
delta1.lm <- lm(delta1.24~ all0.1 +zygosity + as.factor(clade) , data=fitFlow.ag)
d1 <- tidy(anova(delta1.lm))
d1$varExp <- d1[,3]/sum(d1[,3])*100
# term                df  sumsq meansq statistic p.value varExp$sumsq
# <chr>            <int>  <dbl>  <dbl>     <dbl>   <dbl>        <dbl>
#   1 all0.1               1 0.541  0.541      7.99   0.0143        33.5
# 2 zygosity             1 0.0184 0.0184     0.271  0.611          1.14
# 3 as.factor(clade)     4 0.173  0.0433     0.640  0.643         10.8
# 4 Residuals           13 0.880  0.0677    NA     NA             54.6


#mixed-effect treating strain as a random effect - exact same result (woo!)
#not in MS
delta1.lmer <- lmer(delta1.24 ~ all0.1+ as.factor(clade) +zygosity + (1|strain), data = fitFlow)
anova(delta1.lmer, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)
# all0.1           0.74992 0.74992     1 25.173 11.2827 0.002494 **
# as.factor(clade) 0.17111 0.04278     4 13.646  0.6436 0.640561
# zygosity         0.03481 0.03481     1 13.672  0.5237 0.481477

delta1.72.aov <- lm(delta1.72~ all0.1.72 +zygosity + clade , data=fitFlow.ag)
d1.72 <- tidy(delta1.72.aov)
d1.72$varExp <- d1.72[,3]/sum(d1.72[,3])*100
# term         df   sumsq  meansq statistic       p.value varExp$sumsq
# <chr>     <dbl>   <dbl>   <dbl>     <dbl>         <dbl>        <dbl>
# all0.1.72     1 1.37    1.37      159.     0.0000000116       88.3
# zygosity      1 0.00215 0.00215     0.248  0.626               0.138
# clade         4 0.0681  0.0170      1.97   0.159               4.38
# Residuals    13 0.112   0.00865    NA     NA                   7.23



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
delta1.variance24.lm <- lm(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1) ~ fitFlow.ag$all0.1 + fitFlow.sd$zygosity + as.factor(fitFlow.sd$clade))
d1.sd <- tidy(anova(delta1.variance24.lm))
d1.sd$varExp <- d1.sd[,3]/sum(d1.sd[,3])*100
# term                           df  sumsq  meansq statistic    p.value varExp$sumsq
# <chr>                       <int>  <dbl>   <dbl>     <dbl>      <dbl>        <dbl>
#   1 fitFlow.ag$all0.1               1 0.244  0.244      40.9    0.0000235        68.9
# 2 fitFlow.sd$zygosity             1 0.0132 0.0132      2.22   0.160             3.74
# 3 as.factor(fitFlow.sd$clade)     4 0.0195 0.00488     0.817  0.537             5.50
# 4 Residuals                      13 0.0776 0.00597    NA     NA                21.9


delta1.variance72.lm <- lm(c(fitFlow.sd$all10.1.72-fitFlow.sd$all0.1.72) ~ fitFlow.ag$all0.1.72 +zygosity + clade , data=fitFlow.ag)
d1.sd.72 <- tidy(anova(delta1.variance72.lm))
d1.sd.72$varExp <- d1.sd.72[,3]/sum(d1.sd.72[,3])*100
# term                    df   sumsq  meansq statistic p.value varExp$sumsq
# <chr>                <int>   <dbl>   <dbl>     <dbl>   <dbl>        <dbl>
#   1 fitFlow.ag$all0.1.72     1 0.0344  0.0344       6.35  0.0256        25.9
# 2 zygosity                 1 0.00614 0.00614      1.14  0.306          4.63
# 3 clade                    4 0.0220  0.00549      1.01  0.436         16.5
# 4 Residuals               13 0.0703  0.00541     NA    NA             53.0

# There is also a significant negative correlation between the degree of fitness change in low drug and the variance in evolved fitness among replicates – replicate lines were more variable in strain backgrounds that increased in fitness the most.
cor.test(fitFlow.ag$delta1.24, fitFlow.sd$all10.ag.1, paired=TRUE, method="spearman") #S = 694, p-value = 0.03451, rho = 0.48
cor.test(fitFlow.ag$delta1.72, fitFlow.sd$all10.ag.1.72, paired=TRUE, method="spearman") #S = 476, p-value = 0.00286, rho = 0.64

#################################
#Higher drug
#################################
# After ~100 generations of evolution at a low level of fluconazole, only 18 replicates from seven strain backgrounds increased in MIC beyond 1ug/mL measured through a broth microdilution assay (BMD). In most of these replicates, MIC increased only slightly, to 4 ug/mL (five replicates from strain A3, four replicates from A8, and one replicate from A4 and A7). Clinical resistance (MIC >= 8 at 24 hours, Fothergill 2014) was acquired in only eight replicates, six from strain A5 (which had an initial MIC of 4), and one replicate each from strains A17 and A19, two backgrounds with a very low initial MIC (Figure 1A).

nrow(subset(fitFlow, MIC24 < 1)) #144
nrow(subset(fitFlow, MIC24 < 1 & MIC24.10 == 1)) #110 of 144 reps increased to 1,
nrow(subset(fitFlow, MIC24 < 1 & MIC24.10 > 1)) #6 increased above 1
nrow(subset(fitFlow, MIC24 == 1 & MIC24.10 > 1)) #0
nrow(subset(fitFlow, MIC24 == 1 & MIC24.10 < 1)) #0

nrow(subset(fitFlow, MIC24 > 1 & MIC24.10 == 1)) #16  A5, A20
nrow(subset(fitFlow, MIC24 > 1 & change24 > 1)) #4 A5

subset(fitFlow, change24 > 1)

 # Interestingly, the MIC of some replicates from initially-high MIC strains actually decreased to 1ug/mL  (two replicates from A5, all 12 replicates from A20).


nrow(subset(fitFlow, change24 > 1)) #32
MICdec <- subset(fitFlow, change24 < 0)

#SMG all replicates - this is using the range of ancestral values as the limits for up/down in evolved relicates
SMG24up <- subset(fitFlow, SMG24.10 > rep(SMG24.up, each=12)) #79
SMG24down <- subset(fitFlow, SMG24.10 < rep(SMG24.down, each=12)) #66
SMG48up <- subset(fitFlow, SMG48.10 > rep(SMG48.up, each=12)) #47
SMG48down <- subset(fitFlow, SMG48.10 < rep(SMG48.down, each=12)) #118
SMG72up <- subset(fitFlow, SMG72.10 > rep(SMG72.up, each=12)) #55
SMG72down <- subset(fitFlow, SMG72.10 < rep(SMG72.down, each=12)) #112

SMG72up.strains <- subset(fitFlow, SMG72.10 > rep(SMG72.up, each=12)) #55
length(unique(SMG72up.strains$strain)) #15
SMG72down.strains <- subset(fitFlow, SMG72.10 < rep(SMG72.down, each=12)) #122
length(unique(SMG72down.strains$strain)) #19

# SMG72up.MIC1 <- subset(fitFlow, SMG72.10 > rep(SMG72.up, each=12) & MIC24.10 == 1) #29 in  14 strains (3, 12 unique)
# SMG72down.MIC1 <- subset(fitFlow, SMG72.10 < rep(SMG72.down, each=12) & MIC24.10 == 1) #93 in 15  strains (14, 16, 8 unique)
# # A1, A12 & A18 had no replicates MIC=1, A5 started with very high tolerance and did not change
#
# table(SMG72up.MIC1$strain)
# # 10 11 13 15 17 19  2 20  3  4  6  7  8  9
# # 2  2  1  1  1  2  1  4  5  3  1  4  1  1
#
# table(SMG72down.MIC1$strain)
# # 10 11 13 14 15 16 17 19  2 20  4  6  7  8  9
# # 4  1  9 10  8 12  9  4  5  5  7  7  4  4  4

cor.test(fitFlow.ag$SMG72, fitFlow.ag$MIC24, method="spearman")
#S = 1141.8, p-value = 0.5518
cor.test(fitFlow.ag$SMG72.10, fitFlow.ag$MIC24.10, method="spearman")
#S = 1205.2, p-value = 0.694
cor.test(fitFlow.sd$SMG72.10, fitFlow.ag$all0.1, method="spearman")
#S = 2226, p-value = 0.001508
#rho = -0.67
cor.test(fitFlow.sd$SMG72.10, fitFlow.ag$all0.1.72, method="spearman")
#S = 1986, p-value = 0.02866
# rho = -0.49
cor.test(fitFlow.sd$SMG72.10, fitFlow.sd$all10.1, method="spearman")
#S = 450, p-value = 0.001941, rho = 0.66
cor.test(fitFlow.sd$SMG72.10, fitFlow.sd$all10.1.72, method="spearman")
#S = 658, p-value = 0.02458, rho = 0.51

plot(fitFlow.sd$SMG72.10, fitFlow.ag$all0.1)
plot(fitFlow.sd$SMG72.10, fitFlow.ag$all0.1.72)
cor.test(fitFlow_MIC1.sd$SMG72.10, fitFlow_MIC1.ag$all0.1)
plot(fitFlow_MIC1.sd$SMG72.10, fitFlow_MIC1.ag$all0.1.72)

# fitFlow_MIC1 <- subset(fitFlow, MIC24.10 == 1)
# fitFlow_MIC1$clade <- as.factor(fitFlow_MIC1$clade)
# fitFlow_MIC1.ag <- aggregate(fitFlow_MIC1[c("MIC24", "MIC24.10", "all0.1","all0.1.72", "all10.1", "all10.1.72", "delta1.72", "SMG72", "SMG72.10", "change72.SMG")], fitFlow_MIC1[c("line", "clade", "zygosity", "col72")], mean)
# fitFlow_MIC1.ag <- fitFlow_MIC1.ag[order(as.numeric(fitFlow_MIC1.ag$line)),]
#
# fitFlow_MIC1.sd <- aggregate(fitFlow_MIC1[c("MIC24","MIC24.10", "all0.1","all0.1.72", "all10.1", "all10.1.72", "delta1.72", "SMG72", "SMG72.10", "change72.SMG")], fitFlow_MIC1[c("line", "clade", "zygosity", "col72")], sd)
# fitFlow_MIC1.sd <- fitFlow_MIC1.sd[order(as.numeric(fitFlow_MIC1.sd$line)),]
#
# cor.test(fitFlow_MIC1.ag$SMG72, fitFlow_MIC1.ag$MIC24, method="spearman")
# #S = 796.88, p-value = 0.9289
# cor.test(fitFlow_MIC1.ag$SMG72.10, fitFlow_MIC1.ag$MIC24.10, method="spearman")
# # MIC are all 1. Duh,
# cor.test(fitFlow_MIC1.sd$SMG72.10, fitFlow_MIC1.ag$all0.1, method="spearman")
# #S = 1192, p-value = 0.06454
# #rho = -0.46
# cor.test(fitFlow_MIC1.sd$SMG72.10, fitFlow_MIC1.ag$all0.1.72, method="spearman")
# #S = 996, p-value = 0.3934
# # rho = -0.22
#
#
# cor.test(fitFlow_MIC1.sd$SMG72.10, fitFlow_MIC1.sd$all10.1, method="spearman")
# #S = 738, p-value = 0.7155
# cor.test(fitFlow_MIC1.sd$SMG72.10, fitFlow_MIC1.sd$all10.1.72, method="spearman")
# #S = 546, p-value = 0.1944

#Full model
# SMGevol.lmer <- lmer(change72.SMG~ delta1.72 + zygosity+as.factor(clade)+ (1|strain), data=fitFlow_MIC1, REML = FALSE)
# anova(SMGevol.lmer, type = 3)
# # Type III Analysis of Variance Table with Satterthwaite's method
# #                   Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)
# # delta1.72        0.07038 0.070382     1 48.850  2.6657 0.10897
# # zygosity         0.14625 0.146255     1 12.477  5.5393 0.03574 *
# # as.factor(clade) 0.46026 0.115065     4 12.477  4.3580 0.01984 *
# summary(SMGevol.lmer)

SMGevol.lmer_all <- lmer(change72.SMG~ delta1.72 + zygosity+as.factor(clade)+ (1|strain), data=fitFlow, REML = FALSE)
anova(SMGevol.lmer_all, type = 3)
# Type III Analysis of Variance Table with Satterthwaite's method
#                    Sum Sq  Mean Sq NumDF   DenDF F value Pr(>F)
# delta1.72        0.000710 0.000710     1 132.478  0.0269 0.8699
# zygosity         0.017113 0.017113     1  18.774  0.6494 0.4304
# as.factor(clade) 0.209510 0.052378     4  19.051  1.9877 0.1372


#####FLOW

fit10_flowCV <- cor.test(fitFlow.sd$all10.1, fitFlow.sd$t10.G1.1, method="spearman") #S = 384, p-value = 0.0006306, rho = 0.71
fit10_flowCV <- cor.test(fitFlow.sd$all10.1.72, fitFlow.sd$t10.G1.1, method="spearman") #S = 554, p-value = 0.007973, rho = 0.0080
fit10_flowCV.variable <- cor.test(fitFlow.sd.variable$all10.ag.1, fitFlow.sd.variable$t10.G1.1) #t = 2.3115, df = 14, p-value = 0.03655

cor.test(fitFlow.sd$t10.G1.1, fitFlow.ag$all0.1, method="spearman") #S = 2076, p-value = 0.01131, rho = -0.56
cor.test(fitFlow.sd$t10.G1.1, fitFlow.ag$all0.1.72, method="spearman") #S = 1908, p-value = 0.05692, rho = -0.43

cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$delta1.24, method="spearman") #S = 362, p-value = 0.0004083, rho = 0.73
cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$delta1.72, method="spearman") #S = 440, p-value = 0.00166, rho = 0.67

cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$SMG72.10, method="spearman")
#S = 538, p-value = 0.006555, rho = 0.60

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

