##########################################################
### Simulations for validating the CladeDate algorithm ###
##########################################################


library(ape)

library(CladeDate)



setwd("DIR") # set the directory where simulated dataset reside



#################
### MAIN LOOP ###


REPS <- 1:100

# matrix for storing results

BDages <- matrix(ncol=9, nrow=max(REPS))

CHages <- matrix(ncol=9, nrow=max(REPS))

# Create a big loop that repeats the script 

for(i in REPS) {


### Load trees and calibrations ###

# read B-D and ML tree

tr <- read.tree(file=paste0("Simulated",i,".tre"))

MLtr <- read.tree(file=paste0("EstimatedML",i,".tre"))

# load fossil record

load(file=paste0("FossilRecordClade1.",i,".R")) # fr.clade1
load(file=paste0("FossilRecordClade2.",i,".R")) # fr.clade2


### Estimate tehe age of clade from the fossil record ###

## Clade 1 

# first estiamte using the Strauss-Sadler method

date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, method="StraussSadler", PDFfitting="lognormal", KStest=TRUE)

# But if fossil record is not uniformly distributed, use the Robson-Whitlock method
if(date.clade1$KStest$p.value < 0.2) {
date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, method="RobsonWhitlock", PDFfitting="lognormal") }

## Clade 2

# first estiamte using the Strauss-Sadler method

date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, method="StraussSadler", PDFfitting="lognormal", KStest=TRUE)

# But if fossil record is not uniformly distributed, use the Robson-Whitlock method

if(date.clade2$KStest$p.value < 0.2) {
date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, method="RobsonWhitlock", PDFfitting="lognormal") }



### Estimate with Chronos ###

# Point calibration using the median 

Calib <- makeChronosCalib( phy = MLtr,
node = c(fr.clade1$mrca.node, fr.clade2$mrca.node),
age.min = c(date.clade1$Quantiles[2], date.clade2$Quantiles[2]) )

# Multiple branch lengths by sequence length

MLtr$edge.length <- MLtr$edge.length*1000

Chrono <- chronos( phy = MLtr, model="clock", calibration=Calib)

cat(attr(Chrono, "message"))

# If false convergence, don't write results and jumpt to next replicate
if(attr(Chrono, "message") == "false convergence (8)") {
	cat(paste("\nFalse convergence in replicate",i,": results skept\n"))
	next
	}


BDages[i,] <- branching.times(tr)

CHages[i,] <- branching.times(Chrono)

cat(paste("\nReplicate",i,"completed\n"))

}

### END OF LOOP ###
###################




#########################################
### Calculating Acuracy and Precision ###


### Regression Method ###
# Slope (acuracy, bias)
# Coefficient of Determination (precision)

M1 <- lm(estages ~ trueages -1)

summary(M1)

Call:
lm(formula = estages ~ trueages - 1)

Residuals:
     Min       1Q   Median       3Q      Max 
-12.1600  -0.8985   0.0312   0.7716  16.9952 

Coefficients:
         Estimate Std. Error t value Pr(>|t|)    
trueages 0.970067   0.004359   222.5   <2e-16 ***
---
Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 2.444 on 755 degrees of freedom
Multiple R-squared:  0.985,	Adjusted R-squared:  0.985 
F-statistic: 4.952e+04 on 1 and 755 DF,  p-value: < 2.2e-16

confint(M1, "trueages")
             2.5 %    97.5 %
trueages 0.9615097 0.9786251

# Presicion

cor(estages, trueages)

# 0.9735038

########################################




#############
### PLOT  ###
#############

### Simulated vs. estiamted ###

Max <- max(c(CHages, BDages), na.rm=TRUE)

col.cd <- rgb( .8, .3, 0)
col.cd <- rgb( .7, .2, 0)
col.cd.bg <- rgb( .8, .4, 0, 0.5)

quartz("trees",5,5); par(mgp=c(1.8,0.4,0), tcl=-0.3, las=1, cex.axis=0.8)

plot(0:24, 0:24, type="n", xlim=c(0, Max), ylim=c(0, Max), xlab="True node age", ylab="Estimated node age", log="")

points(BDages, CHages, col= col.cd.bg, bg= col.cd.bg, pch=21, cex=0.8) 

abline(a=0, b=1)
abline(M1, lty="dashed")

###########
### END ###
###########
