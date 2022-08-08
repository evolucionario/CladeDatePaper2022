##########################################################
### Simulations for validating the CladeDate algorithm ###
##########################################################

# 1 - teste te ability to estimate the correct date: measure acuracy and presicion

# 2 - Compare with other method:
# 	a. calculate an average relative age (clade age/oldest fossil age) and apply to both calibrations.
#	b. compare with CladeAge

library(ape)


setwd("~/Documents/Avian Time Trees/Simulations Chronos/")
setwd('~/Google Drive/BigBirdTree Calibration/Simulations Chronos/')
setwd('~/Google Drive/BigBirdTree Calibration/Simulations/')


REPS <- 1:50


# matrix for storing results

BDages <- matrix(ncol=9, nrow=max(REPS))

CHages <- matrix(ncol=9, nrow=max(REPS))



#################
### MAIN LOOP ###
#################

# Create a big loop that repeats the script 

for(i in REPS) {


###################################
### Load trees and calibrations ###
###################################

# read B-D and ML tree

tr <- read.tree(file=paste0("Simulated",i,".tre"))

MLtr <- read.tree(file=paste0("EstimatedML",i,".tre"))

# load fossil record

load(file=paste0("FossilRecordClade1.",i,".R")) # fr.clade1
load(file=paste0("FossilRecordClade2.",i,".R")) # fr.clade2

# load calibration information

load(file=paste0("Clade1date.",i,".R"))
load(file=paste0("Clade2date.",i,".R"))


#############################
### Estimate with Chronos ###
#############################

Calib <- makeChronosCalib( phy = MLtr,
node = c(fr.clade1$mrca.node, fr.clade2$mrca.node),
age.min = c(date.clade1$Quantiles[2], date.clade2$Quantiles[2]) )

Chrono <- chronos( phy = MLtr, model="clock", calibration=Calib)

save(Chrono, file=paste0("Chronogram.",i,".R"))

BDages[i,] <- branching.times(tr)

CHages[i,] <- branching.times(Chrono)

paste("Replicate ",i," finished")

}

### END OF MAIN LOOP ###






#########################
### COMPARING RESULTS ###
#########################

quartz("Fig", 6, 3.5)
par(mfcol=c(1,2), mar=c(2,0,1,0), oma=c(2,5,0,1), las=1, cex.axis=0.8, tcl=-0.3, mgp=c(3,.3,0))


plot(0:24, 0:24, type="n", xlab="True node age", ylab="Estimated node age"); abline(a=0, b=1)

points(BDages, CHages, bg=gray(1), pch=21) 

for(i in 1:30) { points(BDages[i,], CHages[i,], bg=i, pch=21)  }

for(i in 1:30) { text(BDages[i,], CHages[i,], lab=i, cex=0.5, col=i)  }
