#################################################
### Comparison Between CladeDate and CladeAge ###
#################################################

### Simulation that creates a dataset (tree, fossils, DNA alignment) to be used with CladeAge and BEST2 so the results can be compared to the new CladeDate + Chronos.


### Load Required Packages ###

library(ape)
library(phangorn)
library(CladeDate)
library(TreeSim)
library(FossilSim)


setwd('~/Documents/Avian Time Trees/CladeDatePaper2022Simulations/Simulations for CladeDate Validation')

# Set a seed for the random number

SEED <- 1


#####################
### Simulate tree ###

# Simulate a birth-death tree. The function doesn't return the same tree with the provided seed so the original tree si provides so the same xml files for BEAST can be used.
#set.seed(SEED)
#tr <- sim.bd.taxa.age(age=20, n=10, numbsim=1, lambda=0.1, mu=0, frac = 1, mrca = TRUE)[[1]]

tr <- read.tree(text="((((t4:1.925737249,t9:1.925737249):0.5303991474,t8:2.456136396):6.566717035,(t6:2.306491955,t1:2.306491955):6.716361476):10.97714657,(((t10:4.999615538,t7:4.999615538):6.866952518,(t5:5.918453114,t3:5.918453114):5.948114941):4.260349859,t2:16.12691791):3.873082085);")

plot(tr); nodelabels()


# delete root edge (creates problems in FossiSim)

tr$root.edge <- NULL

write.tree(tr)

write.tree(tr, file="Simulated0.tre")

### Add outgroup to tree ###

tr1 <- tr

# Create a 2 Ma root edge to attach the outgroup
tr1$root.edge <- 2 

# Create the outgroup linage: a tree with a single specie and a 22 Ma branch
og.edge <- read.tree(text = "(og:22);")

tr2 <- bind.tree(tr1, og.edge, position=2) 

#reroot (cause bind.tree unroots the tree)
tr2 <- root(tr2, outgroup="og", resolve.root = TRUE)

plot(tr2); axisPhylo()

#is.ultrametric(tr2); is.binary(tr2)

### END OF TREE SIMULATION ###
##############################




###################################
### GENERATION OF DNA SEQUENCES ###

# Basic method with high stochasticity due to short sequence length
# The idea is that the sequences, which are simulated in a strict clock fashion, do not dominate the results

# JC model with 

RATE <- 0.02

set.seed(SEED)
DNA <- simSeq(tr2, l = 1000, type = "DNA", rate = RATE)

write.FASTA(as.DNAbin(DNA), file="SimulatedDNA0.fas")

### DNA sequences ready ###


######################################################################
### Infer a ML tree with DNA sequences and fixed original topology ###

MLfit <- pml(tr2, DNA, rate = RATE)

MLtree <- optim.pml(MLfit)

# the resultant tree is unrooted
#MLtree <- optim.pml(MLfit, optRooted=TRUE) # But produced an ultrametric tree

plot(MLtree)

# Root the tree #

MLtree2 <- root(MLtree$tr, outgroup="og")

# then delete the outgroup
MLtree3 <- drop.tip(MLtree2, tip="og")

plot(MLtree3)

write.tree(MLtree3, file="EstimatedML0.tre")

### ML tree ready ###




####################################
### Simulate Fossils on Branches ###

### Different fossilization rates in different clades ###

# Find the two descendant nodes of the root node

calib.nodes <- Children(tr, node=11)

# Count the number of branches in the second subclade, which equals the sum of tips and nodes descendant from the node

branches.clade2 <- length(Descendants(tr, node=calib.nodes[2], type="all"))


# Create vector of rates (one per branch)
# first value for the root edge that is of length zero but exists

Rates <- c(rep(0.1, Nedge(tr)-branches.clade2), rep(0.5, branches.clade2))

# Simulate fossil finds with heterogeneous rates
# conditional repeat ensures there is more than 1 fossil in each clade, otherwise the simulation is repeated

tr$root.edge <- NULL

set.seed(SEED)
Fos2 <- sim.fossils.poisson(rate=Rates, tree=tr)

plot(Fos2, tr); nodelabels(); tiplabels(tr$tip.label)

# Make sure there is more than one fossil in each basal clade

# Obtain fossil recod

fr.clade1 <- fossil.record(calib.nodes[1], tr, Fos2)

fr.clade2 <- fossil.record(calib.nodes[2], tr, Fos2)
 

### END ###




###################################
### Prepare xml file for BEAST2 ###

# Delete the outgroup and load the FASTA file of DNA alignments into BEAUti and setup a basic CladeAge analysis.

# Specify a JC69 model and a strict clock model.

# Specify a starting tree using the simulated true tree: write.tree(tr)
# Otherwise, Beast2 may have a hard time finding a tree compatible with all clade and time constraints specified below.

# Identify the oldest fossil on each branch to set up clades and constraitns.
# See: https://taming-the-beast.org/tutorials/CladeAge-Tutorial/

cbind(Fos2[,"edge"], Fos2[,"hmin"])

plot(Fos2, tr); nodelabels(pos=2, frame="n", font=2); tiplabels(frame="n", pos=2)

tiplabels(tr$tip.label, pos=4, frame="n", font=3, cex=1.3)

# Note that CladeAge uses the stem age of clades so use the oldest fossil in the stem of the clade if present.

# Set boundaries for the fossil sampling rate parameters:
 
#                    <parameter id="minSamplingRate" spec="parameter.RealParameter" name="minSamplingRate">0.05</parameter>
#                    <parameter id="maxSamplingRate" spec="parameter.RealParameter" name="maxSamplingRate">5.00</parameter>

# Set the net diversificationr rate parameter at the empirical rate: 10 species in 20 million years log(10/2)/20 = 0.8 

# Disable the operators that change the topology so only branch length are sampled using the same original topology, if provided as initial tree.

# A CladeAge.xml is provided as a reference but clade dates may need to be modified if new simulated fossil are used.

# Run the xml file in BEAST2






################################################################
### Estimate the Age of the Two Basal Clades Using CladeDate ###

date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, KStest=TRUE, repvalues=FALSE)

summary.clade.date(date.clade1)

	Exact one-sample Kolmogorov-Smirnov test

data:  Mages
D = 0.46449, p-value = 0.4207
alternative hypothesis: two-sided


Quantiles:
    0%    50%    95% 
 6.846  8.501 18.230 

Parameters of the gamma function:
offset  shape   rate 
6.8460 0.7401 0.2295 
 



date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, KStest=TRUE, repvalues=FALSE)
 
summary.clade.date(date.clade2)

	Exact one-sample Kolmogorov-Smirnov test

data:  Mages
D = 0.17469, p-value = 0.2678
alternative hypothesis: two-sided


Quantiles:
   0%   50%   95% 
15.20 15.53 16.65 

Parameters of the gamma function:
 offset   shape    rate 
15.2023  0.9791  2.0355 




#############################
### Estimate with Chronos ###
#############################

# Point calibration using the median 

Calib <- makeChronosCalib( phy = MLtree3,
node = c(fr.clade1$mrca.node, fr.clade2$mrca.node),
age.min = c(date.clade1$Quantiles[2], date.clade2$Quantiles[2]) )

Chrono <- chronos( phy = MLtree3, model="clock", calibration=Calib)

cat(attr(Chrono, "message"))

plot(Chrono); axisPhylo()



### Source compareBTs function to compare correct branching times
# (Even without changing the topology, the trees obtained from BEAST may not be directy comparable.)

source('Compare Branching Times.R')




#############
### Plot ####
#############


### Read CladeAge-BEAST Posterior ###

# Load results from CladeAge + BEAST2 analysis

ptrees <- read.nexus("BeastPosterior.trees") 


### Color Definitions ####

col.ca <- rgb( 0, .45, .7)
col.ca.bg <- rgb( 0, .6, .5, .1)
col.cd <- rgb( .8, .3, 0)
col.cd.bg <- rgb( .8, .4, 0)


### Plot ####

quartz("trees",4.5,4.5); par(mgp=c(1.8,0.4,0), tcl=-0.3, las=1, cex.axis=0.8, xaxs="i", yaxs="i")

plot(compareBTs(tr, Chrono)$BT, type="n", xlab="True branching times", ylab="Estimated branching times", xlim=c(0,23), ylim=c(0,23))

abline(0,1)

for(i in seq(from=100, to=1000, by=1)) {
points(compareBTs(tr, ptrees[[i]])$BT, col= col.ca.bg, bg= col.ca.bg, cex=1.2, pch=21)
}

points(compareBTs(tr, Chrono)$BT, col=col.cd, bg=col.cd.bg, cex=1.2, pch=23)

legend("bottomright", legend=c("CladeDate + chronos", "CladeAge + BEAST2"), cex=0.8, bty="n", pch=c(23,21), pt.cex=1.4, col=c(col.cd, col.ca), pt.bg=c(col.cd.bg, rgb( 0, .6, .5, .7)))

### END ###
###########
