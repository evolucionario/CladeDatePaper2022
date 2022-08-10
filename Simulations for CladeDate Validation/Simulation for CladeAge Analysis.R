##############################################
### Simulation for the CladeAge Comparison ###
##############################################

### Simulation that creates a dataset (tree, fossils, DNA alignment) to be used with CladeAge and BEST2 so the results can be compared to the new CladeDate + Chronos.


library(TreeSim)


#####################
### Simulate tree ###

SEED <- 99

# using TreeSim function in a 'repeat' loop to make sure basal sister taxa both have more than 2 species (so at least one internal node)
# balance calculates the numer of descendants for each dougther clade and the first entry is the root node

#tr <- sim.bd.taxa.age(age=AGE, n=10, numbsim=1, lambda=0.1, mu=0, frac = 1, mrca = TRUE)[[1]]

tr <- read.tree(text="((((t7:11.06730802,(t3:5.30770791,t9:5.30770791):5.759600115):4.673794634,t8:15.74110266):0.5778341325,t6:16.31893679):3.681063209,((t2:4.813437322,t1:4.813437322):5.261159269,((t4:2.856380757,t5:2.856380757):5.879688258,t10:8.736069015):1.338527576):9.925403409):0;")

plot(tr); nodelabels()


# delete root edge (creates problems in FossiSim)

tr$root.edge <- NULL

write.tree(tr, file="Simulated0.tre")
tr <- read.tree(file="Simulated0.tre")

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

# JC model:

set.seed(SEED)
DNA <- simSeq(tr2, l = 1000, type = "DNA", rate = RATE)

write.FASTA(as.DNAbin(DNA), file=paste0("SimulatedDNA0.fas"))

### END ###


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

write.tree(MLtree3, file=paste0("EstimatedML0.tre"))

MLtree3 <- read.tree(file="EstimatedML0.tre")

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
 


save(fr.clade1, file=paste0("FossilRecordClade1.0.R"))
save(fr.clade2, file=paste0("FossilRecordClade2.0.R"))

save(Fos2, file=paste0("FossilSim0.R"))

#load(file="FossilRecordClade1.0.R")
#load(file="FossilRecordClade2.0.R")

### END ###


#############################################
### Estimate Clade Age from Fossil Record ###

# Using the StraussSadler method and the gamma function to simplify

date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, PDFfitting="lognormal", KStest=TRUE)

if(date.clade1$KStest$p.value < 0.05) {

date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, method="RobsonWhitlock", PDFfitting="lognormal")

	}


save(date.clade1, file=paste0("Clade1date0.R"))

#$Quantiles
#      0%      50%      95% 
#15.21580 16.39652 21.34264 

#$KStest

#	Exact one-sample Kolmogorov-Smirnov test

#data:  Mages
#D = 0.21748, p-value = 0.7712
#alternative hypothesis: two-sided


#$PDFfit.model
#[1] "lognormal"

#$PDFfit
#     meanlog        sdlog   
#  0.003230185   1.331636473 
# (0.013316365) (0.009416092)
 



date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, PDFfitting="lognormal", KStest=TRUE)

if(date.clade2$KStest$p.value < 0.05) {

date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, method="RobsonWhitlock", PDFfitting="lognormal")

	} 
 
save(date.clade2, file=paste0("Clade2date0.R"))

#$Quantiles
#       0%       50%       95% 
# 8.603469  9.008276 10.509394 

#$KStest

#	Exact one-sample Kolmogorov-Smirnov test

#data:  Mages
#D = 0.19014, p-value = 0.6255
#alternative hypothesis: two-sided


#$PDFfit.model
#[1] "lognormal"

#$PDFfit
#     meanlog         sdlog    
#  -1.089524008    1.307049413 
# ( 0.013070494) ( 0.009242235)



### Get oldest fossil of each clade to use in CladeAge ###

# Use fossil.record() to obtain the ages of the fossils on each branch.
# Use this information together with the FASTA file of DNA alignments and the tree to build the xml file for BEAST2.
# See: https://taming-the-beast.org/tutorials/CladeAge-Tutorial/



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

Chrono <- ladderize(Chrono)

tr <- ladderize(tr)

plot(tr$edge.length, Chrono$edge.length, las=1); abline(0,1)


truebt <- branching.times(tr)

Chronobt <- branching.times(Chrono)

plot(truebt, Chronobt, las=1); abline(0,1)



### Source compareBTs function to compare correct branching times
# (ladderize and order are not putting the trees in the same order for some reason)

source('Compare Branching Times.R')


### Read CladeAge-Beast MCC tree ###

# Load the result from CladeAge + BEAST2 analysis (the Maximum Clade Credibility tree)

Btree <- read.nexus("MCCtree.tre")

Btree$root.edge <- NULL

plot(Btree); axisPhylo()

Btree <- reorder(Btree)

Btree <- ladderize(Btree, right=FALSE)

Btreebt <- branching.times(Btree)


quartz("trees",4, 9); par(mfcol=c(3,1), las=1)

plot(tr); axisPhylo()
plot(Chrono); axisPhylo()
plot(Btree); axisPhylo()



###################
### Final Plot ####
###################


quartz("trees",5,4.5);par(las=1)

plot(compareBTs(tr, Chrono)$BT, type="n", xlab="True branching times", ylab="Estimated branching times"); abline(0,1, lty="dashed")

points(compareBTs(tr, Chrono)$BT, col=rgb( .8, .4, 0), cex=1.8, pch=16)

points(compareBTs(tr, Btree)$BT, col=rgb( 0, .45, .7), cex=1.8, pch=18)

legend("bottomright", legend=c("CladeDate + chronos", "CladeAge + Beast2"), cex=0.9, bty="n", pch=c(16, 18), pt.cex=1.3, col=c(rgb( .8, .4, 0), rgb( 0, .45, .7)))
