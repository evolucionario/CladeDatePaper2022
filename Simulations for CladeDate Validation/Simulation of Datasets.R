##############################################################
### SIMULATION OF TREES, DNA SEQUENCES, AND FOSSIL RECORDS ###
##############################################################

# by Santiago Claramunt - claramunt.bio@gmail.com


### LOAD REQUIRED PACKAGES ###

library(ape)
library(phangorn)
library(phytools)
library(TreeSim)
library(FossilSim)
#library(devtools)
#install_github("evolucionario/CladeDate")
library(CladeDate)


### Set the Working Directory ###

#setwd("DIR")

#################
### MAIN LOOP ###


### SIMULATION PARAMETERS ###

RATE <- 0.02

REPS <- 1:100

# Loop that repeats the script 

for(i in REPS) {
	
### SIMULATE TREE ###

# using TreeSim function in a 'repeat' loop to make sure basal sister taxa both have more than 2 species (so at least one internal node to calibrate)
# balance calculates the numer of descendants for each dougther clade and the first entry is the root node

repeat {
	
tr <- sim.bd.taxa(n=10, numbsim=1, lambda=0.1, mu=0, frac=0.2, complete=FALSE)[[1]]

plot(tr)
#nodelabels()

if(!any(balance(tr)[1,] < 2)) break }

# delete root edge (creates problems in FossiSim)

tr$root.edge <- NULL

write.tree(tr, file=paste0("Simulated",i,".tre"))


# Add an outgroup to tree (needed for rooting the tree after phylogenetic analysis) #

tr1 <- tr

# Create a 2 Ma root edge to attach the outgroup

tr1$root.edge <- 2 

AGE <- max(branching.times(tr1))

# Create the outgroup linage: a tree with a single specie and a 22 Ma branch
og.edge <- read.tree(text = paste0("(og:", AGE+2, ");"))

tr2 <- bind.tree(tr1, og.edge, position=2) 

# force the resultant tree to be ultrametric to correct for rounding errors

tr2 <- force.ultrametric(tr2, method=c("nnls"))


### SIMULATE DNA SEQUENCES ###

# Basic method with high stochasticity due to short sequence length
# The idea is that the sequences, which are simulated assuming a strict molecular clock, do not dominate the results

# JC model:

DNA <- simSeq(tr2, l = 1000, type = "DNA", rate = RATE)

write.FASTA(as.DNAbin(DNA), file=paste0("SimulatedDNA.",i,".fas"))


### INFER MAXIMUM LIKELIHOOD PHYLOGENY ###

# Infer a ML tree using simulated DNA sequences and the original topology fixed
	
MLfit <- pml(tr2, DNA, rate = RATE)

MLtree <- optim.pml(MLfit)

# the resultant tree is unrooted
#MLtree <- optim.pml(MLfit, optRooted=TRUE) # But produced an ultrametric tree

plot(MLtree)

# Reroot the tree #

MLtree2 <- root(MLtree$tr, outgroup="og")

# then delete the outgroup
MLtree3 <- drop.tip(MLtree2, tip="og")

plot(MLtree3)

write.tree(MLtree3, file=paste0("EstimatedML",i,".tre"))

### ML tree ready ###


### SIMULATE FOSSILS ON BRANCHES ###

## Different fossilization in different clades ##

# Find the two descendant nodes of the root node

calib.nodes <- Children(tr, node=11)

# Count the number of branches in the second subclade, which equals the sum of tips and nodes descendant from the node

branches.clade2 <- length(Descendants(tr, node=calib.nodes[2], type="all"))


# Create vector of rates (one per branch)
# first value for the root edge that is zero but exists

Rates <- c(rep(0.1, Nedge(tr)-branches.clade2), rep(0.5, branches.clade2))

# Simulate fossil finds with heterogeneous rates
# conditional repeat ensures there is more than 1 fossil in each clade, otherwise the simulation is repeated

tr$root.edge <- NULL

repeat {

Fos2 <- sim.fossils.poisson(rate=Rates, tree=tr)

# Obtain fossil recod

fr.clade1 <- fossil.record(calib.nodes[1], tr, Fos2)

fr.clade2 <- fossil.record(calib.nodes[2], tr, Fos2)
 
if(all(c(fr.clade1$n.fos, fr.clade2$n.fos) > 2)) break

}

save(fr.clade1, file=paste0("FossilRecordClade1.",i,".R"))
save(fr.clade2, file=paste0("FossilRecordClade2.",i,".R"))

save(Fos2, file=paste0("FossilSim",i,".R"))


cat(paste("\nReplicate",i,"completed\n"))

}


### END OF MAIN LOOP ###
########################

	
###########
### END ###
###########
