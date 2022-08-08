##########################################################
### Simulations of Trees, Fossils, etc ###
##########################################################
 
# Use unconstrained BD tree instead of age = 20
# Use a sampled/incomplete tree (frac=0.2) to distribute branching times more evenly


library(ape)
library(phangorn)
library(TreeSim)
library(FossilSim)
#library(devtools)
#install_github("evolucionario/CladeDate")
library(CladeDate)

setwd("~/Documents/Avian Time Trees/Simulated Datasets BD2")


# Ancillary function #

### Function that collects the fossil record (ages) of a particular clade indicated by the MRCA node number (x)
# x the node number
# phy a phylogenetic tree of class 'phylo'
# fossils an object of class 'fossils'
# stem	if TRUE, fossils along the basal stem lienage are also reported

fossil.record <- function(x, phy, fossils, stem=FALSE) {

# Phangorn function Descendants identifies the descendant nodes and tips
	
	des <- Descendants(phy, node=x, type = "all")

# Then use the descendants indices to retrieve the fossil ages from the fossil object
	fages <- numeric()

	for(i in 1:length(des)) {
	fages <- c(fages, fossils$hmin[which(fossils$edge == des[i])])
	}
	
	# Get tip names to faciliate the generation of the clade constraint in MrBayes

	tnames <- phy$tip.label[Descendants(phy, node=x, type = "tips")[[1]]]

	if(stem==TRUE) {
		
		stemfage <- fossils$hmin[which(fossils$edge == x)]
		
		fages <- c(stemfage, fages) 
	}

	fages <- sort(fages, decreasing=TRUE)
	
	RES <- list(mrca.node = x, tip.names=tnames, fossil.ages=fages, n.foss=length(fages))
	
	return(RES)
}





#############################
### Simulation Parameters ###
#############################

RATE <- 0.02


#################
### MAIN LOOP ###
#################

# Create a big loop that repeats the script 

REPS <- 1:100

for(i in REPS) {


#####################
### Simulate tree ###
#####################

# using TreeSim function in a repeat loop to make sure basal sister taxa both have more than 2 species (so at least one internal node)
# balance calculates the numer of descendants for each dougther clade and the first entry is the root node

repeat {
	
tr <- sim.bd.taxa(n=10, numbsim=1, lambda=0.1, mu=0, frac=0.2, complete=FALSE)[[1]]

plot(tr)
#nodelabels()

if(!any(balance(tr)[1,] < 2)) break }

# delete root edge (creates problems in FossiSim)

tr$root.edge <- NULL

write.tree(tr, file=paste0("Simulated",i,".tre"))


### Add outgroup to tree ###

tr1 <- tr

# Create a 2 Ma root edge to attach the outgroup
tr1$root.edge <- 2 

AGE <- max(branching.times(tr1))

# Create the outgroup linage: a tree with a single specie and a 22 Ma branch
og.edge <- read.tree(text = paste0("(og:", AGE+2, ");"))

tr2 <- bind.tree(tr1, og.edge, position=2) 

# force to be ultrametric to correct rounding errors

tr2 <- force.ultrametric(tr2, method=c("nnls"))

plot(tr2); axisPhylo()

#is.ultrametric(tr2); is.binary(tr2)


####################################
### Simulate molecular sequences ###
####################################

# Basic method with high stochasticity due to short sequence
# The idea is that the sequences, which are simulated in a strict clock fashion, do not dominate the results

# JC model:

DNA <- simSeq(tr2, l = 1000, type = "DNA", rate = RATE)

write.FASTA(as.DNAbin(DNA), file=paste0("SimulatedDNA.",i,".fas"))


## Alternative ##
# Add a jiggle to branch lengths to simulate rate heterogeneity
#tr2 <- tr
#tr2$edge.length <- tr$edge.length * rlnorm(length(tr$edge.length), 0.1, 0.3)
#tr2$edge.length <- tr$edge.length * rnorm(length(tr$edge.length), 1, 0.3)
#plot(tr2)
# JC model with long sequences so the data carry the non-clock signal
#DNA <- simSeq(tr2, l = 4000, type = "DNA", rate = 0.02)
# But makes estimation more dificult, involvining the clockrate


######################################################################
### Infer a ML tree with DNA sequences and fixed original topology ###
######################################################################

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




####################################
### Simulate fossils on branches ###
####################################

### Different fossilization in different clades ###


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

plot(Fos2, tr); nodelabels(); tiplabels(tr$tip.label)

# Obtain fossil recod

fr.clade1 <- fossil.record(calib.nodes[1], tr, Fos2)

fr.clade2 <- fossil.record(calib.nodes[2], tr, Fos2)
 
if(all(c(fr.clade1$n.fos, fr.clade2$n.fos) > 2)) break

}

save(fr.clade1, file=paste0("FossilRecordClade1.",i,".R"))
save(fr.clade2, file=paste0("FossilRecordClade2.",i,".R"))

save(Fos2, file=paste0("FossilSim",i,".R"))


#############################################
### Estimate Clade Age from Fossil Record ###
#############################################

# Using the StraussSadler method and the gamma function to simplify

date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, PDFfitting="lognormal", KStest=TRUE)

if(date.clade1$KStest$p.value < 0.05) {

date.clade1 <- clade.date(ages=fr.clade1$fossil.ages, method="Solow", PDFfitting="lognormal")

	}


date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, PDFfitting="lognormal", KStest=TRUE)

if(date.clade2$KStest$p.value < 0.05) {

date.clade2 <- clade.date(ages=fr.clade2$fossil.ages, method="Solow", PDFfitting="lognormal")

	}

save(date.clade1, file=paste0("Clade1date.",i,".R"))

save(date.clade2, file=paste0("Clade2date.",i,".R"))

cat(paste("\nReplicate",i,"completed\n"))

}

##### END OF MAIN LOOP ######
