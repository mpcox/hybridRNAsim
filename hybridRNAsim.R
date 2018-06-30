# simulate ortholog-homeolog states

# notes
# o.de is the actual number of DE genes, while o.up + o.down are the observed number of DE genes
# these values are not comparable if o.up + o.down must be signficant *and* 2-fold different

# load libraries
library(compcodeR)

# load functions
source("hybridRNAsim_functions.R")

# variables
iterations <- 5
n.genes <- 1000
scalar <- 10
bio.reps <- 2

significance <- 0.05
fold.change <- 2

o.n.diffexp <- 50
o.seqdepth <- 2.5e7
o.fraction.upregulated <- 0.5
h.n.diffexp <- 50
h.seqdepth <- 2.5e7
h.fraction.upregulated <- 0.5

# setup dataframe
expression.res <- as.data.frame(matrix(0, ncol=26, nrow=0))
colnames(expression.res) <- c(
	"ngenes", "nreps", "o.nreads", "o.de", "o.frac.up", "h.nreads", "h.de", "h.frac.up",
	"o.up", "o.down", "h.up", "h.down", "o.up.fc", "o.down.fc", "h.up.fc", "h.down.fc",
	"I", "L", "S", "R", "U", "I.fc", "L.fc", "S.fc", "R.fc", "U.fc")

# example run
x <- runsims(iterations, n.genes, scalar, bio.reps, significance, fold.change, o.n.diffexp, o.seqdepth, 
	o.fraction.upregulated, h.n.diffexp, h.seqdepth, h.fraction.upregulated)
	
expression.res <- rbind(expression.res, x)
expression.res

# appended example run
x <- runsims(iterations, n.genes, scalar, bio.reps, significance, fold.change, o.n.diffexp, o.seqdepth, 
	o.fraction.upregulated, h.n.diffexp, h.seqdepth, h.fraction.upregulated)
	
expression.res <- rbind(expression.res, x)
expression.res
	