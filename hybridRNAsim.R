# simulate ortholog-homeolog states

# notes
# o.de is the actual number of DE genes, while o.up + o.down are the observed number of DE genes
# these values are not comparable if o.up + o.down must be signficant *and* 2-fold different

# load libraries
library(compcodeR)

# load functions
source("hybridRNAsim_functions.R")

# variables
iterations				<- 5
n.genes					<- 1000
scalar					<- 10
bio.reps				<- 2
		
significance			<- 0.05
fold.change				<- 2

o.n.diffexp				<- 50
o.seqdepth				<- 2.5e7
o.fraction.upregulated 	<- 0.5

h1.n.diffexp			<- 50
h1.seqdepth				<- 2.5e7
h1.fraction.upregulated	<- 0.5

second.hybrid			<- FALSE
h2.n.diffexp			<- 50
h2.seqdepth				<- 2.5e7
h2.fraction.upregulated	<- 0.5


# set up dataframe
if(second.hybrid == FALSE){
	expression.res <- as.data.frame(matrix(0, ncol=26, nrow=0))
	colnames(expression.res) <- c(
	"ngenes", "nreps", "o.nreads", "o.de", "o.frac.up", "h1.nreads", "h1.de", "h1.frac.up",
	"o.up", "o.down", "h1.up", "h1.down", "o.up.fc", "o.down.fc", "h1.up.fc", "h1.down.fc",
	"I1", "L1", "S1", "R1", "U1", "I1.fc", "L1.fc", "S1.fc", "R1.fc", "U1.fc")
}else if(second.hybrid == TRUE){
	expression.res <- as.data.frame(matrix(0, ncol=43, nrow=0))
	colnames(expression.res) <- c(
	"ngenes", "nreps", "o.nreads", "o.de", "o.frac.up", "h1.nreads", "h1.de", "h1.frac.up", 
	"h2.nreads", "h2.de", "h2.frac.up", "o.up", "o.down", "h1.up", "h1.down", "h2.up", 
	"h2.down", "o.up.fc", "o.down.fc", "h1.up.fc", "h1.down.fc", "h2.up.fc", "h2.down.fc",
	"I1", "L1", "S1", "R1", "U1", "I1.fc", "L1.fc", "S1.fc", "R1.fc", "U1.fc",
	"I2", "L2", "S2", "R2", "U2", "I2.fc", "L2.fc", "S2.fc", "R2.fc", "U2.fc")
}

# example run (one hybrid only)
x <- runsims(iterations, n.genes, scalar, bio.reps, significance, fold.change, o.n.diffexp, o.seqdepth, 
	o.fraction.upregulated, h1.n.diffexp, h1.seqdepth, h1.fraction.upregulated)

expression.res <- rbind(expression.res, x)
expression.res

# appended example run (one hybrid only)
x <- runsims(iterations, n.genes, scalar, bio.reps, significance, fold.change, o.n.diffexp, o.seqdepth, 
	o.fraction.upregulated, h1.n.diffexp, h1.seqdepth, h1.fraction.upregulated)

expression.res <- rbind(expression.res, x)
expression.res
