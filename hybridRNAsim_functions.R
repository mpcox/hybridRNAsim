# functions
# simulate ortholog-homeolog states

# classify genes to ortholog-homeolog classes
classify <- function(vector){
	
	ortholog <- vector[1]
	homeolog <- vector[2]
	
	if     ( ortholog ==  0 && homeolog ==  0 ){ class = 'I' }
	else if( ortholog == +1 && homeolog == +1 ){ class = 'I' }
	else if( ortholog == -1 && homeolog == -1 ){ class = 'I' }
	else if( ortholog == +1 && homeolog ==  0 ){ class = 'L' }
	else if( ortholog == -1 && homeolog ==  0 ){ class = 'L' }
	else if( ortholog ==  0 && homeolog == +1 ){ class = 'S' }
	else if( ortholog ==  0 && homeolog == -1 ){ class = 'S' }
	else if( ortholog == +1 && homeolog == -1 ){ class = 'R' }
	else if( ortholog == -1 && homeolog == +1 ){ class = 'R' }
	else { class = 'U' }
	
	return(class)
}

# count numbers of genes in ortholog-homeolog classes
count <- function(vector){
	
	class.I <- length( which(vector == 'I') )
	class.L <- length( which(vector == 'L') )
	class.S <- length( which(vector == 'S') )
	class.R <- length( which(vector == 'R') )
	class.U <- length( which(vector == 'U') )
	
	return( c(class.I, class.L, class.S, class.R, class.U) )
}

# wrapper script for running simulations and analyses
runsims <- function(iterations, n.genes, scalar, bio.reps, significance, fold.change, o.n.diffexp, o.seqdepth, o.fraction.upregulated, h1.n.diffexp, h1.seqdepth, h1.fraction.upregulated, second.hybrid=FALSE, h2.n.diffexp=NULL, h2.seqdepth=NULL, h2.fraction.upregulated=NULL){
	
	# check for second hybrid
	if(second.hybrid == TRUE){
		if( is.null(h2.n.diffexp) || is.null(h2.seqdepth) || is.null(h2.fraction.upregulated) ){
			stop("Must define values for h2.n.diffexp, h2.seqdepth and h2.fraction.upregulated")
		}
	}
	
	# make data frame
	if(second.hybrid == FALSE){
		this.res <- as.data.frame(matrix(0, ncol=26, nrow=0))
		colnames(this.res) <- c(
		"ngenes", "nreps", "o.nreads", "o.de", "o.frac.up", "h1.nreads", "h1.de", "h1.frac.up",
		"o.up", "o.down", "h1.up", "h1.down", "o.up.fc", "o.down.fc", "h1.up.fc", "h1.down.fc",
		"I1", "L1", "S1", "R1", "U1", "I1.fc", "L1.fc", "S1.fc", "R1.fc", "U1.fc")
	}else if(second.hybrid == TRUE){
		this.res <- as.data.frame(matrix(0, ncol=43, nrow=0))
		colnames(this.res) <- c(
		"ngenes", "nreps", "o.nreads", "o.de", "o.frac.up", "h1.nreads", "h1.de", "h1.frac.up", 
		"h2.nreads", "h2.de", "h2.frac.up", "o.up", "o.down", "h1.up", "h1.down", "h2.up", 
		"h2.down", "o.up.fc", "o.down.fc", "h1.up.fc", "h1.down.fc", "h2.up.fc", "h2.down.fc",
		"I1", "L1", "S1", "R1", "U1", "I1.fc", "L1.fc", "S1.fc", "R1.fc", "U1.fc",
		"I2", "L2", "S2", "R2", "U2", "I2.fc", "L2.fc", "S2.fc", "R2.fc", "U2.fc")
	}
	
	# record simulation parameters
	if(second.hybrid == FALSE){
		conditions <- c(n.genes, bio.reps, o.seqdepth, o.n.diffexp, o.fraction.upregulated, h1.seqdepth, h1.n.diffexp, h1.fraction.upregulated)
	}else if(second.hybrid == TRUE){
		conditions <- c(n.genes, bio.reps, o.seqdepth, o.n.diffexp, o.fraction.upregulated, h1.seqdepth, h1.n.diffexp, h1.fraction.upregulated, h2.seqdepth, h2.n.diffexp, h2.fraction.upregulated)
	}

	# simulate orthologs
	o.sim <- generateSyntheticData(	dataset = "o.sim", 
									n.vars = n.genes * scalar,
									samples.per.cond = bio.reps,
									n.diffexp = o.n.diffexp * scalar,
									seqdepth = o.seqdepth * scalar,
									fraction.upregulated = o.fraction.upregulated,
									output.file = "o.sim.rds")
	
	runDiffExp(data.file = "o.sim.rds", result.extent = "DESeq2",
	               Rmdfunction = "DESeq2.createRmd",
	               output.directory = ".", fit.type = "parametric",
	               test = "Wald")
	
	o.sim.dat <- readRDS("o.sim_DESeq2.rds")
	
	# simulate homeologs
	h1.sim <- generateSyntheticData(dataset = "h1.sim", 
									n.vars = n.genes * scalar,
									samples.per.cond = bio.reps,
									n.diffexp = h1.n.diffexp * scalar,
									seqdepth = h1.seqdepth * scalar,
									fraction.upregulated = h1.fraction.upregulated,
									output.file = "h1.sim.rds")
	
	runDiffExp(data.file = "h1.sim.rds", result.extent = "DESeq2",
	               Rmdfunction = "DESeq2.createRmd",
	               output.directory = ".", fit.type = "parametric",
	               test = "Wald")
	
	h1.sim.dat <- readRDS("h1.sim_DESeq2.rds")
	
	# simulate homeologs for second hybrid
	if(second.hybrid == TRUE){
		h2.sim <- generateSyntheticData(dataset = "h2.sim", 
										n.vars = n.genes * scalar,
										samples.per.cond = bio.reps,
										n.diffexp = h2.n.diffexp * scalar,
										seqdepth = h2.seqdepth * scalar,
										fraction.upregulated = h2.fraction.upregulated,
										output.file = "h2.sim.rds")
		
		runDiffExp(data.file = "h2.sim.rds", result.extent = "DESeq2",
		               Rmdfunction = "DESeq2.createRmd",
		               output.directory = ".", fit.type = "parametric",
		               test = "Wald")
		
		h2.sim.dat <- readRDS("h2.sim_DESeq2.rds")		
	}
	
	for( i in 1:iterations ){
		
		# make subsampled dataset (note: datasets can differ in size by 1)
		o.values  <- sample(1:length(o.sim.dat@result.table[,1]),  n.genes, replace=F)
		h1.values <- sample(1:length(h1.sim.dat@result.table[,1]), n.genes, replace=F)
		
		this.o.sim.dat  <- o.sim.dat@result.table[o.values,]
		this.h1.sim.dat <- h1.sim.dat@result.table[h1.values,]
		
		if(second.hybrid == TRUE){
			h2.values <- sample(1:length(h2.sim.dat@result.table[,1]), n.genes, replace=F)
			this.h2.sim.dat <- h2.sim.dat@result.table[h2.values,]
		}
		
		# determine differential expression (statistically significant)
		o.sign.up   <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC >= 0 )
		o.sign.down <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC < 0 )
		o.nonsign   <- setdiff(1:n.genes, c(o.sign.up, o.sign.down))
		
		h1.sign.up   <- which( this.h1.sim.dat$adjpvalue <= significance & this.h1.sim.dat$logFC >= 0 )
		h1.sign.down <- which( this.h1.sim.dat$adjpvalue <= significance & this.h1.sim.dat$logFC < 0 )
		h1.nonsign   <- setdiff(1:n.genes, c(h1.sign.up, h1.sign.down))
		
		if(second.hybrid == TRUE){
			h2.sign.up   <- which( this.h2.sim.dat$adjpvalue <= significance & this.h2.sim.dat$logFC >= 0 )
			h2.sign.down <- which( this.h2.sim.dat$adjpvalue <= significance & this.h2.sim.dat$logFC < 0 )
			h2.nonsign   <- setdiff(1:n.genes, c(h2.sign.up, h2.sign.down))
		}
		
		if(second.hybrid == FALSE){
			sim.nde <- c(length(o.sign.up), length(o.sign.down), length(h1.sign.up), length(h1.sign.down))
		}else if(second.hybrid == TRUE){
			sim.nde <- c(length(o.sign.up), length(o.sign.down), length(h1.sign.up), length(h1.sign.down), length(h2.sign.up), length(h2.sign.down))
		}
	
		# determine differential expression (statistically significant and >=x-fold change)	
		o.sign.up.fc   <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC >= log2(fold.change) )
		o.sign.down.fc <- which( this.o.sim.dat$adjpvalue <= significance & this.o.sim.dat$logFC <= -log2(fold.change) )
		o.nonsign.fc   <- setdiff(1:n.genes, c(o.sign.up.fc, o.sign.down.fc))
		
		h1.sign.up.fc   <- which( this.h1.sim.dat$adjpvalue <= significance & this.h1.sim.dat$logFC >= log2(fold.change) )
		h1.sign.down.fc <- which( this.h1.sim.dat$adjpvalue <= significance & this.h1.sim.dat$logFC <= -log2(fold.change) )
		h1.nonsign.fc   <- setdiff(1:n.genes, c(h1.sign.up.fc, h1.sign.down.fc))

		if(second.hybrid == TRUE){
			h2.sign.up.fc   <- which( this.h2.sim.dat$adjpvalue <= significance & this.h2.sim.dat$logFC >= log2(fold.change) )
			h2.sign.down.fc <- which( this.h2.sim.dat$adjpvalue <= significance & this.h2.sim.dat$logFC <= -log2(fold.change) )
			h2.nonsign.fc   <- setdiff(1:n.genes, c(h2.sign.up.fc, h2.sign.down.fc))
		}
		
		if(second.hybrid == FALSE){
			sim.nde.fc <- c(length(o.sign.up.fc), length(o.sign.down.fc), length(h1.sign.up.fc), length(h1.sign.down.fc))
		}else if(second.hybrid == TRUE){
			sim.nde.fc <- c(length(o.sign.up.fc), length(o.sign.down.fc), length(h1.sign.up.fc), length(h1.sign.down.fc), length(h2.sign.up.fc), length(h2.sign.down.fc))
		}
		
		# build change vectors (statistically significant)
		orthologs <- vector(length=n.genes)
		orthologs[o.sign.up]   <- +1
		orthologs[o.sign.down] <- -1
		orthologs[o.nonsign]   <- 0
		
		homeologs1 <- vector(length=n.genes)
		homeologs1[h1.sign.up]   <- +1
		homeologs1[h1.sign.down] <- -1
		homeologs1[h1.nonsign]   <- 0
		
		if(second.hybrid == TRUE){
			homeologs2 <- vector(length=n.genes)
			homeologs2[h2.sign.up]   <- +1
			homeologs2[h2.sign.down] <- -1
			homeologs2[h2.nonsign]   <- 0			
		}
		
		simulation1 <- matrix(data=cbind(orthologs, homeologs1), ncol=2)
		
		if(second.hybrid == TRUE){
			simulation2 <- matrix(data=cbind(orthologs, homeologs2), ncol=2)
		}
		
		# build change vectors (statistically significant and >=x-fold change)	
		orthologs.fc <- vector(length=n.genes)
		orthologs.fc[o.sign.up.fc]   <- +1
		orthologs.fc[o.sign.down.fc] <- -1
		orthologs.fc[o.nonsign.fc]  <- 0
		
		homeologs1.fc <- vector(length=n.genes)
		homeologs1.fc[h1.sign.up.fc]   <- +1
		homeologs1.fc[h1.sign.down.fc] <- -1
		homeologs1.fc[h1.nonsign.fc]  <- 0
		
		if(second.hybrid == TRUE){
			homeologs2.fc <- vector(length=n.genes)
			homeologs2.fc[h2.sign.up.fc]   <- +1
			homeologs2.fc[h2.sign.down.fc] <- -1
			homeologs2.fc[h2.nonsign.fc]  <- 0
		}
		
		simulation1.fc <- matrix(data=cbind(orthologs.fc, homeologs1.fc), ncol=2)
		
		if(second.hybrid == TRUE){
			simulation2.fc <- matrix(data=cbind(orthologs.fc, homeologs2.fc), ncol=2)
		}
		
		# calculate ortholog-homeolog expression classes
		outcome1    <- apply(simulation1, 1, classify)
		outcome1.fc <- apply(simulation1.fc, 1, classify)
		
		if(second.hybrid == TRUE){
			outcome2    <- apply(simulation2, 1, classify)
			outcome2.fc <- apply(simulation2.fc, 1, classify)		
		}
		
		# summarize results
		if(second.hybrid == FALSE){
			this.res[i,] <- c( conditions, sim.nde, sim.nde.fc, count(outcome1), count(outcome1.fc) )
		}else if(second.hybrid == TRUE){
			this.res[i,] <- c( conditions, sim.nde, sim.nde.fc, count(outcome1), count(outcome1.fc), count(outcome2), count(outcome2.fc) )
		}
			
	}
	
	return( this.res )
}
