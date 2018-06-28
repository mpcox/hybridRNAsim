# hybridRNAsim

### Simulating 'neutral' gene expression patterns in hybrids and their progenitors

When comparing genome-wide patterns of gene expression in hybrid species vs their parents, the changes can be divided into four main classes. 

1. The homeologs match the pattern observed in the parents (**I**, inheritance). 
2. The homeologs show a difference that is not seen in the parents (**S**, bias). 
3. The parents have an expression difference, but this is not seen in the homeologs (**L**, blending). 
4. The parents have an expression difference, and the opposite is seen in the homeologs (**R**, reversal).

These patterns are described in more detail in this paper:

Cox, M.P., T. Dong, G. Shen, Y. Dalvi, D.B. Scott and A.R.D. Ganley. 2014. An interspecific fungal hybrid reveals cross-kingdom rules for allopolyploid gene expression patterns. *PLoS Genetics* 10: e1004180.
[https://doi.org/10.1371/journal.pgen.1004180](https://doi.org/10.1371/journal.pgen.1004180)

An interesting question is whether these patterns are due to selection on functional changes, or are instead just neutral. That is, when you compare two random sets of expression profiles (such as with gene expression from two sets of homeologs), you are always going to see some distribution of I, S, L and R classes just through chance effects.  If the observed distribution does not look different from random simulations of this process, then you should favor neutrality over selection. Alternately, by looking at how the observed and simulated distributions differ, you can estimate what proportion of genes follow the neutral pattern vs possibly being under selection.  (You cannot tell which genes are under selection using this information, of course).

These simulations are all done in R, leveraging the Bioconductor package, [compcodeR](https://bioconductor.org/packages/release/bioc/html/compcodeR.html).  I have generalized the code, which contains the main run commands in the primary file (hybridRNAsim.R).  For clarity, the more complex bits of code are kept in a separate file of functions (hybridRNAsim_functions.R).

Note that the code is not parallelized.  (I tried, but for various reasons, it is not as easy as you might think it should be).  The code is also moderately slow (i.e., for some analyses, it may run for a couple of days).  It is, of course, well suited to being run in an ‘embarrassingly parallel’ setting (i.e., multiple manual sub-jobs).

Like most simulations, it is necessary to condition on the study conditions as much as possible.  As written, the code currently uses the differential expression engine of [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), so all parent ortholog and hybrid homeolog differences that are given to the simulation code should be estimated with the same package. (Other differential expression packages can be subbed in, but DESeq2 seems like a reasonable one for right now).

The key parameters that need to be given to the simulation code are:

-    The number of iterations (i.e., the number of simulation tests you want to perform)
-    The number of genes (which needs to be the same for parents and hybrids)
-    The number of biological replicates (also the same for parents and hybrids)
-    The number of differentially expressed genes between the parent orthologs
-    The fraction of orthologs that are upregulated (vs downregulated)
-    The number of sequence reads in the parent comparison (say, mean of both parents)
-    The number of differentially expressed genes between the hybrid homeologs
-    The fraction of homeologs that are upregulated (vs downregulated)
-    The number of sequence reads in the homeolog comparison (say, the mean)

You can also set:

-    The significance cut-off value (default: 0.05)
-    The fold-change cut-off (default: 2x)

There is also a variable called ‘scalar’.  Because the code draws from a random negative binomial, there is no need to generate an entirely new dataset for every simulation (slow). Instead, you can make a very large number of calls from a random negative binomial, and then randomly draw from this distribution to make each simulation (fast).  The scalar variable sets the size of this distribution by multiplying the number of genes.  The current setting (default: 10x) seems to work well.

All input variables can be controlled as loops (e.g., expression difference in the parents ranging from 0-60% in 2% steps, with 2, 5, 10 or 25 biological replicates).

The output is given on a single row, with multiple rows listing different simulations.  For tracking purposes, all of the input variables are reported in the output.  For sanity checking, the number of up-regulated and down-regulated genes are reported (and should match the input values, within simulation error).  The I, S, L and R classes are also given.

The output variables are given in two forms.  First, all genes that are significantly different (e.g., ‘o.up’).  Second, all genes that are significantly different and match the fold change threshold (e.g., ‘o.up.fc’).  Depending on the question you’re asking, one of these outputs may be more useful than another.  Note: all the input variables assume that the genes are only significantly different and do not have a fold change threshold.  (Pay attention to this, as these different definitions of ‘differentially expressed’ are confusing and constantly get me into trouble).

The following figure shows the example suggested above: expression difference in the parents ranging from 0-60% in 2% steps, with 2, 5, 10 or 25 biological replicates.  You can see some key trends:

-  When Inheritance is high, the other classes are near zero.  As Inheritance decreases, the other classes become more common.
-  You get more accurate estimates of expression difference with more biological replicates.  (We cannot change this for our experiments, but this figure shows how important it is to set that parameter right).

The obvious test to do for our data is to match the observed settings as closely as possible, then see if the number of observed genes in each class fall within the confidence intervals of these neutral simulations. A mock-up of this is shown in the following figure.

![Example Simulation Figure](https://github.com/mpcox/hybridRNAsim/blob/master/SimulationExample.jpg)

Murray Cox

22 May 2017
[Updated: 28 June 2018]
