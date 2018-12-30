# hybridRNAsim

When comparing genome-wide patterns of gene expression in hybrid species vs their parents, the changes can be divided into four main classes:

1. The homeologs match the pattern observed in the parents (**I**, inheritance). 
2. The homeologs show a new difference that is not seen in the parents (**S**, bias). 
3. The parents have an expression difference, but this lost in the homeologs (**L**, blending). 
4. The parents have an expression difference, but the opposite is seen in the homeologs (**R**, reversal).

These patterns are described in more detail in this paper:

Cox, M.P., T. Dong, G. Shen, Y. Dalvi, D.B. Scott and A.R.D. Ganley. 2014. An interspecific fungal hybrid reveals cross-kingdom rules for allopolyploid gene expression patterns. *PLoS Genetics* 10: e1004180.
[https://doi.org/10.1371/journal.pgen.1004180](https://doi.org/10.1371/journal.pgen.1004180)

The question is whether these patterns occur due to selection on functional changes, or are instead just a neutral outcome. That is, when you compare two random sets of expression profiles, you always see some distribution of I, S, L and R classes just through chance effects.  Thus, if the observed distribution is not statistically different from simulations of random homeolog creation, then neutrality should be favored over selection. Further, by looking at how the observed and simulated distributions differ, it is possible to estimate what proportion of genes follow the neutral pattern vs divergent from it (and thus possibly under selection).  As this is a statistical test based on random simulations, it is not possible to tell exactly which genes are under selection using the information from this model.

### Description

The simulation code is written in R, and leverages the Bioconductor package, [compcodeR](https://bioconductor.org/packages/release/bioc/html/compcodeR.html).  The main run commands of the code, which is generalized as much as possible, is contained in the primary code file (hybridRNAsim.R).  For clarity, all functions and the more complex parts of the code are placed in a separate file (hybridRNAsim_functions.R).

Note that the code is not parallelized.  The simulations are well suited to being run in an ‘embarrassingly parallel’ setting (i.e., multiple manual sub-jobs).  However, other parallelization was inefficient, due - I think - to the nature of the compcodeR package.  Consequently, the simulations are moderately slow; for some analyses, the code may run for a couple of days.

### Running

As with most simulations, it is important to condition on the study conditions as much as possible.  As written, the code currently uses the differential expression engine of [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), so all parent ortholog and hybrid homeolog differences that are fed into the simulation code should be estimated with DESeq2 as well. (Other differential expression packages can be subbed into the code, but DESeq2 seems like a reasonable option for now).

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

If you wish to simulate patterns in two hybrids, you must set second.hybrid=TRUE (default: FALSE).  You then also need to assign values to the following parameters:

-    The number of differentially expressed genes between the homeologs in hybrid 2
-    The fraction of homeologs that are upregulated (vs downregulated) in hybrid 2
-    The number of sequence reads in the homeolog comparison (say, the mean) in hybrid 2

There is also a ‘scalar’ variable.  Because the code draws from a random negative binomial, there is no need to generate an entirely new dataset for every simulation (slow). Instead, a very large number of calls can be made from a random negative binomial, and then a subset can be randomly draws from this distribution to create each simulation (fast).  The scalar variable sets the size of this distribution by multiplying the number of genes.  The current setting (default: 10x) seems to work well.  (Note that if the scalar value is 1, all random subsets will necessarily be identical).

All input variables can be controlled as loops (e.g., expression difference in the parents ranging from 0-60% in 2% steps, with 2, 5, 10 or 25 biological replicates).

### Output

Each simulation output is given on a single row, with multiple rows listing the results from different simulations.  For tracking purposes, all of the input variables are reported in the output.  For sanity checking, the number of up-regulated and down-regulated genes are reported (and should match the input values, within simulation error).  The I, S, L and R classes are also given.

The output results are listed in two forms.  First, all genes that are significantly different (e.g., ‘o.up’).  Second, all genes that are significantly different *and* match the fold change threshold (e.g., ‘o.up.fc’).  Depending on the question being asked, one of these outputs may be more useful than the other.  Note: all of the *input* variables assume that the genes are only significantly different, but do not have a fold change threshold.  Pay attention to this, as these different definitions of ‘differentially expressed’ are confusing and a constant source of trouble (at least to me...)

If two hybrids are being simulated, the output results are doubled with 1 appended to values for the first hybrid and 2 for the second hybrid.

### Worked examples

Worked examples for one and two hybrids can be found in the file hybridRNAsim.R

### Example figure

The following figure shows the loop example suggested above: expression difference in the parents ranging from 0-60% in 2% steps, with 2, 5, 10 or 25 biological replicates.  Some key trends can be seen:

-  When Inheritance is high, the other classes are near zero.  As Inheritance decreases, the other classes become more common.
-  More accurate estimates of expression difference are found with more biological replicates.  (This variable often cannot be changed, but this figure shows how important it is to set that parameter right).

The obvious test to perform is to match the observed empirical setting as closely as possible, and then see whether the number of genes observed in each class fall within the confidence intervals of these 'neutral' simulations.

![Example Simulation Figure](https://github.com/mpcox/hybridRNAsim/blob/master/SimulationExample.jpg)

Murray Cox

22 May 2017
[Updated: 2 July 2018]
