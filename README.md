# AtaCNV
This is the R package implementation for AtaCNV, a CNV detection algorithm for scATAC-seq datasets The method is presented in paper "Detecting copy number variations from single-cell chromatin accessibility sequencing data by AtaCNV". 

# Workflow

![Workflow](example/workflow.png)

AtaCNV takes a single-cell read count matrix over genomic bins of 1 million base pairs (1 Mbp) as input. Firstly, cells and genomic bins are filtered according to bin mappability and the number of zero entries. To reduce the extreme noisiness, AtaCNV then smooths the count matrix by fitting a one-order dynamic linear model for each cell. AtaCNV also performs cell-wise local regression to remove the potential biases caused by GC content. If normal cells are available, AtaCNV normalizes the smoothed count data against those of normal cells to deconvolute copy number signals from other confounding factors like chromatin accessibility. If normal cells are not available, observing that tumor single-cell data often contain many non-tumor cells, AtaCNV clusters the cells and identifies a group of high confidence normal cells and normalizes the data against their smoothed depth data. Then, AtaCNV applies the multi-sample BIC-seq algorithm to jointly segment all single cells and estimates copy ratios of the obtained segments for each cell. Using the copy ratio data, AtaCNV calculates a burden score for each cell group and classifies cell groups with high CNV burden as malignant cells. For tumor cells, AtaCNV further infers their discrete copy number states using a Bayesian method.

# Installation
```
install.packages("devtools")
devtools::install_github("Aelita-Stone/AtaCNV")
```

# Usage
Based on the workflow, the R package primarily includes three modules: normalization, segmentation, and copy number state inference.
## Normalization
Using count matrix as input, AtaCNV removes confounding factors like GC content to deconvolute CNV signals.
'''
norm_re <- AtaCNV::normalize(count,
                             genome = "hg19", 
                             mode = "normal cells",
                             normal_cells = cell_type=="normal",
                             output_dir = "./",
                             output_name = "norm_re.rds"
                             gc_correction = TRUE
                             )
'''
## Segmentation
AtaCNV utilizes BICseq to infer potential breakpoints and calculate segment-wise copy ratios.
'''
seg_re <- AtaCNV::calculate_CNV(norm_count = norm_re$norm_count,
                                baseline = norm_re$baseline,
                                output_dir = "./",
                                output_name = "CNV_re.rds")
'''
## Copy state inference
For tumor cells, AtaCNV further infers their discrete copy number states.
'''
CN_state <- AtaCNV::estimate_cnv_state_cluster(count = count,
                                               genome = "hg19",
                                               copy_ratio = norm_re$copy_ratio,
                                               bkp = seg_re$bkp,
                                               label = cell_type)
'''
More detailed instructions can be found in the R package documentation and the examples in {tutorial.R}.

