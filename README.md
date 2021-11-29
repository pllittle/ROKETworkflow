# ROKET Workflow

The workflow steps to analyze TCGA cancer types are contained within this package.

## Analysis Steps

1. Install the workflow package

* `smartr` can be found [here](https://github.com/pllittle/smartr)
* `ROKET` can be found [here](https://github.com/pllittle/ROKET)
* `ggh4x` can be found [here](https://github.com/teunbrand/ggh4x)

```
# Dependencies
all_packs = as.character(installed.packages()[,1])
req_packs = c("data.table","Rcpp","smartr",
	"ggplot2","MiRKAT","readxl","ROKET",
	"GOSemSim","ggh4x","devtools")

for(pack in req_packs){
	if( pack %in% all_packs ) next
	stop(sprintf("Install R package = %s",pack))
}

# Main package
if( !("ROKETworkflow" %in% all_packs) )
	devtools::install_github("pllittle/ROKETworkflow")

```

Also clone the repo.

```
git clone https://github.com/pllittle/ROKETworkflow.git
```

2. Download data

```
library(ROKETworkflow)

# Pick a cancer type
cancer_type = "BLCA"

# change to the parent directory of cloned ROKETworkflow
git_dir = "."

# change to preferred working directory
work_dir = "."

# Create directory hierarchy
my_dirs = ROKETworkflow::setdirs(git_dir = git_dir,
	work_dir = work_dir,
	dataset = cancer_type)

# Run code below, follow instructions until no errors
ROKETworkflow:::get_REF(my_dirs = my_dirs)
ROKETworkflow:::down_PanCan(my_dirs = my_dirs)
ROKETworkflow:::down_SPMs(my_dirs = my_dirs)
```

3. Pre-process data

First, filter somatic point mutation and convert into gene mutation status (1 = mutated gene, 0 = unmutated gene) by sample matrices. Then calculate gene ontology (GO-based), canonical pathways (PATH-based), and mutual exclusivity (ME-based) gene similarity matrices for various mutation frequency thresholds (e.g. include genes mutated in at least 5 samples would require `min_gene_mut = 5`).

```
mSPM = ROKETworkflow:::make_SPM_mat(my_dirs = my_dirs)
dim(mSPM); mSPM[1:5,1:4]

tab = ROKETworkflow:::make_tb1(my_dirs = my_dirs)

# Thresholds used by the manuscript
gene_muts = c(2,5,10,15,20,30)

gsim = ROKETworkflow::run_full_geneSim(my_dirs = my_dirs,
	min_gene_mut = gene_muts[4])
dim(gsim); gsim[1:5,1:5,]
```

4. Calculate optimal transport distances

This code should be run on a cluster with multiple threads. This step could take several hours or days, depending on the distribution of gene mutation frequency. For one or multiple `LAMBDA` penalties, calculate optimal transport based distances based on the three gene similarities.

```
LAMBDAs = c(0.5, 1.0, 5.0, Inf)

# Increase this if more threads are available
ncores = 1

ROKETworkflow::calc_full_DIST(my_dirs = my_dirs,
	min_gene_mut = min_gene_mut,
	LAMBDA = LAMBDAs[1],
	ncores = ncores)
```

5. Determine null models per cancer type and outcome and run kernel regressions with hypothesis testing

```
# Specify how many permutations
nPERM = 1e5

# Run regressions
ROKETworkflow::ANA_TEST_final(my_dirs = my_dirs,
	nPERM = nPERM)
```

6. Summarize over multiple LAMBDAs, cancer types, gene mutation frequencies, clinical outcomes. Compare and contrast findings between Euclidean and optimal transport. The function below assumes all 17 cancer types were processed and analyzed at the 4 LAMBDA penalties, 6 gene mutation frequency thresholds.

```
ROKETworkflow::new_ANA_AGG(my_dirs = my_dirs)
```

