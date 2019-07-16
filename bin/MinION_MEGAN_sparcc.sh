#!/bin/bash
#SBATCH -t 01:30:00

####################################################################
## SPARCC MICROBIOME NETWORKS
####################################################################
## 0. Create .conda environment for sparcc installation: https://bitbucket.org/yonatanf/sparcc/src/default/
#conda create --prefix $STORE/ENVS/conda/py_sparcc python=2.6.9
#conda install numpy=1.9.2 pandas=0.16.2
#conda clean -t # clean tarballs from pkgs (cache)
#conda clean -p # clean non-used packages from pkgs

## 1. SPARCC CORRELATION MATRIX
module load miniconda2/4.5.11 parallel
source activate py_sparcc

sparCCdir=$STORE/ENVS/yonatanf-sparcc-3aff6141c3f1
WDIR=$LUSTRE/Results/Metagenome/METALGEN/7_networks

# Input: .txt (tsv) table formatted in RStudio from rma_taxa_count_filtered.csv
python SparCC.py -h

# OPTIONS:
# -a allows to choose the type of correlation analysis to carry out: sparcc (default), pearson, spearman

python $sparCCdir/SparCC.py $WDIR/sparcctable.txt -i 5 --cov_file=$WDIR/basis_corr/cov_sparcc.out --cor_file=$WDIR/basis_corr/cor_sparcc.out

python $sparCCdir/MakeBootstraps.py $WDIR/sparcctable.txt -n 100 -t permutation_#.txt -p $WDIR/pvals/

for f in $(ls $WDIR/pvals/permutation*.txt); do
outname=$(echo $WDIR/pvals/perm_cor_${f##*permutation_})
sem -j 100 python $sparCCdir/SparCC.py $f -i 5 --cor_file=$outname
done

python $sparCCdir/PseudoPvals.py $WDIR/basis_corr/cor_sparcc.out $WDIR/pvals/perm_cor_#.txt 20 -o $WDIR/pvals/pvals.one_sided.txt -t one_sided
python $sparCCdir/PseudoPvals.py $WDIR/basis_corr/cor_sparcc.out $WDIR/pvals/perm_cor_#.txt 20 -o $WDIR/pvals/pvals.two_sided.txt -t two_sided

### 2. PCIT R PACKAGE: GETTING EDGELIST FOR CYTOSCAPE
## See https://academic.oup.com/bioinformatics/article/26/3/411/215002 and https://cran.r-project.org/web/packages/PCIT/PCIT.pdf
## For R table reading, we need to change '[SK] ' to 'SK.'
module load cesga/2018 gcc/6.4.0 R/3.5.1
R
library(PCIT)

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

WD <- "/mnt/lustre/scratch/home/otras/ini/alg/Results/Metagenome/METALGEN/7_networks/basis_corr/"

cormat <- read.delim(paste0(WD, "cor_sparcc.out"), dec = ".", stringsAsFactors = F, row.names = 1)
cormat <- data.matrix(cormat)
cormat2 <- round2(cormat, 10)
rownames(cormat2) <- colnames(cormat2)
isSymmetric(cormat2)
pcit <- pcit(cormat2)
signif <- idx(pcit)
nonsignif <- idxInvert(nrow(cormat2), signif)
cormat2[nonsignif] <- 0
cormatadj <- abs(cormat2)

edgel <- getEdgeList(cormatadj)
edgel$sign[edgel$Weight<0] <- '-'
edgel$sign[edgel$Weight>0] <- '+'
edgel$Weight <- abs(edgel$Weight)

write.table(edgel, file = paste0(WD,"../","edgelist.txt"), quote = F, row.names = F, sep = "\t")

quit()
n
module purge
