#!/bin/bash

PLINK="path_to_plink"
HIGHLD="high-LD-regions.txt"
VCF="ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
NTHREADS=4

#Convert VCF into bfiles
$PLINK --vcf $VCF --threads $NTHREADS --make-bed --out input

#Filter rsids; remove EBI structural variants (esv)
awk '{print $2}' input.bim | grep -v 'rs' > non_rsids.txt
$PLINK --bfile input --exclude non_rsids.txt --threads $NTHREADS --make-bed --out input

#Exclude high LD regions
$PLINK --bfile input --exclude range $HIGHLD --make-bed --threads $NTHREADS --out temp1

#Filter variants with call rate >= 0.9 and MAF >= 5%
$PLINK --bfile temp1 --geno 0.1 --maf 0.05 --threads $NTHREADS --make-bed --out temp2

#--indep-pairwise [window size]<kb> [step size (variant ct)] [r^2 threshold] produces
#a pruned subset of markers that are in approximate linkage equilibrium
$PLINK --bfile temp2 --indep-pairwise 50 5 0.2 --threads $NTHREADS --out temp2

#Extract pruned variants
$PLINK --bfile temp2 --extract temp2.prune.in --threads $NTHREADS --make-bed --out pruned

#Get allelic dosages for samples
$PLINK --bfile pruned --recodde A --threads $NTHREADS --out pruned