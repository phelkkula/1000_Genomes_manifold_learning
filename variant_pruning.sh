#!/bin/bash

PLINK="path_to_plink"
HIGHLD="high-LD-regions.txt"
VCF="ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
NTHREADS=4

#Convert VCF into bfiles
$PLINK --vcf $VCF --threads $NTHREADS --make-bed --out input