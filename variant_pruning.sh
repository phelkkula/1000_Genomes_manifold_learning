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