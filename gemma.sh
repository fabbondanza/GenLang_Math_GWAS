#!/bin/bash

###
# The current master script creates scripts to run a GWAS with GEMMA. It is based on the work from Else Eising. Currently it does not handle chrX
# Filippo Abbondanza
# 18.10.2021
###

#---------------------------------------
# Defining directories
#---------------------------------------

VCFdir="/home/fabbonda/scratch/private/York_Else"
WorkingDir="/home/fabbonda/scratch/private/GWAS_Math/YORK"
GenotypeDir="/home/fabbonda/scratch/private/York_Else"

#---------------------------------------
# Make scripts to submit to cluster
#---------------------------------------

source activate geno_utils # This depends if you have a cond aenv with all the tools you need

echo "Prep scripts"

cd $WorkingDir

cat > Script_template_York_Gemma_association_analysis.txt <<EOF1
#!/bin/bash
#SBATCH --job-name="York_gwas"
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G

#---------------------------------------
# calculate relatedness
#---------------------------------------

source activate geno_utils

# For the calculation of the relation matrix (using all autosomes)
# use genotyped, pruned data
# for each chromosome, make a 'leave-1-chromosome-out' SNP list for the GRM

# Script expects the following folder: $WorkingDir/SNPlist_for_GRMs

cd $WorkingDir

echo 'Create leave-one-out data'
awk '\$1!=__CHR__ {print \$2}' $GenotypeDir/york.full.autosomes.pruned.bim > $WorkingDir/SNPlist_for_GRMs/SNPlist_chr__CHR__.txt

if [ -f output/York_all_relatedness_chr__CHR__ ];
  then
    echo -E "Relatedness matrix for chr$i already exists"
  else
  echo 'Calculate relatedness matrix'

  # use the genotyped, pruned summary statistics to calculate the relatedness matrix
  gemma \
  -bfile $GenotypeDir/york.full.autosomes.pruned \
  -gk 1 \
  -snps $WorkingDir/SNPlist_for_GRMs/SNPlist_chr__CHR__.txt \
  -o York_all_relatedness_chr__CHR__.LOCO
fi

#---------------------------------------
# Association analysis for total cohort
#---------------------------------------

# Gemma usage:
# -a -> specifies annotations file for LOCO
# -k -> specifies the relatedness matrix
# -n -> specifies the number of column to use for the analysis

# Examples bimbam (which seems to have some parsing issues)
#-g $GenotypeDir/York_chr__CHR__.plink.bimbam \

# Note: I've removed header from math_york.sample so that thte wc -l matches the one from GREML matrix


#Phenotype 1. In the -n specify the column of Phenotype 1
gemma \
-bfile $GenotypeDir/fulldataset.plink.chr__CHR__ \
-p $WorkingDir/math_york.sample \
-n 4 \
-a $GenotypeDir/York_chr__CHR__.plink.v2.map \
-k $WorkingDir/output/York_all_relatedness_chr__CHR__.cXX.txt \
-lmm 1 \
-o OMA_z_age_adj_chr__CHR__.LOCO \
-loco __CHR__

#Phenotype 2. In the -n specify the column of Phenotype 2
gemma \
-bfile $GenotypeDir/fulldataset.plink.chr__CHR__ \
-p $WorkingDir/math_york.sample \
-n 5 \
-a $GenotypeDir/York_chr__CHR__.plink.v2.map \
-k $WorkingDir/output/York_all_relatedness_chr__CHR__.cXX.txt \
-lmm 1 \
-o OMS_z_age_adj_chr__CHR__.LOCO \
-loco __CHR__

EOF1

chmod 777 Script_template_York_Gemma_association_analysis.txt

# Generate a job per sample nd submit job
for i in {1..22}
  do
  echo -e "Generate script for chr$i"
  sed "s/__CHR__/$i/g" Script_template_York_Gemma_association_analysis.txt > Script_York_Gemma_association_analysis_chr${i}_LOCO.sh
  chmod 777 Script_York_Gemma_association_analysis_chr${i}_LOCO.sh
  sbatch Script_York_Gemma_association_analysis_chr${i}_LOCO.sh
done

