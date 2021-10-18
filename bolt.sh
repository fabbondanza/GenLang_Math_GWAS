#!/bin/bash

###
# The current master script creates scripts to run a GWAs and SNP-heritability with BOLT-LMM and BOLT-REML respectively. The current script considers a single phenotype name which has to be passed as an argument
# Filippo Abbondanza
# 18.10.2021
###

pheno=$1
phenofile="/home/fabbonda/scratch/private/GWAS_Math/ALSPAC/math_alspac.sample"
path="/home/fabbonda/scratch/private/GWAS_Math"

set -e # With this command the script will be interrupted at the first error

cd $path/$pheno/ # The folder needs to be created in advance

echo -e "Make scripts for $pheno"

for i in {1..22}

do

cat > $pheno.chr${i}.sh <<EOF1
#!/bin/bash
#SBATCH --job-name="gwas"
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

~/scratch/apps/conda/envs/geno_utils/bin/bolt \
--bfile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/ARRAY_2019-05-13/all1/data/data \
--bgenFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/HRC_2019-05-13/all1/data/bgen_8bit/data_8bit_$i.bgen \
--sampleFile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/HRC_2019-05-13/all1/data/data.sample \
--bgenMinMAF=0.01 \
--phenoFile=$phenofile \
--phenoCol=$pheno \
--covarFile=$phenofile \
--covarCol=sex \
--lmm \
--LDscoresFile=/home/fabbonda/scratch/tools/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/home/fabbonda/scratch/tools/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
--numThreads=8 \
--statsFile=results/$pheno.chr${i}.stats.gz \
--statsFileBgenSnps=results/$pheno.chr${i}.bgen.stats.gz \
--verboseStats

EOF1

chmod +x $pheno.chr${i}.sh
#sbatch $pheno.chr${i}.sh # This will depend on the cluster manager system


done

cat > $pheno.heritability.sh <<EOF2
#!/bin/bash

#SBATCH --job-name="heritability"
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

~/scratch/apps/conda/envs/geno_utils/bin/bolt \
--bfile=/mnt/shared/projects/uosa/Silvia_Paracchini/shared/shared_data/new_alspac_external_500G_4/genetics/ARRAY_2019-05-13/all1/data/data \
--phenoFile=${phenofile} \
--phenoCol=$pheno \
--covarFile=${phenofile} \
--covarCol=sex \
--reml \
--numThreads=8 \
--geneticMapFile=/home/fabbonda/scratch/tools/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
> $pheno.h2

EOF2

chmod +x $pheno.heritability.sh
#sbatch $pheno.heritability.sh

# The script below will be used to combine all the chromosome and 
cat > results/${pheno}.combine.results.sh <<EOF3
#!/bin/bash

#SBATCH --job-name="combine.results"
#SBATCH --mem=12G
#SBATCH --cpus-per-task=24

if [ -f sumstats.${pheno} ];
  then
    echo "Sumstat already exists"
  else
    for i in {1..22}
    do
      echo -e "Opening files: chr$i"
      gunzip ${pheno}.chr\${i}.bgen.stats.gz
      awk 'NR > 1' ${pheno}.chr\${i}.bgen.stats > ${pheno}.chr\${i}.temp
      gzip ${pheno}.chr\${i}.bgen.stats
    done

    echo "Concatenating files"
    cat ${pheno}.chr1.temp ${pheno}.chr2.temp ${pheno}.chr3.temp ${pheno}.chr4.temp ${pheno}.chr5.temp ${pheno}.chr6.temp ${pheno}.chr7.temp ${pheno}.chr8.temp ${pheno}.chr9.temp ${pheno}.chr10.temp ${pheno}.chr11.temp ${pheno}.chr12.temp  ${pheno}.chr13.temp ${pheno}.chr14.temp ${pheno}.chr15.temp ${pheno}.chr16.temp ${pheno}.chr17.temp ${pheno}.chr18.temp ${pheno}.chr19.temp ${pheno}.chr20.temp ${pheno}.chr21.temp ${pheno}.chr22.temp > sumstats.${pheno}

    rm  *.temp
fi

# EXCLUDE FOR NOW

# Make file for Manhattan plot
# 1: SNP
# 2: CHROM
# 3: BP
# 14: P_BOLT_LMM_INF

# Header for sumstats
sed -i $'1 i\\\ID\tCHR\tBP\tGENPOS\tALLELE1\tALLELE0\tA1FREQ\tINFO\tCHISQ_LINREG\tP_LINREG\tBETA\tSE\tCHISQ_BOLT_LMM_INF\tP_BOLT_LMM_INF' sumstats.${pheno}

if [ -f manhattan.${pheno}.gz ];
  then
    echo "Manhattan file already exists"
  else
    echo "Preparing file for Manhattan plot"
    cut -d$'\t' -f 1,2,3,5,6,7,11,12,14 sumstats.${pheno} > manhattan.${pheno}
    sed -i $'1 i\\\ID\tCHR\tBP\tP'  manhattan.${pheno}

    gzip manhattan.${pheno}
fi

echo "Finished"

EOF3

chmod +x results/${pheno}.combine.results.sh
cd results
#sbatch ${pheno}.combine.results.sh

cd $path

