#Transfer of data to the PSU Cluster using Cyberduck program
Each file was placed in a separate folder corresponding to the sample ID

#Setting up conda environment and installing quality analysis tool FastQC

#FastQC analysis on all sequences
fastqc sample_*/*.fastq.gz

#Installing MultiQC into a conda environment(Python2.7-env) since it required Python version 2.7-3.6
multiqc fastqc/
#Generates a collated report of all fastqc analysis metrics into a single text file for a quick view
cat multiqc_fastqc.txt
#Observation: Several sample sequences failed QC checks in terms of per base sequence quality and other checks
#Installed quality processing tool Trimmomatic and the format for the analysis command is as follows

trimmomatic PE -phred33 ../../sequences/sample_number/sample_ID_R1.fastq.gz ../../sequences/sample_number/sample_ID_R2.fastq.gz  sample_ID_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#FastQC and MultiQC were run again to ensure appropriate quality for downstream analysis
multiqc ../../trimmed/

#hg38 reference genome file downloaded from UCSC genome browser
#STAR aligner was used to perform alignment
#Generation of STAR genome indices
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir star_indices --genomeFastaFiles ref/hg38.fa

bwa mem ../ref/hg38/hg38.fa ../../umi_removed_sample_1.fastq > sample_1_umi_removed.sam


cat sample_1_counts.txt | awk '($2!=0) {print}'

#Extraction of UMIs using UMI Tools
 umi_tools extract --stdin=rna_analysis/trimmed/sample_1/CR2000_S1_L001_trimmed.fastq --bc-pattern=NNNNNNNNN --log=processed.log --stdout umi_removed_sample_1.fastq



for a in nt.*.tar.gz; do tar xzf $a; done


#Using HTSeq to get feature counts
The appropriate Ensembl annotation file was downloaded from UCSC Genome browser

htseq-count -s yes -t exon ../aligned/sample_1/sample_1Aligned.out.sam ../ref/hg38.ensGene.gtf > sample_1_counts.txt

htseq-count -s yes -t exon ../bwa_index/sample_1_umi_removed.sam ../ref/hg38/hg38.ensGene.gtf > sample_1_counts.txt

htseq-count -s no ../bwa_index/sample_1_umi_removed.sam ../ref/hg38/hg38.ncbiRefSeq.gtf > sample_1_counts.txt











#Installation of PLINK and PRSice using conda and tarball respectively
#Tutorial from https://choishingwan.github.io/PRS-Tutorial/target/
#Quality control of Base data
#Step 1: Checking intactness of files
md5sum summary_stats/PGC_data/pgc-bip2021-all.vcf.tsv
Output:
6d00aab1094c69ee668f3c04edfc031f
Indicates that the file is uncorrupted and intact
#Step 2:
Checking the Genome Build
Base data: NCBI Build 37/UCSC hg19
Target data: NCBI Build 37
#Step 3:
Removing SNPs with MAF < 1% and INFO < 0.8 (with very large base sample sizes these thresholds could be reduced if sensitivity checks indicate reliable results)

cat pgc-bip2021-all.vcf.tsv |\
awk 'NR==1 || ($11 > 0.01) && ($10 > 0.8) {print}' |\
gzip  > bip_pgc_qual.gz

#Step 4:
Removal of duplicate SNPs

gunzip -c bip_pgc_qual.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip > bip_pgc_nodup.gz

#Step 5:
Removal of Ambiguous SNPs

gunzip -c bip_pgc_nodup.gz|\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
gzip > bip_pgc_qc.gz

#Quality control of Target data
#Step 1:
General QC checks: Removing SNPs with low genotyping rate, low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing individuals with low genotyping rate

plink --bfile target_data/RG_CR_Final_3 --maf 0.05 --mind 0.1 --geno 0.1 --hwe 1e-6 --make-just-bim --make-just-fam --write-snplist--out target.qc

#Step 2:
Removal of highly correlated SNPs

plink --bfile target_data/RG_CR_Final_3 --keep qc_results/target.qc.fam --extract qc_results/target.qc.snplist --indep-pairwise 200 50 0.25 --out target.QC

OUTPUT:
Total genotyping rate is 0.994837

#Step 3:
Compute heterozygosity

plink --bfile target_data/RG_CR_Final_3 --extract qc_results/target.QC.prune.in --keep qc_results/target.qc.fam --het --out target.QC

Output:
Total genotyping rate is 0.994194.
192701 variants and 118 people pass filters and QC

#Step 4:
Exclusion of samples with high heterozygosity rate based on F coefficient estimates

data <- read.table("/storage/home/sma6401/work/BD_mito/gwas/qc_results/target.QC.het",header=T)
m <- mean(data$F)
s <- sd(data$F)
filtered <- subset(data, F <= m+3*s & F >= m-3*s)
write.table(filtered[,c(1,2)], "BIP.PGC.Filtered.sample", quote=F, row.names=F)

# Step 5:
Removal of Mismatching SNPs

setwd("/storage/home/sma6401/work/BD_mito/gwas/summary_stats/PGC_data")
# Load the bim file, the summary statistic and the QC SNP list into R
# Read in bim file
bim <- read.table("/storage/home/sma6401/work/BD_mito/gwas/qc_results/target.qc.bim") 
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
# Read in QCed SNPs
qc <- read.table("/storage/home/sma6401/work/BD_mito/gwas/qc_results/target.qc.snplist", header = F, stringsAsFactors = F)
# Read in the GWAS data
base_data_qc_snp_mismatch_removal <- read.table(gzfile("/storage/home/sma6401/work/BD_mito/gwas/summary_stats/PGC_data/bip_pgc_qc.gz"),header = T,stringsAsFactors = F,sep="\t")
colnames(base_data_qc_snp_mismatch_removal) <- c("CHROM","POS","ID","A1","A2","BETA","SE","PVAL","NGT","FCAS","FCON","IMPINFO","NEFFDIV2","NCAS","NCON","DIRE")
library(dplyr)
# Identify SNPs that require strand flipping
# Merge summary statistic with target
info <- left_join(bim, base_data_qc_snp_mismatch_removal , by = c("CHR"= "CHROM","BP"="POS","SNP"="ID"))
# Filter QCed SNPs
info <- info[info$SNP %in% qc$V1,]
# Function for finding the complementary allele
complement <- function(x) {
    switch (
        as.character(x),
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
# Get SNPs that have the same alleles across base and target
info.match <- subset(info, "A1" == "B.A1" & "A2" == "B.A2")
# Identify SNPs that are complementary between base and target

info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
# Update the complementary alleles in the bim file
# This allow us to match the allele in subsequent analysis

complement.snps <- bim$SNP %in% info.complement$SNP

bim[complement.snps,]$B.A1 <- sapply(bim[complement.snps,]$B.A1, complement)
bim[complement.snps,]$B.A2 <- sapply(bim[complement.snps,]$B.A2, complement)

# identify SNPs that need recoding
info.recode <- subset(info, "A1" == "B.A2" & "A2" == "B.A1")

# Update the recode SNPs
recode.snps <- bim$SNP %in% info.recode$SNP
tmp <- bim[recode.snps,]$B.A1
bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
bim[recode.snps,]$B.A2 <- tmp

# identify SNPs that need recoding & complement
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update the recode + strand flip SNPs
com.snps <- bim$SNP %in% info.crecode$SNP
tmp <- bim[com.snps,]$B.A1

bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
# Output updated bim file
write.table(bim[,c("SNP", "B.A1")],"target.qc.a1",quote = F,row.names = F,col.names = F,sep="\t")

# Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)
mismatch <-
    bim$SNP[!(bim$SNP %in% info.match$SNP |
                bim$SNP %in% info.complement$SNP | 
                bim$SNP %in% info.recode$SNP |
                bim$SNP %in% info.crecode$SNP)]
                
write.table(mismatch,"bip_pgc.mismatch",quote = F,row.names = F,col.names = F)


#Remove duplicate SNPs
#Identify duplicated SNPs
sort -k3n target.qc.bim |  uniq  -f2 -D | cut -f2 > dupeSNP.txt
plink --noweb --bfile ../target_data/RG_CR_Final_3 --exclude dupeSNP.txt --make-bed --out target.qc

#Sex check on sex chromosomes
plink --bfile ../target_data/RG_CR_Final_3 --extract target.QC.prune.in --keep ../BIP.PGC.Filtered.sample --check-sex --out target.qc

the generated file target.qc.sexcheck containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.

#Remove mismatched sex information
data <- read.table("target.qc.sexcheck", header=T)
valid <- read.table("../BIP.PGC.Filtered.sample", header=T)
valid <- subset(data, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "target.qc.valid", row.names=F, col.names=F, sep="\t", quote=F)

#Assuming there is no sample overlap
#Checking for relatedness
plink --bfile ../target_data/RG_CR_Final_3 --extract target.QC.prune.in --keep target.qc.valid --rel-cutoff 0.125 --out target.qc
#Remove related individuals
12 individuals were removed

#Generate final QC'ed target data file
plink --bfile target_data/RG_CR_Final_3 --make-bed --keep qc_results/target.qc.rel.id --out target.qc --extract qc_results/target.qc.snplist --exclude target.qc.mismatch --a1-allele target.qc.a1
