#!/usr/bin/env bash
# This is the pipeline script for Ext-BSA
# Bam filenames need to modify
# Sequencing depth pruning need to modify in variant calling

# Step0: Setting
DPlow=30   # lowest depth for variant pruning
DPhigh=200 # highest depth for variant pruning
REF=../Zea_mays.AGPv3.27.dna.genome.fa.mod.fa # define reference genome file


# Step1: Qulity control
#-----------------------------------------------------------------------------------
echo -e "\n-----------------------STEP 1------------------------------------\n"
echo -e "Starting quality control...\n"
echo -e "-----------------------------------------------------------------\n\n"

mkdir -p  01.QC && cd 01.QC
ln -sf ../*.gz  ./

echo -e "Fastqc processing...\n"

for file in *.gz

# FastQC process
do

  echo ${file}
  fastqc ${file} > ${file%%.*}.fastqc.log 2>&1

done

# Trimmomatic process
echo -e "\nTrimmomatic processing...\n"

for line in *.R1.clean.fastq.gz

do

  java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog trimlog.txt \
    ${line%%.*}.R1.clean.fastq.gz    ${line%%.*}.R2.clean.fastq.gz\
    ${line%%.*}.pair.R1.clean.fastq.gz  ${line%%.*}.unpair.R1.clean.fastq.gz\
    ${line%%.*}.pair.R2.clean.fastq.gz  ${line%%.*}.unpair.R2.clean.fastq.gz\
    LEADING:20  TRAILING:20  MINLEN:100  > ${line%%.*}.trimmomatic.log 2>&1

done

echo -e "Step 1: Quality control finished!\n\n"


# Step2: Sequence alignmnet
#-----------------------------------------------------------------------------------
echo -e "-----------------------STEP 2------------------------------------\n"
echo -e "Starting short PE read extension and sequence alignment...\n"
echo -e "-----------------------------------------------------------------\n\n"

cd ../ && mkdir -p  02.Alignment  && cd 02.Alignment
ln -sf ../ref ./
ln -sf ../01.QC/*.pair.R*.clean.fastq.gz   ./


# Flash to merge paired-end data
echo -e "FLASH processing...\n"

for file in *.pair.R1.clean.fastq.gz

do

  flash -m 30 -M 80 -O -r 150 -f 250 -z -o ${file%%.*}\
    ${file%%.*}.pair.R1.clean.fastq.gz  ${file%%.*}.pair.R2.clean.fastq.gz\
    > ${file%%.*}.flash.log 2>&1

done


# Check if alignment index exists
echo -e "BBMap alignment and Picard Markduplicates processing...\n"

if [ ! -d ref ]

  then

    echo -e "BBMap index start creation ...\n"
    bbmap.sh ref=${REF}

  else
    echo -e "BBMap index has existed!\n"

fi


for LINE in *.extendedFrags.fastq.gz

do

# Alignment with BBMap
  bbmap.sh  -Xmx100g  threads=10 in=$LINE  fast=t outm=${LINE%%.*}.sam  > ${LINE%%.*}.bbmap.log 2>&1

# Convert BAM and sort
  samtools sort -@ 10 ${LINE%%.*}.sam -o  ${LINE%%.*}.sorted.bam 

# Mark duplicate
  java  -jar   /home/zoucheng/bin/picard-2.10.10/picard.jar  MarkDuplicates I=${LINE%%.*}.sorted.bam \
    O=${LINE%%.*}.rmdup.sorted.bam METRICS_FILE=${LINE%%.*}.sort.metrics   REMOVE_DUPLICATES=true\
    ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT > ${LINE%%.*}.rmdup.log  2>&1

done

# Merge and sort bams belonging to the same pool
# This part need to change according to your files
samtools merge -@ 10 K.rmdup.merge.bam  K_L4_I334.rmdup.sorted.bam  K_L4_I335.rmdup.sorted.bam 
samtools merge -@ 10 G.rmdup.merge.bam  G_L5_I330.rmdup.sorted.bam  G_L5_I331.rmdup.sorted.bam

echo -e "Add read group processing...\n"

i=1

for LINE in G.rmdup.merge.bam  K.rmdup.merge.bam;

do

# Add read group
  bamaddrg\
     -b  ${LINE%%.*}.rmdup.merge.bam\
     -s  ${LINE%%.*}\
     -r  group"$i" > ${LINE%%.*}.rmdup.addrg.sort.bam

# Index bam file 
  samtools index ${LINE%%.*}.rmdup.addrg.sort.bam
  i=`expr $i + 1`

done
echo -e "Step 2: Sequence alignment finished!\n\n"

# Step3: Variant calling
#-----------------------------------------------------------------------------------
echo -e "-----------------------STEP 3------------------------------------\n"
echo -e "Variant calling for merged BAM files...\n"
echo -e "-----------------------------------------------------------------\n\n"

# Soft link reference genome
ln -sf ${REF}  ./

# Create bam list file
ls  *.rmdup.addrg.sort.bam > bam_list

# Index fasta reference file
samtools faidx Zea_mays.AGPv3.27.dna.genome.fa.mod.fa

# Call variants using bam
freebayes -f Zea_mays.AGPv3.27.dna.genome.fa.mod.fa -i -X -u  -K -L bam_list > F00.freebayes.BSA.vcf


# Remove missing genotypes in BSA(bulked sample analysis)
awk '$10 != "." && $11 != "."'  F00.freebayes.BSA.vcf > F01.ensure.gt.vcf

# Generate text format for next step
bcftools  query\
     -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%DP[\t%SAMPLE;%GT;%DP;%AD]\n' F01.ensure.gt.vcf\
     > F02.bcftools.format.vcf.txt

# Remove scaffold, Mt, Pt lines and calculate allele frequency
# Define the total sequencing depth threshold

egrep -v "scaffold|Mt|Pt" F02.bcftools.format.vcf.txt |\
     awk -v dplow=${DPlow} -v dphigh=${DPhigh} '$5 >= dplow && $5 <= dphigh' |\
     sed 's/[,;]/\t/g' |\
     awk '{printf"%s\t%s\t%s\t%s\t%s:%s\t%.4f\t%s:%s\t%.4f\n",$1,$2,$3,$4,$6,$7,$10/$8,$11,$12,$15/$13}'\
     > F03.allele.frequency.txt

# Format allele frequency text for Rscript testing
cut -f1,2,6,8 F03.allele.frequency.txt | sort -n -k 1 -k 2 > F04.ratio.for.test.sort.txt

echo -e "Step 3: Variant calling and allele frequency calculation have finished!\n\n"

# Step4: Ext-BSA testing
#-----------------------------------------------------------------------------------
echo -e "-----------------------STEP 4------------------------------------\n"
echo -e "Starting Ext-BSA...\n"
echo -e "-----------------------------------------------------------------\n\n"

cat << EOF > Ext-BSA.script.R
## Ext-BSA -- step 1
## read file from stdin
## load qqman package for plot

Args <- commandArgs()
file <- Args[6]
require(qqman) 
t <- read.table(file, header = F)
ouf_data <- "Ext_BSA_extension_result.txt"
ouf_pdf <- "Ext_BSA_extension_result.pdf"

## Ext-BSA -- step 2
## define data format and significant level

diffr <- abs(t\$V3 - t\$V4)
s.diffr <- sort(diffr)
r1 <- rank(diffr)
y <- -log10(1 - (r1-0.01)/length(r1))
index.y <- which(r1 > length(r1)*0.99)
threshold0 <-  -log10(1 -  length(r1)*0.99/length(r1))
print(threshold0)

## Ext-BSA -- step 3
## construct test function and run it

walking <- function(start, step = 5, threshold = threshold0) {
  score = y[start]
  left = list()
  left[[1]] = c(start,score)
  i = 1
  while(score > threshold) {
    tinterval = y[(start - step*i):(start - step*(i - 1) - 1)]
    score = score + sum(tail(sort(tinterval), 2)) - step
    i = i+1
    inp = c(start - step*(i - 1), score)
    left[[i]] = inp
  }
  
  start_score = y[start]
  right = list()
  right[[1]] = c(start,score)
  i = 1
  while(score > threshold) {
    tinterval = y[(start + step*(i - 1) + 1):(start + step*i)]
    score = score + sum(tail(sort(tinterval), 2)) - step
    i = i + 1
    inp = c(start-step*(i - 1), score)
    right[[i]] = inp
    
  }
  return(c(left, right))
}

index.y <- subset(index.y, index.y > 5)
lout <- lapply(index.y, walking)
outm <- matrix(unlist(lout), ncol = 2, byrow = T)
colnames(outm) <- c("ind", "value")
df <- data.frame(outm)
df2 <- data.frame( t[df\$ind, 2], t[df\$ind, 1], t[df\$ind, 2], df\$value )
colnames(df2) <- c("SNP", "CHR", "BP", "P")


## Ext-BSA -- step 4 
# part 1 -- remove replicate
final_result <- aggregate(df2, by = list(df2\$CHR, df2\$BP), FUN = max)[,3:6]
write.table(final_result, file = ouf_data, col.names = T, row.names = F, quote = F, sep = "\t")

# part 2 -- plot
cutoff <- 0.05/(dim(df2)[1])
pdf(file = ouf_pdf, width = 20, height = 5)
manhattan(df2, logp = F, genomewideline = -log10(cutoff), suggestiveline = F)
dev.off()

## T-test -- step5
## bulked sample number
n1 <- 30
n2 <- 30

raw.data <- t
t.stat  <- matrix(NA, ncol=1, nrow=nrow(raw.data))
p.value <- matrix(NA, ncol=1, nrow=nrow(raw.data))

## Calculation formula
# ExpFreq=(ObsFreq1*SampleSize1+ObsFreq2*SampleSize2)/&
#   (SampleSize1+SampleSize2)
# 
# TTest=(ObsFreq1-ObsFreq2)/&
#   SQRT((1.0/(2.0*SampleSize1)+1.0/(2.0*SampleSize2))*&
#          ExpFreq*(1.0-ExpFreq))

# Calculate the t.stat
raw.data[,5] <- (raw.data[,3] * n1 + raw.data[,4] * n2)/(n1 + n2)
t.stat[, 1] <- (raw.data[,3] - raw.data[,4])/sqrt((1.0/(2.0 * n1) + 1.0/(2.0 * n2)) * raw.data[,5] * (1.0 - raw.data[,5]))

## Calculate the pvalue and define the DF of Test 

p.value <- 2 * pt(abs(t.stat),  n1 - 1, lower.tail = F)
t_result <- data.frame(raw.data[,c(2,1:2)], p.value)
colnames(t_result) <- c("SNP", "CHR", "BP", "P")
write.table(t_result, file = "Ext_BSA_T_test.txt", quote=F, row.names=F, col.names = T, sep = "\t")
# plot T-test results
pdf(file="Ext_BSA_T_test_result.pdf", width = 20, height = 5)
cut <- -log10(0.05/(dim(t_result)[1]))
manhattan(t_result, genomewideline=F, suggestiveline=cut)
dev.off()

EOF

# Run R script

Rscript  Ext-BSA.script.R  F04.ratio.for.test.sort.txt  

echo -e "Step 4: Ext-BSA finished!\n\n"
echo -e "The pipeline finished!\n\n"
#-----------------------------------------------------------------------------------
