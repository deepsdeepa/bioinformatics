# #!/bin/bash

NORMAL_FASTQ_R1="raw_data/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq.gz"
NORMAL_FASTQ_R2="raw_data/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq.gz"
CANCER_FASTQ_R1="raw_data/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq.gz"
CANCER_FASTQ_R2="raw_data/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq.gz"
REFERENCE_GENOME="/home/admin/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"   

OUTPUT_DIR="output"

mkdir -p $OUTPUT_DIR

# # Step 1: Quality Control using FastQC
echo "Running FastQC for Quality Control..."
mkdir -p $OUTPUT_DIR/fastqc_results
fastqc $NORMAL_FASTQ_R1 $NORMAL_FASTQ_R2 -o $OUTPUT_DIR/fastqc_results
fastqc $CANCER_FASTQ_R1 $CANCER_FASTQ_R2 -o $OUTPUT_DIR/fastqc_results
echo "FastQC analysis completed. Check $OUTPUT_DIR/fastqc_results for details."

# Step 2: Trim low quality sequences and adapter sequences using Cutadapt or Fastp
echo "Trimming adapter sequences and low-quality bases..."
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $OUTPUT_DIR/normal_trimmed_R1.fastq.gz -p $OUTPUT_DIR/normal_trimmed_R2.fastq.gz $NORMAL_FASTQ_R1 $NORMAL_FASTQ_R2
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $OUTPUT_DIR/cancer_trimmed_R1.fastq.gz -p $OUTPUT_DIR/cancer_trimmed_R2.fastq.gz $CANCER_FASTQ_R1 $CANCER_FASTQ_R2
echo "Adapter trimming completed."

# Step 3: Re-run FastQC after trimming to check the quality improvement
echo "Re-running FastQC after trimming..."
fastqc $OUTPUT_DIR/normal_trimmed_R1.fastq.gz $OUTPUT_DIR/normal_trimmed_R2.fastq.gz -o $OUTPUT_DIR/fastqc_results
fastqc $OUTPUT_DIR/cancer_trimmed_R1.fastq.gz $OUTPUT_DIR/cancer_trimmed_R2.fastq.gz -o $OUTPUT_DIR/fastqc_results
echo "Re-run FastQC analysis completed. Check $OUTPUT_DIR/fastqc_results for details."

# # Step 4: Alignment using BWA
# echo "Aligning trimmed samples to the reference genome using BWA..."
bwa mem $REFERENCE_GENOME $OUTPUT_DIR/normal_trimmed_R1.fastq.gz $OUTPUT_DIR/normal_trimmed_R2.fastq.gz > $OUTPUT_DIR/normal_sample.sam
bwa mem $REFERENCE_GENOME $OUTPUT_DIR/cancer_trimmed_R1.fastq.gz $OUTPUT_DIR/cancer_trimmed_R2.fastq.gz > $OUTPUT_DIR/cancer_sample.sam

# Convert SAM to BAM, sort and index BAM files
samtools view -bS $OUTPUT_DIR/normal_sample.sam | samtools sort -o $OUTPUT_DIR/normal_sample_sorted.bam
samtools view -bS $OUTPUT_DIR/cancer_sample.sam | samtools sort -o $OUTPUT_DIR/cancer_sample_sorted.bam

samtools index $OUTPUT_DIR/normal_sample_sorted.bam
samtools index $OUTPUT_DIR/cancer_sample_sorted.bam
echo "Alignment completed."

# Step 5: Create FASTA dictionary for reference genome using Picard
echo "Creating FASTA dictionary for the reference genome..."
picard CreateSequenceDictionary \
  R=$REFERENCE_GENOME \
  O=/home/admin/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.dict
    
  echo "FASTA dictionary created."

# Step 6: Add read groups using Picard
echo "Adding read groups to BAM files..."
picard AddOrReplaceReadGroups \
 I=$OUTPUT_DIR/normal_sample_sorted.bam \
  O=$OUTPUT_DIR/normal_sample_sorted_with_RG.bam \
    RGID=id1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=normal_sample

picard AddOrReplaceReadGroups \
  I=$OUTPUT_DIR/cancer_sample_sorted.bam \
  O=$OUTPUT_DIR/cancer_sample_sorted_with_RG.bam \
  RGID=id2 RGLB=lib2 RGPL=illumina RGPU=unit2 RGSM=cancer_sample
echo "Read groups added."

# Step 7: Somatic Mutation Calling using Mutect2 (GATK)
echo "Running Mutect2 for somatic mutation calling..."
samtools faidx $REFERENCE_GENOME
samtools index $OUTPUT_DIR/normal_sample_sorted_with_RG.bam
samtools index $OUTPUT_DIR/cancer_sample_sorted_with_RG.bam

mkdir -p $OUTPUT_DIR/mutect2
gatk Mutect2 \
  -R $REFERENCE_GENOME \
  -I $OUTPUT_DIR/normal_sample_sorted_with_RG.bam \
  -I $OUTPUT_DIR/cancer_sample_sorted_with_RG.bam \
  -O $OUTPUT_DIR/mutect2/somatic_mutations.vcf
echo "Mutect2 mutation calling completed."

# Step 8: Calculate Background Mutation Level in Normal Tissue
echo "Calculating background mutation level in the normal tissue..."
bcftools mpileup -f $REFERENCE_GENOME $OUTPUT_DIR/normal_sample_sorted_with_RG.bam | \
    bcftools call -mv -Ov -o $OUTPUT_DIR/normal_sample_variants.vcf

TOTAL_VARIANTS=$(grep -v '^#' $OUTPUT_DIR/normal_sample_variants.vcf | wc -l)
TOTAL_READS=$(samtools flagstat $OUTPUT_DIR/normal_sample_sorted_with_RG.bam | grep "in total" | awk '{print $1}')
MUTATIONS_PER_MILLION=$(echo "scale=2; $TOTAL_VARIANTS / $TOTAL_READS * 1000000" | bc)

echo "Background mutation level: $MUTATIONS_PER_MILLION mutations per million reads"

# Step 9: Detect somatic mutations
echo "Detecting somatic mutations present in cancer sample but not in normal sample..."
grep -v '^#' $OUTPUT_DIR/normal_sample_variants.vcf | awk '{print $3}' > $OUTPUT_DIR/normal_variants.ids
grep -v '^#' $OUTPUT_DIR/mutect2/somatic_mutations.vcf | awk '{print $3}' > $OUTPUT_DIR/cancer_variants.ids
comm -23 $OUTPUT_DIR/cancer_variants.ids $OUTPUT_DIR/normal_variants.ids > $OUTPUT_DIR/somatic_mutations.ids

echo "Somatic mutations detected. Check $OUTPUT_DIR/somatic_mutations.ids for details."
echo "NGS Data Analysis pipeline completed."
