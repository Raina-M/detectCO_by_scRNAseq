# Pipeline
### 1. GMRs selection on phased reference genome
Genotyping markers on reference (GMRs) are able to distinguish two haplotypes, and thus they will also be used in gametes for genotyping.
```
# ----- Mapping NGS reads to haplotype 1 of phased reference genome -----
# index reference
bowtie2 build --threads $CPU reference_hap1.fa reference_hap1_index
# alignments
bowtie2 -p $CPU -x reference_hap1_index -1 read1.fq.gz -2 read2.fq.gz | samtools sort -@ $CPU > aln.sorted.bam

# ----- SNP calling -----
bcftools mpileup -Oz --threads $CPU -o aln.mpileup.gz -f reference_hap1.fa aln.sorted.bam

echo -e "*\t*\t*\t*\t2" > ploidy.txt

bcftools call -Avm -Ov --threads $CPU --ploidy-file ploidy.txt aln.mpileup.gz > aln.vcf
```
