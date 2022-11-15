# Pipeline
### 1. GMRs selection on phased reference genome
Genotyping markers on reference (GMRs) are able to distinguish two haplotypes, and thus they will also be used for gametes genotyping. You can choose other tools to do alignments and SNP calling. Here I show an example by `bowtie2` and `bcftools`.
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

# ----- Converting VCF -----
SHOREmap convert --marker aln.vcf --folder out_dir -runid myID > convert.log

# ----- Select allelic SNPs as GMRs -----
# mapping quality over 50, alternative allele coverage between 5 and 30, and allele frequency between 0.4 nd 0.6
awk '$6>50 && $7>=5 && $7<=30 && $8>=0.4 && $8<=0.6' myID_converted_variant.txt > GMRs.txt
```
The threshold of selecting GMRs could vary case by case, so it is suggested to plot the distribution of important features and the make a decision by the observation.

### 2. Preprocessing scRNA-seq data
Raw scRNA-seq data usually include barcode errors and contaminations such as doublets, ambient RNA, etc. So we need to correct those errors beforehand.
```
bcctools correct --alts 16 --spacer 12 \            # 10x v3 library: 16 barcode + 12 UMI
                  cellranger_whitelist_10Xv3.txt \  # 10x v3 whitelist, can be downloaded from 10x website
                  R1.fastq.gz R2.fastq.gz \         # 10x raw reads
                  > sc_reads_corrected.tsv          # corrected sequences with extracted barcodes and UMIs

# filter barcode:
    # record the number of occurence of each barcode
    # then remove the cell with multiple barcodes
    # or wihtout determined barcodes
cut -f2 sc_reads_corrected.tsv | sort | uniq -c | awk -F"," 'NF==1' | grep -v "*" > viable_BC.stats    

# select cells with over a certain amount of reads
awk '$1>5000 {print $2}' vaible_BC.stats > viable_BC_gt5k_reads.list
awk 'FNR==NR{a[$1]=$1; next}; $2 in a {print $0;}' viable_BC_gt5k_reads.list sc_reads_corrected.tsv > sc_reads.tsv
```
