# Crossover calling pipeline

N.B.: Please check the original paper for more deatails.
> The manuscript is still under preparation

### 1. GMRs selection on phased reference genome
Genotyping markers on reference (GMRs) are able to distinguish two haplotypes, and thus they will also be used for gametes genotyping. You can choose other tools to do alignments and SNP calling. Here I show an example by using `bowtie2` and `bcftools`. The derived SNPs were converted by [SHOREmap](http://bioinfo.mpipz.mpg.de/shoremap/).
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
The threshold of selecting GMRs could vary case by case, so it is suggested to plot the distribution of important features and then make a decision based on the observation. It might be worth mentioning that reference genome used here only contains chromosomal scaffolds becasue we can only detect crossover on chromosomes.

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

# select cells with over a certain amount of reads (e.g. 5k reads)
awk '$1>5000 {print $2}' vaible_BC.stats > viable_BC_gt5k_reads.list
awk 'FNR==NR{a[$1]=$1; next}; $2 in a {print $0;}' viable_BC_gt5k_reads.list sc_reads_corrected.tsv > sc_reads.tsv
```
After correction of barcodes and obtaining viable barcodes, we then split the one single raw sequencing file into cell-specific RNA seqeunces based on the barcodes, i.e., all sequences with the same cell barcodes belong to the same cell. This step can also be called demultiplxing.
```
# split one single tsv file from 'bcctools correct' into cell/barcode-specific files
# every file will be named by the 2nd field with .tsv as suffix, i.e. corrected barcode + .tsv
[ ! -d demultiplexed_cells ] && mkdir -p demultiplexed_cells
cd demultiplexed_cells
awk '{print>$2".tsv"}' sc_reads.tsv

# convert .tsv to .fastq format
for file in *.tsv
do
    BC=${file%.tsv*}
    [ ! -d demultiplexed_cells/$BC ] && mkdir -p demultiplexed_cells/$BC
    cd demultiplexed_cells/$BC
    mv ../$file .
    
    # awk '{print "@"$1"_"$4"\n"$5"\n+\n"$9}' $file  > ${BC}_R1.fastq
    # awk '{print "@"$1"_"$4"\n"$6"\n+\n"$10}' $file > ${BC}_R2.fastq
    awk '{print "@"$1"_"$4"\n"$5$6"\n+\n"$9$10}' $file  > ${BC}.fastq
    # rm $file
done
```
In the last step, we need to include UMI in the read header as this information is necessary for deduplication later. As for converting `.tsv` into `.fastq` file, you can also convert into read pairs. Here we converted to a single-end like read because `bcftools` does not call SNPs on reads whose mated reads are not mapped but most scRNA-seq read 1 cannot be mapped. So if you used a different tool for SNP calling in gametes, please check out this point as well.

### 3. Mapping scRNA-seq reads from gametes to reference
If annotations of the genome is not available, hisat2 can be used to map the RNA reads of each cell. If you have annotations, you can of course try STAR or any other RNA aligners.
```
# build reference genome index files
hisat2-build -p $CPU reference_hap1.fa reference_hap1_index

# map scRNA reads to reference for each cell
cd demultiplexed_cells
for BC in *
do
    cd demultiplexed_cells/$BC
    hisat2 -p $CPU -x reference_hap1_index \
                   -U ${BC}.fastq \
                   -S ${BC}.sam \
                   --summary-file aln.stdout
    samtools sort -@ $CPU ${BC}.sam > ${BC}.sorted.bam
    rm ${BC}.sam
done
```
Next, eliminate PCR-related quantification biases using [UMIcollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse), which makes use of the UMIs that we extracted during barcode correction.
```
cd demultiplexed_cells
for BC in *
do
    cd demultiplexed_cells/$BC
    samtools index ${BC}.sorted.bam
    
    umicollapse bam -i ${BC}.sorted.bam  \
                    -o ${BC}.sorted.dedup.bam \
                    > umi_collapse.stdout
done
```
#### - Optional: Separation of cells by species origin 
Before moving forward, it is noteworthy that our single-cell library was prepared for mixed pollens from two different species - *R. breviuscula* and *R. tenuis*, so we need to assign each cell with a species identity. To this end, we mapped reads across all cells to both species, and then compare the alignment rates between two species to determine which species that a certain cell exactly comes from. The alignment rate can be read from `hisat2` log file `aln.stdout`, and read number kept can be read from `umi_collapse.stdout`. If a cell was assigned to a ceratin species but the alignment rate was below 25%, this cell would be discarded.

### 4. SNP calling and selection of markers (GMGs) across gametes
SNP calling for gametes were also done by `bcftools` but you can of course choose other tools. The following code is just for one cell, identified by `$BC` (barcode). In practice, you need loop all valid cells. In the last step, `get_subset.pl` was originally from [TIGER](https://github.com/Imoteph/TIGER_Whole-Genome_Genotyping-by-Sequencing), a tool for genotyping and CO detection for F2 offspring. You can check the original TIGER paper for details, but this script is also included in this github page.
```
# ----- SNP calling -----
# mpileup to generating genotype likelihood
bcftools mpileup -Oz -o ${BC}.mpileup.gz -f reference_hap1.fa ${BC}.sorted.dedup.bam

echo -e "*\t*\t*\t*\t2" > ploidy.txt

# call variants
bcftools call -Am -Ov --ploidy-file ploidy.txt ${BC}.mpileup.gz > ${BC}.vcf

# ----- Converting variants to desired format -----
# extract DP4 lines
bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ${BC}.vcf -o ${BC}.comma.txt

# replace comma separators with tabs
tr ',' '\t' < ${BC}.comma.txt > ${BC}.tabbed.txt
rm ${BC}.comma.txt

# get rid of mitochondria and chloroplast reads;
# create columns for read counts for ref allele and alt allele 
# by adding the first two of the DP4 fields and the second two of the DP4 fields.
awk '{if ($1 !="ChrC" && $1!="ChrM") print $1 "\t" $2 "\t" $3 "\t" $8+$9 "\t"  $4 "\t"  $10+$11}' \
    ${BC}.tabbed.txt > ${BC}.input
rm ${BC}.tabbed.txt

# compare with GMRs and find the overlap between gamete SNPs and GMRs
perl $TIGER/get_subset.pl ${BC}.input 1,2 GMRs.txt 2,3 0 > ${BC}_input_corrected.txt
rm ${BC}.input
```
You are supposed to have genotying markers for each gametes (GMGs) after finishing this step. which would then be used for CO calling. However, not all gametes are viable for CO calling due to contaminations or insufficient markers. Thus, some filtering is needed before identification of COs.

### 5. Gamete filtering by GMGs
Cells with few GMGs are not reliable for CO detection and can be discarded. You need to set a threshold by considering your genome size, GMG numbers across all gametes as well as what CO resolution you expect to get in the end.
To remove doublets, count the times of switches of GMGs??? genotype across gametes. Cells with frequent switches, i.e., switching rate (genotype switching times/number of GMGs) greater than a cutoff (we used 0.07, but you need to find your own cutoff), are doublets.
```
while read BC
do
    awk '$4!=0 || $6!=0' ${BC}_input_corrected.txt > ${BC}_input.tmp
    marker_num=`wc -l ${BC}_input.tmp | cut -d" " -f1`
    # count genotype switching times
    awk '{
        if      ( $4 / ($4+$6) < 0.2) print 0;   # ref genotype
        else if ( $4 / ($4+$6) > 0.8) print 1;   # alt genotype
        else print 0.5;                          # undecided
    }' ${BC}_input.tmp > ${BC}_smt_genotypes.txt # smoothed genotype of this cell
    head -n-1 ${BC}_smt_genotypes.txt > tmp1
    tail -n+2 ${BC}_smt_genotypes.txt > tmp2
    switch_num=`paste tmp1 tmp2 | awk '$2-$1!=0' | wc -l`   # genotype switching times
    echo -e $BC"\t"$marker_num"\t"$switch_num >> switches.stats 

    rm ${BC}_input.tmp
    rm ${BC}_smt_genotypes.txt
    rm tmp*
done < barcode.list
```
```
# cells whose switch rate > 0.07 are doublets, others (<=0.07) are singleton
awk '$2>=400' switches.stats | awk '$3/$2<=0.07 {print $1}' > barcode_gt400markers_no_doublets.list
# also print the singletons whose switch times >= 100
# ! manual check if those cells are doublets after calling CO
# awk '$2>=400' switches.stats | awk '$3/$2<=0.07 && $3>=100 {print $1}' > singletons_gt100switch.list
```

### 6. Crossover calling
Crossovers are identified based on the genotype conversion events. Due to the noisz signlas in scRNA-seq data, smoothing is necessary before define genotype of a certain region to avoid artifact genotype conversion.
```
while read BC
do
    Rscript hapCO_identification.R -i ${BC}_input_corrected.txt \
                                   -p $BC -g reference_hap1.genome -o outdir
done < barcode_gt400markers_no_doublets.list
```
