# Separate *R. breviuscula* scRNA-seq data from mixed pollen nuclei sequences 
For multiplexing purpose, we prepared a scRNA-seq library for mixed *R. breviuscula* and *R. tenuis* pollen nuclei. Therefore, we need to separate *R. breviuscula* data out before doing any further analysis.

## Step 1 - Preprocessing of raw scRNA sequencing data
To correct cell barcodes
```
bcctools correct --spacer 12 \    # 10x v3 library: 16 barcode + 12 UMI
                  cellranger_whitelist_10Xv3.txt \    # downloaded from 10X Genomics
                  R1.fastq.gz R2.fastq.gz \    # scRNA-seq paired reads
                  > corrected.tsv
```
```
# filter barcode
cut -f2,4 a4984_merged_corrected.tsv > barcode_umi.list
cut -f1 barcode_umi.list | sort | uniq -c > BC.stats  # record the number of occurence of each barcode
awk '{print $2}' BC.stats | awk -F"," 'NF==1' | sed '/^*/d' > uniq_BC.list
awk 'FNR==NR{a[$1]=$1; next}; $2 in a {print $0;}' uniq_BC.list BC.stats > uniq_BC.stats
awk '$1>5000 {print $2}' uniq_BC.stats > uniq_BC_gt5k_reads.list  # 365,771,748 reads

awk 'FNR==NR{a[$1]=$1; next}; $2 in a {print $0;}' uniq_BC_gt5k_reads.list a4984_merged_corrected.tsv > a4984_merged_cells_gt5k_reads.tsv
```

## Step 2 - Split of cells by barcodes / demultiplexing
split one single tsv file from bcctools correct into cell/barcode-specific files
```
cd ${WD}/2_split_cells
awk '{print>$2".tsv"}' ${WD}/1_BC_correction/a4984_merged_cells_gt5k_reads.tsv

# convert .tsv to .fastq format
for f in *.tsv
do
    BC=${f%.tsv*}
    cd $WD
    [ ! -d ${WD}/2_split_cells/$BC ] && mkdir -p ${WD}/2_split_cells/$BC
    mv $f ${WD}/2_split_cells/$BC
    cd ${WD}/2_split_cells/$BC

    awk '{print "@"$1"_"$4"\n"$5"\n+\n"$9}' $f  > ${BC}_R1.fastq
    awk '{print "@"$1"_"$4"\n"$6"\n+\n"$10}' $f > ${BC}_R2.fastq
done
```

## Step 3 - ALIGNMENTS & SPIECIS ASSIGNMENT
build reference genome index files
```
cd $WD
# R. brevi. genome
mkdir rb_index
cd rb_index
ln -s $rbHap1 .
hisat2-build -p16 Rbrevi_nuclear_hap1_chrs.fasta rb_chrs_genome_index
```
R. tenuis genome
extract only chromosomal seqs
```
cd $WD
mkdir rt_index
cd rt_index
ln -s $rtHap1 .
samtools faidx Rtenuis.hap1.all.sequences.fasta `echo -e "Chr1_h1\nChr2_h1"` > Rtenuis_hap1_chrs.fasta
unlink Rtenuis.hap1.all.sequences.fasta

hisat2-build -p16 Rtenuis_hap1_chrs.fasta rt_chrs_genome_index

# mapping to R. breviuscula
cd ${WD}/2_split_cells
for BC in *
do
    cd ${WD}/2_split_cells/$BC
    # bcftools does not call SNPs on reads whose mated reads are not mapped
    # most scRNA-seq read1 cannot be mapped, so merge two reads and map as SE mode
    awk '{if(NR%4==1) print $0"\tR1"; else print $0}' ${BC}_R1.fastq > ${BC}_R1.tmp.fastq
    awk '{if(NR%4==1) print $0"\tR2"; else print $0}' ${BC}_R2.fastq > ${BC}_R2.tmp.fastq
    cat ${BC}_R1.tmp.fastq ${BC}_R2.tmp.fastq > ${BC}.fastq
    rm ${BC}_R1.tmp.fastq
    rm ${BC}_R2.tmp.fastq

    hisat2 -p 8 -x ${WD}/rb_index/rb_chrs_genome_index \
                -U ${BC}.fastq \
                -S ${BC}_Rb.sam \
                --summary-file hisat_map2rb_stdout.log
    samtools sort -@4 ${BC}_Rb.sam > ${BC}_Rb.sorted.bam
    rm ${BC}_Rb.sam
done

# mapping to R. tenuis
cd ${WD}/2_split_cells
for BC in *
do
    cd ${WD}/2_split_cells/$BC
    hisat2 -p 8 -x ${WD}/rt_index/rt_chrs_genome_index \
                -U ${BC}.fastq \
                -S ${BC}_Rt.sam \
                --summary-file hisat_map2rt_stdout.log
    samtools sort -@4 ${BC}_Rt.sam > ${BC}_Rt.sorted.bam
    rm ${BC}_Rt.sam
done

# compare alignment read for each cell
# to determine which species should this cell be assigned
cd ${WD}/2_split_cells
for BC in *
do
    cd ${WD}/2_split_cells/$BC
    rb_rate=`tail -n 1 hisat_map2rb_stdout.log | cut -d"%" -f1`
    rt_rate=`tail -n 1 hisat_map2rt_stdout.log | cut -d"%" -f1`
    readnum=`head -n 1 hisat_map2rb_stdout.log | cut -d" " -f1`
    echo -e $BC"\t"$readnum"\t"$rb_rate"\t"$rt_rate >> ${WD}/stats_and_plots/aln_rates.stats
done
```

## Step 4 - COLLAPSE UMI & MORE FILTERING
```
cd ${WD}/2_split_cells
for BC in *
do
    cd ${WD}/2_split_cells/$BC
    samtools index ${BC}_Rb.sorted.bam
    samtools index ${BC}_Rt.sorted.bam
    
    cd $UMIcollapse
    ./umicollapse bam -i ${WD}/2_split_cells/$BC/${BC}_Rb.sorted.bam  \
                      -o ${WD}/2_split_cells/$BC/${BC}_Rb.sorted.dedup.bam \
                      > ${WD}/2_split_cells/$BC/umi_collapse_rb.stdout
    ./umicollapse bam -i ${WD}/2_split_cells/$BC/${BC}_Rt.sorted.bam  \
                      -o ${WD}/2_split_cells/$BC/${BC}_Rt.sorted.dedup.bam \
                      > ${WD}/2_split_cells/$BC/umi_collapse_rt.stdout
done
# count number of reads after collasping UMIs
cd ${WD}/2_split_cells
for BC in *
do
    cd ${WD}/2_split_cells/$BC
    rb_input_rnum=`cut -f2 umi_collapse_rb.stdout | sed -n '3p'`
    rb_rnum_dedup=`cut -f2 umi_collapse_rb.stdout | sed -n '9p'`
    
    rt_input_rnum=`cut -f2 umi_collapse_rt.stdout | sed -n '3p'`
    rt_rnum_dedup=`cut -f2 umi_collapse_rt.stdout | sed -n '9p'`
    
    echo -e $BC"\t"$rb_input_rnum"\t"$rb_rnum_dedup"\t"$rt_input_rnum"\t"$rt_rnum_dedup >> ${WD}/stats_and_plots/readnum_dedup.stats
done

# see ${WD}/stats_and_plots R script - determine cutoffs to split species
awk '$3/($3+$5)>=0.67 {print $1}' readnum_dedup.stats > Rb_BC_dedup.list
awk '$5/($3+$5)>=0.6 {print $1}' readnum_dedup.stats > Rt_BC_dedup.list
# output the cell barcodes that were assigned to "R. breviuscula"
# and the alingment rate is no less than 25%
grep -wFf Rb_BC_dedup.list aln_rates.stats | awk '$3>=25 {print $1}' > Rb_BC_dedup_alnrate_gt25.list
```
