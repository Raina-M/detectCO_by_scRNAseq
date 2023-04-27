# Separate *R. breviuscula* scRNA-seq data from mixed pollen nuclei sequences 
For multiplexing purpose, we prepared a scRNA-seq library for mixed *R. breviuscula* and *R. tenuis* pollen nuclei. Therefore, we need to separate *R. breviuscula* data out before doing any further analysis.

### Step 1: Preprocessing of raw scRNA sequencing data
ScRNA-seq
```
bcctools correct --spacer 12 \     # 10x v3 library: 16 barcode + 12 UMI
                  cellranger_whitelist_10Xv3.txt \
                  a4984_merged_R1.fastq.gz a4984_merged_R2.fastq.gz \
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
