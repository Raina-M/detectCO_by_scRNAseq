# Separate *R. breviuscula* scRNA-seq data from mixed pollen nuclei sequences 
For multiplexing purpose, we prepared a scRNA-seq library for mixed *R. breviuscula* and *R. tenuis* pollen nuclei. Therefore, we need to separate *R. breviuscula* data out before doing any further analysis.

Preprocessing of raw scRNA sequencing data, pipeline

build reference genome index files

R. tenuis genome
extract only chromosomal seqs
```
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
