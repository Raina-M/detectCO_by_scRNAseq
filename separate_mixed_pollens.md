# Separate *R. breviuscula* scRNA-seq data from mixed pollen nuclei sequences 
For multiplexing purpose, we prepared a scRNA-seq library for mixed *R. breviuscula* and *R. tenuis* pollen nuclei. Therefore, we need to separate *R. breviuscula* data out before doing any further analysis. To this end, we firstly aligned the scRNA-seq reads of all demultiplexed cells to both *R. breviuscula* and *R. tenuis* reference genome and then removed PCR-related quantification biases (AKA deduplication) using [UMIcollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse). This part has been delineated in the [main pipeline](https://github.com/Raina-M/detectCO_by_scRNAseq/blob/main/pipeline.md). Now, we resume from obtaining alingment results of each cell to *R. breviuscula* and *R. tenuis* genome. For a single cell, if the fraction of aligned reads to *R. breviuscula* is much higher than that to *R. tenuis*, then we say this cell is from *R. breviuscula* and vice versa.

The standard output of UMIcollapse derives some informative statistics, among which we can extract the number of reads mapping to provided reference genome.
```
# count number of reads after collasping UMIs / deduplication
cd demultiplexed_cells
for BC in *
do
    cd /path/demultiplexed_cells/$BC
    rb_rnum_dedup=`cut -f2 umi_collapse_rb.stdout | sed -n '9p'`  # Number of reads mapped to R. breviuscula after deduplicating
    rt_rnum_dedup=`cut -f2 umi_collapse_rt.stdout | sed -n '9p'`  # Number of reads mapped to R. tenuis      after deduplicating
    
    echo -e $BC"\t"$rb_rnum_dedup"\t"$rt_rnum_dedup >> $outDir/readnum_dedup.stats
done
```
Now we can analyze the alignments and try to find out the origin of each cell. It is suggested to plot the distributions of alignment fraction, i.e., for each cell, alignment fraction = number of reads aligned to R. breviuscula / (number of reads aligned to R. breviuscula + number of reads aligned to R. tenuis). Please refer to (our paper)[^1] for more details. We took 0.67 and 0.6 as cutoff for *R. breviuscula* and *R. tenuis* respectively.
```
awk '$3/($3+$5)>=0.67 {print $1}' readnum_dedup.stats > Rb_BC_dedup.list  # cells from R. breviuscula
awk '$5/($3+$5)>=0.6  {print $1}' readnum_dedup.stats > Rt_BC_dedup.list  # cells from R. tenuis
```
Following this, we can see that not all cells were assigned to a certain species, there are cells with similar alignment reads to both species. They are potentially doublets that captured cells from different species, which should be discarded.

[^1]: Castellani, Marco, et al. "Meiotic recombination dynamics in plants with repeat-based holocentromeres shed light on the primary drivers of crossover patterning." bioRxiv (2023): 2023-04.
