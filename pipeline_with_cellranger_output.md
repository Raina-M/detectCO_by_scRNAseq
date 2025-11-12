# Crossover calling pipeline using the output from 10X cellranger

This is an alternative pipeline to the `pipeline.md` when you use the output from 10X cellranger to call crossovers.

### 1. Selection of markers on the phased reference genome
Defining the reference markers has no difference from the previous pipeline.

### 2. Preprocessing scRNA-seq data
This step has the main change if you would like to use the alignment results from 10X cellranger. I would show one example run of cellranger but you need to adapt the commands based on your own data, for example, the 10X kit used for preparing the scRNA library.
```
# run cellranger
cellranger count --id $ID \                               # ID name
                 --transcriptome $reference_dir/genome \  # reference genome directory where the outputs from 'cellranger mkref' located
                 --fastqs $seqDir \                       # directory of your raw sequencing data, usually containing index, R1, R2, R3
                 --sample $library_name \                 # normally the library name at the start of your sequence file
                 --include-introns true \
                 --expect-cells 11600 \                   # can be a rough estimation from cell counter
                 --chemistry SC5P-R2 \                    # change this based on the kit you used to prepare for scRNA-seq library
                 --localcores=32 \
                 --localmem=64
```
After the cellranger run, the rest of the steps are the same as the previous pipeline. We need to extract the cells/barcodes that have enough number of reads. The barcode information is included in the output of cellranger alignments. The alignment file is located at `./$ID/outs/possorted_genome_bam.bam`

```
# record the number of occurrences of each barcode
BMA="./$ID/outs/possorted_genome_bam.bam"
samtools view  $BAM | awk '{for(i=1;i<=NF;i++) if($i ~ "^CB:Z:") printf $i"\n"}' | cut -d ':' -f 3 | sort --parallel 16 | uniq -c > read_num.stats
```
We need to determine a threshold to distinguish the real cells from the background. This is based on the occurrence of barcodes. One occurrence of a barcode means one read is contained in this cell. A real cell should have enough reads. Broken or ambient cells normally have quite low read counts. It is suggested to determine this threshold by plotting a "barcode rank plot" as in the cellranger summary report.

I also provided a customised R script to display this plot (`barcode_rank_plot.R`). The only required input is the `read_num.stats` file. The default cutoff between real cells and background is 10k. Change this value in the script based on your case. Feel free to play with different cutoffs to find the optimal "cliff" boundary. 10X cellranger has a detailed interpretation about the [barcode rank plot](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-barcode-rank-plot).

Once the threshold is determined, you can split the `bam` file for each barcode, a step we usually call "demultiplexing". Here shows one example of using cells with at least 10k reads:
```
awk '$1>=10000 {print $2}' ./read_num.stats | while read barcode
do  
  [ ! -d ${WD}/demultiplex/$barcode ] && mkdir -p "${WD}/demultiplex/$barcode"
  cd ${WD}/demultiplex/$barcode

  # extract read names with this barcode
  samtools view -h $BAM | egrep "^@|CB:Z:${barcode}" | samtools sort -o ${barcode%-1}.sorted.bam
  
  # total read number
  read_num=`awk -v bc=$barcode '$2==bc {print $1}' $WD/read_num.stats`

  # unique mapping read number
  uniq_map=`samtools view ${barcode%-1}.sorted.bam | awk '$5>3' | wc -l`
  
  echo $barcode" "$read_num" "$uniq_map >> $WD/aligned_num.stats
done
```
This step is usually time-consuming, so you can write the commands in the `while` loop as a shell script and run it in parallel. After finishing this step, you can do more filtering based on the aligned or uniquely aligned read number across all cells. This is optional and flexible. You can also change the information you would like to note down for the subsequent analysis. Here I use mapping quality more than 3 (`awk '$5>3'`) to define the uniquely mapped reads based on the scoring scheme of mapping quality in `STAR`.

The remaining analysis should be the same as the `pipeline.md`. If you have any questions or suggestions, welcome to discuss with me.
