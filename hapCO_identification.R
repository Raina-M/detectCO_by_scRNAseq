# CO identification in gamete nuclei (Note: only homozygous change of alleles would be captured).
# writen by Hequan Sun, MPIPZ Email: sunhequan@gmail.com

# Modified by Meng Zhang (mzhang@mpipz.mpg.de)
# on May 25th, 2022
# Adapted to Rhynchospora CO identification with scRNA-seq from pollens

# Last modification: Sep 28th, 2022, change smooth function 1 by AF sooth, not genotype


library("optparse", lib.loc="/home/mzhang/R/x86_64-pc-linux-gnu-library/4.0/")

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help="Input file with single cell SNP information", metavar="file"),
  make_option(c("-p", "--prefix"), type = "character", default="out", 
              help = "Prefix of output file names [default: %default]",
              metavar = "string"),
  make_option(c("-g", "--genome"), type = "character", default = NULL,
              help = "Chromosme size file with chromosme name and size in each line, tab-delimited",
              metavar = "file"),
  make_option(c("-o", "--outpath"), type="character", default = getwd(),
              help = "Output file path [default: current working directory]",
              metavar = "file_path"),
  make_option(c("-c","--cell_markers"), type = "integer", default = 400,
              help = "Minimum number of markers required in the cell [default: %default]",
              metavar = "int"),
  make_option(c("-s","--block_size"), type = "integer", default = 2000000,
              help = "Minimum block size to define genotype [default: %default]",
              metavar = "int"),
  make_option(c("-n","--marker_num"), type = "integer", default = 5,
              help = "Minimum number of markers required in blocks to determine block genotypes [default: %default]",
              metavar = "int"),
  make_option("--window", type = "integer", default = 5,
              help = "Window size for smoothing [default: %default] (Cannot be modified currently, will be updated later)",
              metavar = "int"),
  make_option("--baseAF", type = "double", default = 0.3,
              help = "Base allele frequency threshold in 1st-round smoothing [default: %default]",
              metavar = "num"),
  make_option("--windowAF", type = "double", default = 0.4, metavar = "num",
              help = "Window allele frequency threshold in 2nd-round smoothing [default: %default]"),
  make_option("--genotype", type = "double", default = 0.2, metavar = "num",
              help = "Threshold to define block genotypes [default: %default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$genome)){
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input and genome file).", call.=FALSE)
}

### ========== functions =======================================================
# function for smoothing the allele counts
allele_count_smoother<-function(acnt_chr)
{
  # note: this only works at chr-wise; otherwise you need to split junctions between chrs.
  this_dim <- dim(acnt_chr)
  if(this_dim[2] != 6)
  {
    stop("Wrong data of acnt. \n")
  }else
  {
    # add two new cols: 1) rawAF=count_of_ref)/count_of_ref+count_of_alt; and 2) smtAF
    acnt_plus_af <- cbind(acnt_chr, acnt_chr[, 4]/(acnt_chr[, 4]+ acnt_chr[, 6]), rep(0, this_dim[1]) )
    colnames(acnt_plus_af) <- c(paste("V", 1:6, sep=""), "rawAF", "smtAF")
  }
  # smoothing: two iteration: with first iteration, there are problems at borders.
  # 1st-round smoothing
  adj_marker_num = 2
  for(row in 1:this_dim[1])
  {
    surrounding_case <- max(row-adj_marker_num, 1):min(row+adj_marker_num, this_dim[1])
    if(row==1)
    {
      rigth_case       <- (row+1):min(row+adj_marker_num, this_dim[1])
      # only sum(presence/absence) <= no read count effect: the site itself no in sum
      right_sum        <- sum(acnt_plus_af[rigth_case, 7]>0)
      window_af        <- right_sum / (length(surrounding_case)-1)        
    }else if(row==this_dim[1])
    {
      left_case        <- max(row-adj_marker_num, 1):(row-1)
      # only sum(presence/absence) <= no read count effect: the site itself no in sum
      left_sum         <- sum(acnt_plus_af[left_case, 7]>0)
      window_af        <- left_sum  / (length(surrounding_case)-1)        
    }else
    {
      left_case        <- max(row-adj_marker_num, 1):(row-1)
      rigth_case       <- (row+1):min(row+adj_marker_num, this_dim[1])
      # only sum(presence/absence) <= no read count effect: the site itself no in sum
      left_sum         <- sum(acnt_plus_af[left_case, 7])
      right_sum        <- sum(acnt_plus_af[rigth_case, 7])
      window_af        <- (left_sum + right_sum) / (length(surrounding_case)-1)       
    }
    #acnt_plus_af[row, 8] <- ifelse(window_af>0.5, 1, 0)
    if(window_af > (1-opt$baseAF))
    {
      acnt_plus_af[row, 8] <- 1 # smoothed
    }else if(window_af < opt$baseAF)
    {
      acnt_plus_af[row, 8] <- 0 # smoothed
    }else
    {
      acnt_plus_af[row, 8] <- acnt_plus_af[row, 7] # no change
    }
  }
  # 2nd round smoothing
  for(row in 1:this_dim[1])
  {
    surrounding_case     <- max(row-adj_marker_num, 1):min(row+adj_marker_num, this_dim[1])
    if(row==1)
    {
      rigth_case       <- (row+1):min(row+adj_marker_num, this_dim[1])
      # only sum(presence/absence) <= no read count effect: the site itself no in sum
      right_sum        <- sum(acnt_plus_af[rigth_case, 8]>0)
      window_af        <- right_sum / (length(surrounding_case)-1)        
    }else if(row==this_dim[1])
    {
      left_case        <- max(row-adj_marker_num, 1):(row-1)
      # only sum(presence/absence) <= no read count effect: the site itself no in sum
      left_sum         <- sum(acnt_plus_af[left_case, 8]>0)
      window_af        <- left_sum  / (length(surrounding_case)-1)        
    }else
    {
      left_case        <- max(row-adj_marker_num, 1):(row-1)
      rigth_case       <- (row+1):min(row+adj_marker_num, this_dim[1])
      # only sum(presence/absence) <= no read count effect: the site itself no in sum
      left_sum         <- sum(acnt_plus_af[left_case, 8]>0)
      right_sum        <- sum(acnt_plus_af[rigth_case, 8]>0)
      window_af        <- (left_sum + right_sum) / (length(surrounding_case)-1)       
    }    
    # only sum(presence/absence) <= no read count effect
    # window_af            <- sum(acnt_plus_af[surrounding_case, 8]>0) / length(surrounding_case) 
    #acnt_plus_af[row, 8] <- ifelse(window_af>0.5, 1, 0)
    if(opt$windowAF<=window_af & window_af<=(1-opt$windowAF)) # check here!  2022-09-06
    {
      acnt_plus_af[row, 8] <- acnt_plus_af[row, 7] # reset as initial
    }
  }  
  #
  return(acnt_plus_af)
}

# make genotypes into blocks
get_genotype_block<-function(acnt_valid_chr_afsmoothed)
{
  this_dim      = dim(acnt_valid_chr_afsmoothed)
  # initiation
  this_block    = matrix(data=NA,nrow=1,ncol=4)
  this_block[1] = acnt_valid_chr_afsmoothed[1, 2] # start
  this_block[2] = acnt_valid_chr_afsmoothed[1, 2] # end
  this_block[3] = -1                              # genotype
  this_block[4] = 1                               # marker number
  #
  blocks        = matrix(data=NA,nrow=0,ncol=4)
  #
  if(acnt_valid_chr_afsmoothed[1, 8]>=(1-opt$genotype))
  {
    this_block[3] = 1 # ref
  }else if(acnt_valid_chr_afsmoothed[1, 8]<=opt$genotype)
  {
    this_block[3] = 0 # alt
  }else
  {
    this_block[3] = -1 # het
  }
  
  if(this_dim[1] > 1)
  {
    for (row in 2:this_dim[1])
    {
      if(acnt_valid_chr_afsmoothed[row, 8] != acnt_valid_chr_afsmoothed[row-1, 8])
      {
        blocks <- rbind(blocks, this_block)
        # update for next block
        this_block[1] = acnt_valid_chr_afsmoothed[row, 2] # start
        this_block[2] = acnt_valid_chr_afsmoothed[row, 2] # end
        this_block[3] = -1                                # genotype  
        this_block[4] = 1                                 # marker number
        #
        if(acnt_valid_chr_afsmoothed[row, 8]>=(1-opt$genotype))
        {
          this_block[3] = 1 # ref
        }else if(acnt_valid_chr_afsmoothed[row, 8]<=opt$genotype)
        {
          this_block[3] = 0 # alt
        }else
        {
          this_block[3] = -1 # het
        }      
      }else
      {
        # update current block end
        this_block[2] = acnt_valid_chr_afsmoothed[row, 2]
        this_block[4] = this_block[4] + 1              # marker number
        
      }
    }
  }
  # last block
  blocks <- rbind(blocks, this_block)
  #
  colnames(blocks) <- c("sta", "end", "genotype", "markernum")
  return(blocks)
  #
}

# filter blocks according to block size and marker number and
# get good blocks (supported with markers sufficiently)
filter_blocks<-function(blocks, min_block_size, min_marker)
{
  # min_block_size: length of blocks to keep a block
  # min_marker    : minimum number of markers to keep a block
  new_blocks        = matrix(data=NA,nrow=0,ncol=4)
  colnames(new_blocks) <- c("sta", "end", "genotype", "markernum")
  #  
  this_dim =  dim(blocks)
  for(row in 1:this_dim[1])
  {
    if(blocks[row, 2] - blocks[row, 1] >= min_block_size & blocks[row, 4]>=min_marker)
    {
      new_blocks <- rbind(new_blocks, blocks[row, ])
    }
  }
  #
  return(new_blocks)
}

# find break points between good blocks
get_breakpoints<-function(blocks, chr)
{
  # work on merged blocks
  breakpoints           <- matrix(data=NA,nrow=0,ncol=6)
  #colnames(breakpoints) <- c("chr", "sta", "end", "genotype", "markernum")
  #
  this_dim = dim(blocks)
  tmp_breakpoint        <- matrix(data=NA,nrow=0,ncol=6)
  if(this_dim[1] > 1)
  {
    marker_accumulated    <- blocks[1, 4]
    for(row in 2:this_dim[1])
    {
      if(blocks[row, 3] != blocks[row-1, 3])    # genotype change
      {
        tmp_breakpoint[1]  = chr
        tmp_breakpoint[2]  = blocks[row-1, 2]   # bp start = end   of pre block
        tmp_breakpoint[3]  = blocks[row,   1]   # bp end   = start of next block
        tmp_breakpoint[4]  = blocks[row-1, 3]   # pre  genotype
        tmp_breakpoint[5]  = blocks[row,   3]   # next genotype
        tmp_breakpoint[6]  = marker_accumulated # marker
        breakpoints        = rbind(breakpoints, tmp_breakpoint)
        #
        marker_accumulated = blocks[row, 4]
      }else
      {
        marker_accumulated = marker_accumulated + blocks[row, 4]
      }
    }
  }else
  {
    # no breakpoints
  }
  ##
  return(breakpoints)
}

# finely tuning break points within predicted break point interval
fine_breakpoints<-function(breakpoints, acnt_valid_chr_afsmoothed)
{
  breakpoints_updated <- breakpoints
  this_dim            <- dim(breakpoints)
  for (bp in 1:this_dim[1])
  {
    interval_cnts <- acnt_valid_chr_afsmoothed[acnt_valid_chr_afsmoothed$V2>=breakpoints[bp, 2] & acnt_valid_chr_afsmoothed$V2<=breakpoints[bp, 3], ]
    
    seq           <- interval_cnts$smtAF
    seq_len       <- length(seq)
    max_score     <- 0
    max_pos       <- 1
    for (pos in 1:(length(seq)-1))
    {
      left_cnt_0  <- sum(seq[1:pos]==0)
      left_cnt_1  <- sum(seq[1:pos]==1)
      right_cnt_0 <- sum(seq[(pos+1):seq_len]==0)
      right_cnt_1 <- sum(seq[(pos+1):seq_len]==1)   
      
      af0Lef = left_cnt_0/pos
      af1Lef = left_cnt_1/pos
      
      af0Rig = right_cnt_0/(seq_len-pos)
      af1Rig = right_cnt_1/(seq_len-pos)
      
      score1 = af0Lef * af1Rig
      score2 = af1Lef * af0Rig
      
      scotmp = ifelse(score1>score2, score1, score2)

      if(scotmp > max_score)
      {
        max_score <- scotmp
        max_pos   <- pos
      }
    }
    #
    breakpoints_updated[bp, 2] <- interval_cnts[max_pos,   2]
    breakpoints_updated[bp, 3] <- interval_cnts[max_pos+1, 2]
  }
  return(breakpoints_updated)
}

# 
fine_blocks<-function(breakpoints_updated, final_blocks)
{
  this_dim <- dim(final_blocks)
  for(bl in c(2:this_dim[1]))
  {
    final_blocks[bl-1, 2] <- breakpoints_updated[bl-1, 2]
    final_blocks[bl, 1]   <- breakpoints_updated[bl-1, 3]
  }
  return(final_blocks)
}

# merge smaller blocks into larger ones as final blocks to visualize
make_final_blocks<-function(filtered_blocks, this_chr_size)
{
  # 1       --> bp1_start
  # bp1_end --> bp2_start
  # ...
  #         --> chr_size
  final_blocks <- matrix(data=NA,nrow=0,ncol=4)
  this_block   <- matrix(data=NA,nrow=0,ncol=4)
  this_block[1]<-1 # current block start
  this_block[2]<-1 # current block end
  #
  this_dim = dim(filtered_blocks)
  accumlated_markers <- 0
  row=1
  if(this_dim[1]>1)
  {
    for(row in 2:this_dim[1])
    {
      if(filtered_blocks[row, 3] != filtered_blocks[row-1, 3]) # check geno
      {
        this_block[2] <- filtered_blocks[row-1, 2] # current block end
        this_block[3] <- filtered_blocks[row-1, 3] # current genotype
        this_block[4] <- accumlated_markers        # current markers
        final_blocks  <- rbind(final_blocks, this_block)
        #
        this_block[1]<- filtered_blocks[row, 1] # next block start
        this_block[2]<  filtered_blocks[row, 2] # next block end
      }else
      {
        accumlated_markers = accumlated_markers + filtered_blocks[row, 4]
      }
    }
  }
  if(this_dim[1] > 0)
  {
    this_block[2] <- filtered_blocks[row, 2] # current block end
    this_block[3] <- filtered_blocks[row, 3] # current genotype
    this_block[4] <- accumlated_markers        # current markers
    final_blocks  <- rbind(final_blocks, this_block)
    if(final_blocks[length(final_blocks[, 1]), 2] < this_chr_size)
    {
      final_blocks[length(final_blocks[, 1]), 2] = this_chr_size
    }
  }
  return(final_blocks)
}
### ========== end of function definitions ===================================

min_marker_num <- opt$cell_markers
# wsize  <- 20000 # windows not used at the moment.

# chrsizes
chrsizefile <- opt$genome
df_chrsizes <- read.table(chrsizefile, header = FALSE, sep = "\t")
chr_num <- dim(df_chrsizes)[1]
chrsize <- df_chrsizes$V2 

min_chr_size_index <- which.min(chrsize)
max_chr_size_index <- which.max(chrsize)

# this influences double CO detection and noise filtering
min_block_size <- opt$block_size  
# increase this to avoid double crossovers
# however this also missing good cos at chr ends (to improve further)
min_marker    <- opt$marker_num

in_path     <- opt$input
out_path <- opt$outpath

barcode <- opt$prefix
cat("\nInfo: barcode ", barcode)

# prepare output: update the path if necessary
pdf(paste(out_path, "/", barcode, "_allele_cnts_at_markers_sorted_co.pdf", sep=""),
    family="Helvetica", height=9, width=8.26772)
par(mai = c(0.4, 1, 0.1, 0.7)); # margin: bottom, left, top, right
m <- cbind(1:5)
layout(m, heights =rep(1,5))

#
acnt      <- read.table(opt$input)
informativepos <- acnt$V4>0 | acnt$V6>0
if(sum(informativepos)>=min_marker_num){
  acnt_valid     <- acnt[informativepos, ]
  #
  total_smoothed <- 0
  for (chr in 1:chr_num){
    ##
    acnt_valid_chr            <- acnt_valid[acnt_valid$V1==chr, ]
    acnt_valid_chr_afsmoothed <- allele_count_smoother(acnt_valid_chr)
    smooth_applied_sites      <- acnt_valid_chr_afsmoothed[acnt_valid_chr_afsmoothed[, 7] != acnt_valid_chr_afsmoothed[, 8], ]
    total_smoothed            <- total_smoothed + length(smooth_applied_sites$V1)
    ##
    blocks        <- get_genotype_block(acnt_valid_chr_afsmoothed)
    ##
    filtered_blocks <- filter_blocks(blocks, min_block_size, min_marker)
    ##
    breakpoints   <- get_breakpoints(filtered_blocks, chr)
    ## finely tunning breakpoints within initial breakpoint interval
    if(length(breakpoints[, 1])>0){
      breakpoints_updated <-fine_breakpoints(breakpoints, acnt_valid_chr_afsmoothed)
    }
    else{
      breakpoints_updated <- breakpoints
    }
    ##
    final_blocks2  <- make_final_blocks(filtered_blocks, chrsize[chr])

    ##
    if(length(breakpoints[, 1])>0){
      final_blocks   <-fine_blocks(breakpoints_updated, final_blocks2)
    }
    else{
      final_blocks <- final_blocks2
    }
    ## visualize smoothed result
    ## 
    this_ylim = 30
    plot(acnt_valid_chr[acnt_valid_chr$V4>0, 2], 
         acnt_valid_chr[acnt_valid_chr$V4>0, 4], 
         col="red", type="h", ylim=c(-1*this_ylim,this_ylim), xlim=c(1, max(chrsize)), axes=F, xlab="Position", ylab="Allele count")
    axis(1, at=c(seq(0, chrsize[chr], 1e+07), chrsize[chr]), labels = round(c(seq(0, chrsize[chr], 1e+07), chrsize[chr])/1e+06), cex.axis=1.0)
    axis(2, at = c(-1*this_ylim, 0, this_ylim), labels = c(this_ylim, 0, this_ylim))
    points(acnt_valid_chr[acnt_valid_chr$V6>0, 2], 
           acnt_valid_chr[acnt_valid_chr$V6>0, 6]*-1, 
           col="blue", type="h")
    lines(c(0,chrsize[chr]),c(0,0),col="black",lty=2)
    #
    mtext(paste("(Mb chr", chr, ")\n"), at= chrsize[chr]+5e+06,  side=1, line=2.0, cex=0.7)
    #
    rect(1, this_ylim-3, chrsize[chr], this_ylim+5, col="gray", border = NA) # background
    #
    fb_dim = dim(final_blocks)
    if(fb_dim[1]>0){
      rect(as.vector(final_blocks[, 1]),
           this_ylim-3, 
           as.vector(final_blocks[, 2]),
           this_ylim+5, 
           col=as.vector(ifelse(final_blocks[, 3]==1, "red", ifelse(final_blocks[, 3]==0, "blue", "purple"))), 
           border = NA)
    }
    ##
    if(chr==min_chr_size_index){
      legend(chrsize[chr]+round(0.5*(chrsize[max_chr_size_index]-chrsize[chr])), 
             10,
             pch    = c(15, 15),
             col    = c("red", "blue"),
             legend = c(paste("genotype A",   sep=""), 
                        paste("genotype B", sep="")), 
             horiz  = F,
             border = "gray",
             #bty    = "n",
             #text.width=8,
             cex=1)
    }
    ## collect output
    if(chr==1){
      write.table(breakpoints_updated, file = paste(out_path, "/", barcode, "_allele_cnts_at_markers_sorted_co_pred.txt", sep=""), quote=F, row.names = F, col.names = F, append=F)
      final_blocks_chr <- cbind(rep(chr, length(final_blocks[,1])), final_blocks)
      write.table(final_blocks_chr, file = paste(out_path, "/", barcode, "_allele_cnts_at_markers_sorted_co_block_pred.txt", sep=""), quote=F, row.names = F, col.names = F, append=F)          
    }
    else{
      write.table(breakpoints_updated, file = paste(out_path, "/", barcode, "_allele_cnts_at_markers_sorted_co_pred.txt", sep=""), quote=F, row.names = F, col.names = F, append=T)     
      final_blocks_chr <- cbind(rep(chr, length(final_blocks[,1])), final_blocks)
      write.table(final_blocks_chr, file = paste(out_path, "/", barcode, "_allele_cnts_at_markers_sorted_co_block_pred.txt", sep=""), quote=F, row.names = F, col.names = F, append=T)          
    }
    ##
  } # end of chr
  #
  cat("\nInfo: total smoothed sites: ", total_smoothed, 
      " out of ",                          length(acnt_valid$V1), 
      " markers for barcode ",
      barcode,
      "\n",
      sep="")
}else{
  par(mai = c(1, 1, 0.1, 0.7)); # margin: bottom, left, top, right
  plot(1, 1,
       col="red", 
       type="h", 
       axes=F, 
       xlab="Sorry! Markers too limited and no detection of CO could be performed! Check your makers!", 
       ylab="")      
  cat("\nInfo: warning! for barcode ",barcode,      
      " -- too limited number of marekers:", sum(informativepos), "; cannot continue CO identification. \n", sep="")
}
# close
dev.off()
