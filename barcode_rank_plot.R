#!bin/Rscript

read_num <- read.table("read_num.stats")[,1]
cutoff=10000

##### barcode rank plot #####
valid_reads <- sort(read_num[which(read_num>=cutoff)], decreasing = T)
bkgd_reads <- sort(read_num[which(read_num<cutoff)], decreasing = T)

# prepare output: update the path if necessary
pdf("barcode_rank_plot.pdf",
    family="Helvetica", height=6, width=6)
par(mai = c(1.2, 1, 1, 0.5)); # margin: bottom, left, top, right

# valid barcodes
plot(x=log10(seq(1, length(valid_reads), 1)),
     y=log10(sort(valid_reads, decreasing = T)),
     xlim = c(0,6), ylim=c(0,6),
     cex=0.5, col="blue",
     main = "Barcode Rank Plot",
     xlab = "Barcodes", ylab = 'Read counts',
     border = 'black', axes = F,
     cex.lab=1.5, cex.main=1.5)
# background barcodes
points(x=log10(seq(length(valid_reads)+1, length(read_num), 1)),
       y=log10(bkgd_reads),
       pch=1, cex=0.5, col="darkgrey")
box(lwd=2)
axis(side = 1, col.ticks = 1, cex.axis=1.25,
     at=seq(0,6,1), labels = c(1, 10, 100, 1000, "10k", "100k", "1M"))
axis(side = 2, col.ticks = 1, cex.axis=1.25,
     at=seq(0,6,2), labels = c(1, 100, "10k", "1M"))
legend(4, 5.8, legend=c("Cell", "Background"),
       col=c("blue", "darkgrey"), cex=1,
       box.lty=0, pch=1)
dev.off()
