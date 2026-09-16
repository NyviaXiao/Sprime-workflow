# Draw one chromosome-ideogram landscape from final classification rows.
library(karyoploteR)

classes <- read.delim(snakemake@input[["classification"]], check.names=FALSE,
                      na.strings="NA")
chromosomes <- paste0("chr", as.character(snakemake@params[["chromosomes"]]))
classes$chr <- ifelse(grepl("^chr", classes$chromosome), classes$chromosome,
                      paste0("chr", classes$chromosome))
classes <- classes[classes$chr %in% chromosomes, ]
colors <- c(neanderthal="#2c7bb6", denisovan="#d7191c", ambiguous="#999999")

png(snakemake@output[[0]], width=1800, height=1000, res=180)
kp <- plotKaryotype(genome="hg19", chromosomes=chromosomes, plot.type=1,
                    main=paste(snakemake@wildcards[["population"]], "introgression landscape"))
for (label in names(colors)) {
  rows <- classes[classes$archaic_class == label, ]
  if (nrow(rows)) {
    kpRect(kp, chr=rows$chr, x0=as.numeric(rows$start), x1=as.numeric(rows$end),
           y0=0.15, y1=0.85, data.panel=1, r0=0, r1=1, col=colors[[label]], border=NA)
  }
}
legend("bottom", legend=names(colors), fill=colors, horiz=TRUE, bty="n")
dev.off()
