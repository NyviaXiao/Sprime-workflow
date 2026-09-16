# Draw separate Neanderthal and Denisovan affinity ideograms. Input vectors are
# passed separately by Snakemake, so reference labels can never index the wrong
# affinity file.
library(karyoploteR)

chromosomes <- paste0("chr", as.character(snakemake@params[["chromosomes"]]))
palette <- colorRampPalette(c("#ffffcc", "#41b6c4", "#0c2c84"))(101)

draw_group <- function(output, paths, refs, wanted) {
  png(output, width=1800, height=max(900, 220 * length(refs)), res=180)
  kp <- plotKaryotype(genome="hg19", chromosomes=chromosomes, plot.type=1,
                      main=paste(snakemake@wildcards[["population"]], wanted, "affinity landscape"))
  n <- length(refs)
  for (i in seq_along(refs)) {
    rows <- read.delim(paths[[i]], check.names=FALSE, na.strings="NA")
    rows <- rows[rows$archaic_class == wanted, ]
    if (!nrow(rows)) next
    rows$chr <- ifelse(grepl("^chr", rows$chromosome), rows$chromosome, paste0("chr", rows$chromosome))
    rows <- rows[rows$chr %in% chromosomes & is.finite(rows$match_rate), ]
    if (!nrow(rows)) next
    idx <- pmax(1, pmin(101, 1 + round(as.numeric(rows$match_rate) * 100)))
    kpRect(kp, chr=rows$chr, x0=as.numeric(rows$start), x1=as.numeric(rows$end),
           y0=0.05, y1=0.95, data.panel=1, r0=(i-1)/n, r1=i/n,
           col=palette[idx], border=NA)
    kpAddLabels(kp, labels=refs[[i]], srt=0, pos=3, r0=(i-1)/n, r1=i/n,
                data.panel=1, cex=.7)
  }
  legend("bottom", legend=c("0", "0.5", "1"), fill=palette[c(1, 51, 101)],
         title="match rate", horiz=TRUE, bty="n")
  dev.off()
}

draw_group(snakemake@output[["neanderthal"]], as.list(snakemake@input[["nean_aff"]]),
           as.character(snakemake@params[["neanderthal_refs"]]), "neanderthal")
draw_group(snakemake@output[["denisovan"]], as.list(snakemake@input[["deni_aff"]]),
           as.character(snakemake@params[["denisovan_refs"]]), "denisovan")
