# Reference-specific affinity landscapes drawn as internal chromosome layers.
# Snakemake passes Neanderthal and Denisovan file vectors separately, preserving
# the explicit reference-to-file mapping established in the workflow rule.
library(karyoploteR)

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

if (snakemake@params[["genome_build"]] != "GRCh37") {
  stop("Landscape plotting currently supports genome_build: GRCh37 only")
}

chromosomes <- unique(normalize_chr(snakemake@params[["chromosomes"]]))

add_colorbar <- function(palette) {
  # A true continuous raster colour bar, rather than three discrete swatches.
  old <- par(no.readonly=TRUE)
  on.exit(par(old), add=TRUE)
  par(fig=c(0.36, 0.64, 0.015, 0.09), new=TRUE, mar=c(1.8, 0, 1.4, 0))
  plot.new(); plot.window(xlim=c(0, 1), ylim=c(0, 1))
  rasterImage(as.raster(matrix(palette, nrow=1)), 0, 0.42, 1, 0.72, interpolate=FALSE)
  axis(1, at=c(0, .5, 1), labels=c("0", "0.5", "1"), tick=FALSE, line=-.4, cex.axis=.9)
  mtext("Match rate", side=3, line=.1, cex=.9)
}

draw_group <- function(output, paths, refs, wanted, palette) {
  height <- max(1500, 1250 + 80 * length(refs))
  png(output, width=2600, height=height, res=250)
  kp <- plotKaryotype(genome="hg19", chromosomes=chromosomes, plot.type=1)
  lengths <- kp$chromosome.lengths[chromosomes]
  n <- length(refs)

  # Every reference is a thin horizontal layer *inside* one chromosome body.
  # Grey backgrounds make absent displayed segments explicit without creating
  # external stacked tracks.
  for (i in seq_along(refs)) {
    lower <- (i - 1) / n
    upper <- i / n
    kpRect(kp, chr=chromosomes, x0=0, x1=lengths, y0=.10, y1=.90,
           data.panel=1, r0=lower, r1=upper, col="#DDE3E6", border=NA)
    rows <- read.delim(paths[[i]], check.names=FALSE, na.strings="NA")
    rows <- rows[rows$archaic_class == wanted, ]
    if (nrow(rows)) {
      rows$chr <- normalize_chr(rows$chromosome)
      rows <- rows[rows$chr %in% chromosomes & is.finite(rows$match_rate), ]
      if (nrow(rows)) {
        idx <- pmax(1, pmin(length(palette), 1 + round(as.numeric(rows$match_rate) * (length(palette) - 1))))
        kpRect(kp, chr=rows$chr, x0=as.numeric(rows$start), x1=as.numeric(rows$end),
               y0=.10, y1=.90, data.panel=1, r0=lower, r1=upper,
               col=palette[idx], border=NA)
      }
    }
  }
  title(main=paste(snakemake@wildcards[["population"]], "-", tools::toTitleCase(wanted), "Affinity Landscape"),
        line=1, cex.main=1.4, font.main=2)
  # Show labels only once in a figure-level key, ordered top to bottom exactly
  # as the internal layers are drawn. No reference text is repeated by chr.
  legend("topright", legend=paste(seq_along(refs), refs),
         title="Reference layers - top to bottom", bty="n", cex=.9,
         inset=c(-.02, 0), xpd=NA)
  add_colorbar(palette)
  dev.off()
}

nean_palette <- colorRampPalette(c("#E6F0F7", "#7FA6C8", "#234F7D"))(101)
deni_palette <- colorRampPalette(c("#E5F3F1", "#63A6A2", "#17666B"))(101)
draw_group(snakemake@output[["neanderthal"]], as.list(snakemake@input[["nean_aff"]]),
           as.character(snakemake@params[["neanderthal_refs"]]), "neanderthal", nean_palette)
draw_group(snakemake@output[["denisovan"]], as.list(snakemake@input[["deni_aff"]]),
           as.character(snakemake@params[["denisovan_refs"]]), "denisovan", deni_palette)
