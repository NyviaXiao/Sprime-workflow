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
  # Draw the continuous legend in a margin-free inset so the short panel
  # cannot trigger R's "figure margins too large" error.
  old <- par(no.readonly=TRUE)
  on.exit(par(old), add=TRUE)
  par(fig=c(0.32, 0.68, 0.015, 0.13), new=TRUE,
      mar=c(0, 0, 0, 0), xpd=NA)
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
  rasterImage(as.raster(matrix(palette, nrow=1)),
              0.03, 0.40, 0.97, 0.62, interpolate=FALSE)
  rect(0.03, 0.40, 0.97, 0.62, border="#B8C3C8", lwd=0.8)
  text(x=c(0.03, 0.5, 0.97), y=0.24,
       labels=c("0", "0.5", "1"), cex=0.82)
  text(x=0.5, y=0.82, labels="Match rate", cex=0.88)
}

map_reference_labels <- function(ids, configured) {
  labels <- as.character(ids)
  if (!length(configured)) return(labels)
  for (i in seq_along(labels)) {
    hit <- Filter(function(ref) identical(as.character(ref$id), labels[[i]]), configured)
    if (length(hit) && !is.null(hit[[1]]$tag)) labels[[i]] <- as.character(hit[[1]]$tag)
  }
  labels
}

display_labels <- function(ids) {
  # Tags are presentation labels only; IDs still control all file mapping.
  configured <- list()
  if ("config" %in% slotNames(snakemake)) configured <- snakemake@config[["archaic_references"]]
  map_reference_labels(ids, configured)
}

draw_group <- function(output, paths, refs, wanted, palette) {
  height <- max(1500, 1250 + 80 * length(refs))
  png(output, width=2800, height=height, res=250)
  par(mar=c(7, 5, 8, 13), oma=c(0, 0, 0, 0))
  kp <- plotKaryotype(genome="hg19", chromosomes=chromosomes, plot.type=6)
  lengths <- kp$chromosome.lengths[chromosomes]
  n <- length(refs)

  # Every reference is a thin horizontal layer *inside* one chromosome body.
  # Grey backgrounds make absent displayed segments explicit without creating
  # external stacked tracks.
  for (i in seq_along(refs)) {
    # karyoploteR's y-axis grows bottom-to-top. Reverse the index so the first
    # configured reference is physically the top internal layer, matching the
    # single figure-level legend.
    lower <- (n - i) / n
    upper <- (n - i + 1) / n
    kpRect(kp, chr=chromosomes, x0=0, x1=lengths, y0=0, y1=1,
           data.panel="ideogram", r0=lower, r1=upper, col="#DDE3E6", border=NA)
    rows <- read.delim(paths[[i]], check.names=FALSE, na.strings="NA")
    rows <- rows[rows$archaic_class == wanted, ]
    if (nrow(rows)) {
      rows$chr <- normalize_chr(rows$chromosome)
      rows <- rows[rows$chr %in% chromosomes & is.finite(rows$match_rate), ]
      if (nrow(rows)) {
        idx <- pmax(1, pmin(length(palette), 1 + round(as.numeric(rows$match_rate) * (length(palette) - 1))))
        kpRect(kp, chr=rows$chr, x0=as.numeric(rows$start), x1=as.numeric(rows$end),
               y0=0, y1=1, data.panel="ideogram", r0=lower, r1=upper,
               col=palette[idx], border=NA)
      }
    }
  }
  title(main=paste(snakemake@wildcards[["population"]], "-", tools::toTitleCase(wanted), "Affinity Landscape"),
        line=3.5, cex.main=1.4, font.main=2)
  # A single right-side annotation block identifies the internal layers. It is
  # deliberately separate from the chromosome panel and never repeats by chr.
  legend("right", legend=paste(seq_along(refs), display_labels(refs)),
         title="Reference layers", bty="n", cex=.9,
         inset=c(-.14, 0), xpd=NA)
  add_colorbar(palette)
  dev.off()
}

nean_palette <- colorRampPalette(c("#E6F0F7", "#7FA6C8", "#234F7D"))(101)
deni_palette <- colorRampPalette(c("#E5F3F1", "#63A6A2", "#17666B"))(101)
draw_group(snakemake@output[["neanderthal"]], as.list(snakemake@input[["nean_aff"]]),
           as.character(snakemake@params[["neanderthal_refs"]]), "neanderthal", nean_palette)
draw_group(snakemake@output[["denisovan"]], as.list(snakemake@input[["deni_aff"]]),
           as.character(snakemake@params[["denisovan_refs"]]), "denisovan", deni_palette)
