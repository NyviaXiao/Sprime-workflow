# Publication-style chromosome-body landscape from final classification rows.
library(karyoploteR)

normalize_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

if (snakemake@params[["genome_build"]] != "GRCh37") {
  stop("Landscape plotting currently supports genome_build: GRCh37 only")
}

classes <- read.delim(snakemake@input[["classification"]], check.names=FALSE,
                      na.strings="NA")
chromosomes <- unique(normalize_chr(snakemake@params[["chromosomes"]]))
classes$chr <- normalize_chr(classes$chromosome)
classes <- classes[classes$chr %in% chromosomes, ]

# Muted, print-friendly colors. Grey denotes displayed chromosome positions not
# covered by a final classified segment; it does not imply no introgression.
background <- "#DDE3E6"
colors <- c(neanderthal="#3F6F9F", denisovan="#A95243", ambiguous="#7D8790")

png(snakemake@output[["introgression"]], width=2600, height=1700, res=250)

main_title <- paste(
  snakemake@wildcards[["population"]],
  "- Introgression Landscape"
)

plot.params <- getDefaultPlotParams(plot.type=6)

# 给标题和底部 legend 留出独立空间
plot.params$topmargin <- 80
plot.params$bottommargin <- 75

kp <- plotKaryotype(
  genome="hg19",
  chromosomes=chromosomes,
  plot.type=6,
  plot.params=plot.params,
  main=main_title,
  cex=1.45,
  font=2
)

lengths <- kp$chromosome.lengths[chromosomes]

# A uniform grey rectangle establishes the chromosome body. Classification
# intervals are then drawn directly into the same vertical body, not an
# external track above or below the chromosome.
kpRect(kp, chr=chromosomes, x0=0, x1=lengths, y0=0, y1=1,
       data.panel="ideogram", r0=0, r1=1, col=background, border=NA)
for (label in names(colors)) {
  rows <- classes[classes$archaic_class == label, ]
  if (nrow(rows)) {
    kpRect(kp, chr=rows$chr, x0=as.numeric(rows$start), x1=as.numeric(rows$end),
           y0=0, y1=1, data.panel="ideogram", r0=0, r1=1,
           col=colors[[label]], border=NA)
  }
}

# Titles and the one shared legend sit outside the chromosome panel so neither
# can overlap genomic intervals or chromosome labels.
# title(main=paste(snakemake@wildcards[["population"]], "- Introgression Landscape"),
#       line=1, cex.main=1.45, font.main=2)
legend("bottom", legend=c("Neanderthal", "Denisovan", "Ambiguous"),
       fill=unname(colors), horiz=TRUE, bty="n", inset=c(0, -0.08),
       xpd=NA, cex=1.05)
dev.off()
