# Render the three requested genome-wide views from existing tables.
# This intentionally uses base R so the workflow does not require a plotting
# package on the compute node. Coordinates remain the original segment spans.
read_tab <- function(path) read.delim(path, check.names=FALSE, na.strings="NA")
class_file <- snakemake@input[["classification"]]
class_df <- read_tab(class_file)
build_lengths <- c(`1`=249250621,`2`=243199373,`3`=198022430,`4`=191154276,
 `5`=180915260,`6`=171115067,`7`=159138663,`8`=146364022,`9`=141213431,
 `10`=135534747,`11`=135006516,`12`=133851895,`13`=115169878,`14`=107349540,
 `15`=102531392,`16`=90354753,`17`=81195210,`18`=78077248,`19`=59128983,
 `20`=63025520,`21`=48129895,`22`=51304566)
chroms <- as.character(snakemake@params[["chromosomes"]])
if (length(chroms) == 0) chroms <- sort(unique(as.character(class_df$chromosome)))
lengths <- build_lengths[chroms]
missing <- is.na(lengths)
if (any(missing)) {
  for (i in which(missing)) lengths[i] <- max(as.numeric(class_df$end[class_df$chromosome == chroms[i]]), 1)
}
offset <- cumsum(c(0, head(as.numeric(lengths), -1)))
to_x <- function(chr, pos) offset[match(as.character(chr), chroms)] + as.numeric(pos)
cols <- c(neanderthal="#2c7bb6", denisovan="#d7191c", ambiguous="#999999")
draw_intro <- function(path) {
  png(path, width=1800, height=700, res=150); par(mar=c(5,4,3,1))
  plot(NA, xlim=c(0, sum(lengths)), ylim=c(0,1), xaxt="n", yaxt="n", xlab="Genomic position (chr1–22)", ylab="Segments", main=paste(snakemake@wildcards[["population"]], "introgression landscape"))
  for (i in seq_along(chroms)) rect(offset[i], .1, offset[i]+lengths[i], .9, col="#f4f4f4", border="white")
  for (cl in names(cols)) { z <- class_df[class_df$archaic_class == cl,]; if (nrow(z)) rect(to_x(z$chromosome,z$start), .25, to_x(z$chromosome,z$end), .75, col=cols[cl], border=NA) }
  axis(1, at=offset+lengths/2, labels=paste0("chr", chroms)); legend("top", legend=names(cols), fill=cols, horiz=TRUE, bty="n"); dev.off()
}
draw_affinity <- function(path, refs, wanted) {
  png(path, width=1800, height=max(500, 180*length(refs)), res=150); par(mar=c(5,8,3,1))
  plot(NA, xlim=c(0,sum(lengths)), ylim=c(.5,length(refs)+.5), xaxt="n", yaxt="n", xlab="Genomic position (chr1–22)", ylab="Archaic reference", main=paste(snakemake@wildcards[["population"]], wanted, "affinity landscape"))
  for (i in seq_along(chroms)) rect(offset[i], .5, offset[i]+lengths[i], length(refs)+.5, col=ifelse(i %% 2, "#fafafa", "#f0f0f0"), border="white")
  aff_paths <- as.character(snakemake@input[["aff"]])
  for (i in seq_along(refs)) { d <- read_tab(aff_paths[i]); d <- d[d$archaic_class == wanted,]; if (nrow(d)) { rate <- as.numeric(d$match_rate); rate[!is.finite(rate)] <- NA; cc <- colorRampPalette(c("#ffffcc","#41b6c4","#0c2c84"))(101); rect(to_x(d$chromosome,d$start), i-.35, to_x(d$chromosome,d$end), i+.35, col=cc[pmax(1,pmin(101,1+round(rate*100)))], border=NA) } }
  axis(1, at=offset+lengths/2, labels=paste0("chr",chroms)); axis(2, at=seq_along(refs), labels=refs, las=1); dev.off()
}
draw_intro(snakemake@output[["introgression"]])
nean <- as.character(snakemake@params[["neanderthal_refs"]]); deni <- as.character(snakemake@params[["denisovan_refs"]])
if (length(nean)) draw_affinity(snakemake@output[["neanderthal"]], nean, "neanderthal")
if (length(deni)) draw_affinity(snakemake@output[["denisovan"]], deni, "denisovan")
