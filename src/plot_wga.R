#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(circlize)
library(data.table)

fai_file <- snakemake@input[["ref_fai"]]
query_fai_file <- snakemake@input[["query_fai"]]
paf_file <- snakemake@input[["paf"]]
plot_file <- snakemake@output[["plot"]]
indiv_name <- snakemake@wildcards[["indiv"]]

# dev
# fai_file <- "/home/tom/Containers/assembly-utils/img/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai"
# query_fai_file <- "output/027_split/BB31.fa.fai"
# paf_file <- "output/040_wga/BB31.paf"
# plot_file <- "test/BB31.pdf"
# indiv_name <- "BB31"

# need the names for the FAI
fai_names <- c("name", "length", "offset", "line_bases", "line_bytes")

# only get the main chromosomes from the ref assembly
fai <- fread(fai_file, col.names = fai_names)[
    startsWith(name, "NC_037")][
        order(length, decreasing = TRUE)]
ref_order <- fai[, unique(name)]

# read the query genome fai
query_fai <- fread(query_fai_file,
                   col.names = fai_names)[order(length, decreasing = TRUE)]

# read the paf
pafnames <- c(
    "query",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "target",
    "target_length",
    "target_start",
    "target_end",
    "residue_matches",
    "alignment_block_length",
    "mapping_quality")
mypaf <- fread(paf_file, select = 1:12, col.names = pafnames, fill = TRUE)

# we only want hits > x kb against the ref chroms
block_size_min <- 1e4
score_min <- 30
filtered_paf <- mypaf[alignment_block_length > block_size_min & 
                          mapping_quality > score_min &
                          startsWith(target, "NC_037") ]

# order by query contigs with the highest percent of hits?
query_hitlength <- filtered_paf[, .(
    total_hitlength = sum(alignment_block_length)),
    by = .(query, query_length)]

query_hitlength[, frac_hitlength := total_hitlength / query_length]
query_order <- query_hitlength[order(frac_hitlength, decreasing = TRUE),
                               unique(query)]

# order by ref contig which they have their biggest hit on?
keep <- filtered_paf[, .I[which.max(alignment_block_length)], by = query][, V1]
biggest_hits <- filtered_paf[keep]
biggest_hits[, ref_fact := factor(target, levels = ref_order)]
setorder(biggest_hits, ref_fact)
query_order <- biggest_hits[, unique(query)]

# links
x1 <- filtered_paf[query %in% query_order,
                   .(target, target_start, target_end)]

# subtract the length for query positions to reverse the coordinates
x2 <- filtered_paf[query %in% query_order,
                   .(query,
                     query_length - query_start,
                     query_length - query_end)]


# join the 2 fais
query_fai_subset <- query_fai[name %in% query_order]
joined_fai <- rbind(fai,
                    query_fai_subset,
                    fill = TRUE)

myfai <- joined_fai[, .(name, 0, length)]


# generate chr colours
fai[, chr_col := viridis::viridis(.N)]
LookupChrColour <- function(x){
    mycol <- fai[name == x, chr_col]
    if(!length(mycol == 1)) {
        return("grey50")
    } else {return(mycol)}
}

# generate link colours
# alpha works if border is off
# linkcol <- fai[x1, on=c(name="target"), ggplot2::alpha(chr_col, 1)]
linkcol <- fai[x1, on=c(name="target"), chr_col]

# draw the plot

cairo_pdf(plot_file)
par(xpd = NA)
circos.initializeWithIdeogram(cytoband = myfai[length > 1e3],
                              plotType = c("axis"),
                              chromosome.index = c(ref_order, rev(query_order)))

text(0, 1.05, indiv_name, cex = 1)
text(0, -1.05, "Reference", cex = 1, )

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    # named chromosomes - don't fit
    # chr = CELL_META$sector.index
    chr = ""
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0,
                xlim[2], 1,
                col = LookupChrColour(CELL_META$sector.index))
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)},
    track.height = 0.1, bg.border = NA)

circos.genomicLink(x1, x2,
                   col = linkcol,
                   border = linkcol)
dev.off()

sessionInfo()

# testing - manual ribbons
# i <- x2[, .I[which.max(V2 - V3)]]
# x1[i]
# x2[i]
# 
# circos.genomicLink(x1[i], x2[i], col = linkcol[i], border = NA)
# circos.link("NC_037648.1", c(4952021, 7309184),
#             "contig_00035", c(9070068, 11426699))

