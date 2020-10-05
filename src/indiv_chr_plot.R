#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(circlize)
library(data.table)

# files
fai_file <- snakemake@input[["ref_fai"]]
query_fai_file <- snakemake@input[["query_fai"]]
paf_file <- snakemake@input[["paf"]]
plot_file <- snakemake@output[["plot"]]

# indiv_name <- snakemake@wildcards[["indiv"]]
chr_to_plot <- snakemake@wildcards[["chr"]]

block_size_min <- snakemake@params[["block_size"]]
score_min <- snakemake@params[["score"]]

# dev
# fai_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai"
# query_fai_file <- "output/027_oriented/BB28.fa.fai"
# paf_file <- "output/040_wga/BB28.paf"
# plot_file <- "test/inversion.pdf"

# params
# chr_to_plot <- "NC_037644.1"
# block_size_min <- 1e4
# score_min <- 30
# indiv_name <- "BB28"

# need the names for the FAI
fai_names <- c("name", "length", "offset", "line_bases", "line_bytes")

# only get the main chromosomes from the ref assembly
fai <- fread(fai_file, col.names = fai_names)[
    name %in% chr_to_plot][
        order(length, decreasing = TRUE)]
ref_order <- fai[, unique(name)]


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


# READ QUERY

# read the query genome fai
query_fai <- fread(query_fai_file,
                   col.names = fai_names)[order(length, decreasing = TRUE)]

# read the paf
mypaf <- fread(paf_file, select = 1:12, col.names = pafnames, fill = TRUE)

# we only want hits > x kb against the ref chroms
filtered_paf <- mypaf[alignment_block_length > block_size_min & 
          mapping_quality > score_min &
          target == chr_to_plot]

# order by ref contig which they have their biggest hit on?
keep <- filtered_paf[, .I[which.max(alignment_block_length)], by = query][, V1]
biggest_hits <- filtered_paf[keep]
biggest_hits[, ref_fact := factor(target, levels = ref_order)]
setorder(biggest_hits, ref_fact)
query_order <- biggest_hits[!is.na(ref_fact), unique(query)]

# GENERATE PLOT

# detect inverted links
targetpaf <- filtered_paf[query %in% query_order]
noninvpaf <- targetpaf[strand == "+"]
invpaf <- targetpaf[strand == "-"]

# join the 2 fais
mycols <- viridis::viridis(4)

fai[name == chr_to_plot, chr_col := mycols[[1]]]

query_fai_subset <- query_fai[name %in% query_order]
query_fai_subset[, chr_col := mycols[[2]]]

joined_fai <- rbind(fai,
                    query_fai_subset,
                    fill = TRUE)

myfai <- joined_fai[, .(name, 0, length)]


# generate chr colours
LookupChrColour <- function(x){
    mycol <- joined_fai[name == x, chr_col]
    if(!length(mycol == 1)) {
        return("grey50")
    } else {return(mycol)}
}


cairo_pdf(plot_file,
          family = "Lato")
# par(xpd = NA)
circos.initializeWithIdeogram(cytoband = myfai[length > 1e3],
                              plotType = c("axis"),
                              chromosome.index = c(ref_order,
                                                   rev(query_order)))

text(0, 1.05, indiv_name, cex = 1, col = mycols[[2]])
text(0, -1.05, "Reference", cex = 1, col = mycols[[1]])

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    # named chromosomes - don't fit
    chr = CELL_META$sector.index
    # chr = ""
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0,
                xlim[2], 1,
                col = LookupChrColour(CELL_META$sector.index))
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)},
    track.height = 0.1, bg.border = NA)


circos.genomicLink(targetpaf[, .(target, target_start, target_end)],
                   targetpaf[, .(query,
                                 query_length - query_start,
                                 query_length - query_end)],
                   col = ggplot2::alpha(mycols[[3]], 0.8),
                   border = NA)

# circos.genomicLink(noninvpaf[, .(target, target_start, target_end)],
#                    noninvpaf[, .(query,
#                                  query_length - query_start,
#                                  query_length - query_end)],
#                    col = ggplot2::alpha("grey", 0.5),
#                    border = NA)
# 
# 
# circos.genomicLink(invpaf[, .(target, target_start, target_end)],
#                    invpaf[, .(query,
#                               query_length - query_start,
#                               query_length - query_end)],
#                    col = ggplot2::alpha(mycols[[2]], 0.8),
#                    border = NA)


dev.off()

sessionInfo()
