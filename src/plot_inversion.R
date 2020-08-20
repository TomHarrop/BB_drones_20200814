#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(circlize)
library(data.table)

# fai_file <- snakemake@input[["ref_fai"]]
# query_fai_file <- snakemake@input[["query_fai"]]
# paf_file <- snakemake@input[["paf"]]
# plot_file <- snakemake@output[["plot"]]
# indiv_name <- snakemake@wildcards[["indiv"]]

# dev
fai_file <- "data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai"
query_fai_file <- "output/027_split/BB31.fa.fai"
paf_file <- "output/040_wga/BB31.paf"
query2_fai_file <- "output/027_split/BB55.fa.fai"
paf2_file <- "output/040_wga/BB55.paf"
plot_file <- "test/inversion.pdf"
# indiv_name <- "BB31"

chr_to_plot <- "NC_037644.1"

# need the names for the FAI
fai_names <- c("name", "length", "offset", "line_bases", "line_bytes")

# only get the main chromosomes from the ref assembly
fai <- fread(fai_file, col.names = fai_names)[
    name %in% chr_to_plot][
        order(length, decreasing = TRUE)]
ref_order <- fai[, unique(name)]

###########
# QUERY 1 #
###########

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

# order by ref contig which they have their biggest hit on?
keep <- filtered_paf[, .I[which.max(alignment_block_length)], by = query][, V1]
biggest_hits <- filtered_paf[keep]
biggest_hits[, ref_fact := factor(target, levels = ref_order)]
setorder(biggest_hits, ref_fact)
query_order <- biggest_hits[!is.na(ref_fact), unique(query)]


###########
# QUERY 2 #
###########


# read the query genome fai
query2_fai <- fread(query2_fai_file,
                    col.names = fai_names)[order(length, decreasing = TRUE)]

# read the paf
mypaf2 <- fread(paf2_file, select = 1:12, col.names = pafnames, fill = TRUE)

# we only want hits > x kb against the ref chroms
block_size_min <- 1e3
score_min <- 30
filtered_paf2 <- mypaf2[alignment_block_length > block_size_min & 
                            mapping_quality > score_min &
                            startsWith(target, "NC_037") ]

# order by ref contig which they have their biggest hit on?
keep2 <- filtered_paf2[, .I[which.max(alignment_block_length)], by = query][, V1]
biggest_hits2 <- filtered_paf2[keep2]
biggest_hits2[, ref_fact := factor(target, levels = ref_order)]
setorder(biggest_hits2, ref_fact)
query_order2 <- biggest_hits2[!is.na(ref_fact), unique(query)]


#########
# PLOTS #
#########

# links
x1 <- filtered_paf[query %in% query_order &
                       query == "contig_00023",
                   .(target, target_start, target_end)]

# subtract the length for query positions to reverse the coordinates
x2 <- filtered_paf[query %in% query_order &
                       query == "contig_00023",
                   .(query,
                     query_length - query_start,
                     query_length - query_end)]
# links
x3 <- filtered_paf2[query %in% query_order2 &
                        query == "contig_00048",
                   .(target, target_start, target_end)]

# subtract the length for query positions to reverse the coordinates
x4 <- filtered_paf2[query %in% query_order2 &
                        query == "contig_00048",
                   .(query,
                     query_length - query_start,
                     query_length - query_end)]


# non-inversion links
x5 <- rbind(filtered_paf[query %in% query_order &
                 query != "contig_00023",
             .(target, target_start, target_end)],
    filtered_paf2[query %in% query_order2 &
                  query != "contig_00048",
              .(target, target_start, target_end)])


x6 <- rbind(
    filtered_paf[query %in% query_order &
                     query != "contig_00023",
                 .(query,
                   query_length - query_start,
                   query_length - query_end)],
    filtered_paf2[query %in% query_order2 &
                      query != "contig_00048",
                  .(query,
                    query_length - query_start,
                    query_length - query_end)]
)

# join the 2 fais
mycols <- viridis::viridis(3)


query_fai_subset <- query_fai[name %in% query_order]
query_fai_subset[, chr_col := mycols[[2]]]
query_fai2_subset <- query2_fai[name %in% query_order2]
query_fai2_subset[, chr_col := mycols[[3]]]

joined_fai <- rbind(fai,
                    rbind(query_fai_subset, query_fai2_subset),
                    fill = TRUE)

myfai <- joined_fai[, .(name, 0, length)]


# generate chr colours
LookupChrColour <- function(x){
    mycol <- joined_fai[name == x, chr_col]
    if(!length(mycol == 1)) {
        return("grey50")
    } else {return(mycol)}
}

# generate link colours
# alpha works if border is off
linkcol <- 
linkcol2 <- 
# linkcol <- fai[x1, on=c(name="target"), chr_col]

# draw the plot

cairo_pdf(plot_file)
# par(xpd = NA)
circos.initializeWithIdeogram(cytoband = myfai[length > 1e3],
                              plotType = c("axis"),
                              chromosome.index = c(ref_order, rev(query_order), rev(query_order2)))

# text(0, 1.05, indiv_name, cex = 1)
# text(0, -1.05, "Reference", cex = 1, )

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

circos.genomicLink(x5, x6,
                   col = ggplot2::alpha("grey", 0.5),
                   border = NA)

circos.genomicLink(x1, x2,
                   col = ggplot2::alpha(mycols[[2]], 0.8),
                   border = NA)
circos.genomicLink(x3, x4,
                   col = ggplot2::alpha(mycols[[3]], 0.8),
                   border = NA)

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

