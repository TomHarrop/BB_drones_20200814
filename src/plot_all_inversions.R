

library(cowplot)
library(magick)

inversion_pdfs <- list.files("output/040_wga",
                             pattern = "NC_037644.1.pdf",
                             full.names = TRUE)

names(inversion_pdfs) <- sapply(basename(inversion_pdfs), function(x)
    unlist(strsplit(x, ".", fixed = TRUE))[[1]])


inversion_tbs <- lapply(inversion_pdfs, image_read_pdf)

fh <- grid::convertUnit(grid::unit(227, "pt"), "in", valueOnly = TRUE)
fw <- grid::convertUnit(grid::unit(398, "pt"), "in", valueOnly = TRUE)

cairo_pdf("test/inversions.pdf", width = fw, height = fh)

plot_grid(plotlist = lapply(inversion_tbs, function(x)
                 ggdraw() + draw_image(x)),
          nrow = 2)
dev.off() 
