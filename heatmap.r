#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(dplyr)
require(ggplot2)

source('results.r')
source('clustering.r')

summary_heatmap <- function(wide = FALSE){
    
    pdfname <- ifelse(wide, 'ms-summary-heatmap-w.pdf', 'ms-summary-heatmap.pdf')
    width   <- ifelse(wide, 8, 6.5)
    
    d  <- preprocess0(cluster_proteins = TRUE)
    #cl <- get_protein_ordr(d$c, return_cl = TRUE)
    
    if(wide){
        h <- d$c %>%
            group_by(protein, uhgroup)
    }else{
        h <- d$c %>%
            group_by(protein, hg0)
    }
    
    h <- h %>%
        mutate(
            hg0_lit = any(lit),
            screens = paste0(sort(unique(screen)), collapse = '')
        ) %>%
        mutate(screens = factor(
            screens,
            levels = c('A', 'E', 'AE'),
            ordered = TRUE)
        ) %>%
        summarise_all(first) %>%
        ungroup()
        #mutate(protein = factor(protein, levels = cl, ordered = TRUE))
    
    if(wide){
        p <- ggplot(h, aes(y = protein, x = uhgroup))
    }else{
        p <- ggplot(h, aes(y = protein, x = hg0))
    }
    
    p <- p +
        geom_tile(aes(fill = screens)) +
        geom_point(aes(alpha = lit), color = 'white') +
        scale_alpha_manual(
            guide = guide_legend(title = 'Novelty'),
            values = c(
                'FALSE' = 1.0,
                'TRUE'  = 0.0
            ),
            labels = c(
                'FALSE' = 'Novelty',
                'TRUE'  = 'Known from\nliterature'
            )
        ) +
            scale_fill_manual(
            guide = guide_legend(title = 'Screens'),
            values = c(
                'A' = '#6EA945',
                'E' = '#FCCC06',
                'AE'= '#007B7F'
            ),
            labels = c(
                'A'  = 'HEK cells',
                'E'  = 'E. coli & liposomes',
                'AE' = 'Both'
            )
        ) +
        ggtitle('MS/MS screening of LTP cargoes: summary') +
        xlab('Lipids: main groups') +
        ylab('Proteins') +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro-Medium'),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
            legend.key = element_rect(fill = '#007B7F', color = 'white')
        )
    
    ggsave(pdfname, device = cairo_pdf, width = width, height = 9)
    
    return(h)
    
}
