#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(lazyeval)
require(stringr)
require(ggrepel)
require(purrr)
require(readr)
require(grid)
require(directlabels)

#
# constants
#

infile_sa  <- 'stats-antonella.csv'
infile_se  <- 'stats-enric.csv'
infile_d   <- 'ltp_domains.tsv'
infile_lit_ligands <- 'binding_properties_plain.csv'

#
# functions
#

proteins_preprocess <- function(){
    
    scoln <- c('protein', 'ionm', 'stat', 'cnt')
    
    ae <- bind_rows(
        read_tsv(infile_a),
        read_tsv(infile_e)
    ) %>%
    group_by(protein) %>%
    select(protein, domain) %>%
    summarise_all(first) %>%
    ungroup()
    
    s <- bind_rows(
            read_tsv(infile_sa, col_names = scoln) %>% mutate(screen = 'A'),
            read_tsv(infile_se, col_names = scoln) %>% mutate(screen = 'E')
        ) %>%
        select(protein, screen) %>%
        group_by(protein) %>%
        mutate(
            screen = paste0(unique(sort(screen)), collapse = ''),
            x = 0.0
        ) %>%
        summarise_all(first) %>%
        ungroup() %>%
        left_join(
            read_tsv(infile_d, col_names = c('protein', 'domain')),
            by = c('protein')
        ) %>%
        arrange(screen) %>%
        mutate(
            screen = factor(
                screen,
                levels = c('A', 'E', 'AE'),
                ordered = TRUE
            ),
            protein = factor(
                protein,
                levels = unique(protein)
            ),
            domain = recode(domain, lipocalin = 'Lipocalin', scp2 = 'Scp2')
        )
        
        for(dom in unique(s$domain)){
            cnt <- dim(s %>% filter(domain == dom))[1]
            s[s$domain == dom,'x'] <- seq(-1 * (cnt - 1) / 2, (cnt - 1) / 2)
        }
    
    return(s)
    
}

proteins_plot <- function(){
    
    fname <- 'proteins'
    
    s <- proteins_preprocess()
    
    p <- ggplot(s, aes(y = domain, color = screen, x = x, label = protein)) +
        # geom_vline(aes(xintercept = 0.0), color = 'gray70') +
        # geom_dotplot(binwidth = 1, dotsize = .2, method = 'histodot', binaxis = 'y', stackdir = 'centerwhole') +
        geom_point(position = position_dodge(width = .9), size = 2.7) +
        geom_text(
            angle = 45,
            color = 'black',
            size = 1.66,
            hjust = 'center',
            vjust = 'middle',
            #nudge_x = .3,
            #nudge_y = .1,
            family = 'DINPro'
        ) +
        facet_grid(domain ~., scales = 'free_y') +
        scale_color_manual(
            values = c(
                'A' = '#6EA945',
                'E' = '#FDD73F',
                'AE'= '#49969A'
            ),
            labels = c(
                'A'  = 'HEK cells',
                'E'  = 'E coli & liposomes',
                'AE' = 'Both'
            ),
            guide = guide_legend(
                title = 'Expressed in'
            )
        ) +
        xlab('Proteins') +
        ylab('Lipid transfer domain families') +
        ggtitle('Lipid transfer proteins in the cargo screenings') +
        theme_minimal() +
        theme(
            panel.background = element_blank(),
            text = element_text(family = 'DINPro'),
            axis.text.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(color = 'black'),
            strip.background = element_blank(),
            strip.text = element_blank()
            #panel.spacing = unit(2.0, 'lines')
        )
    
    cairo_pdf(sprintf('%s.pdf', fname), width = 7, height = 3)
        
        gt = ggplot_gtable(ggplot_build(p))
        gt$layout$clip = "off"
        
        grid.draw(gt)
        
    dev.off()
    
    s %>% write_tsv(sprintf('%s.tsv', fname))
    
}
