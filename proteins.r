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

source('results.r')

#
# constants
#

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_sa  <- 'stats-antonella.csv'
infile_se  <- 'stats-enric.csv'
infile_d   <- 'ltp_domains.tsv'
infile_lit_ligands <- 'binding_properties_plain.csv'
infile_master1 <- 'master_part1.csv'

#
# functions
#

proteins_preprocess <- function(){
    
    scoln <- c('protein', 'ionm', 'stat', 'cnt')
    
    ae <- bind_rows(
        suppressMessages(read_tsv(infile_a)),
        suppressMessages(read_tsv(infile_e))
    ) %>%
    group_by(protein) %>%
    select(protein, domain) %>%
    summarise_all(first) %>%
    ungroup()
    
    s <- bind_rows(
            suppressMessages(
                read_tsv(infile_sa, col_names = scoln) %>% mutate(screen = 'A')
            ),
            suppressMessages(
                read_tsv(infile_se, col_names = scoln) %>% mutate(screen = 'E')
            )
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
            suppressMessages(
                read_tsv(infile_d, col_names = c('protein', 'domain'))
            ),
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
            domain = recode(
                domain,
                lipocalin = 'Lipocalin',
                scp2 = 'Scp2',
                IP_trans = 'IP-trans',
                LBP_BPI_CETP = 'LBP-BPI-CETP'
            )
        )
        
        for(dom in unique(s$domain)){
            cnt <- dim(s %>% filter(domain == dom))[1]
            s[s$domain == dom,'x'] <- seq(-1 * (cnt - 1) / 2, (cnt - 1) / 2)
            # s[s$domain == dom,'x'] <- seq(1:cnt)
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
            angle = 35,
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
    
    cairo_pdf(sprintf('%s.pdf', fname), width = 7, height = 2.5)
        
        gt = ggplot_gtable(ggplot_build(p))
        gt$layout$clip = "off"
        
        grid.draw(gt)
        
    dev.off()
    
    s %>% write_tsv(sprintf('%s.tsv', fname))
    
}

master_table <- function(output = TRUE){
    
    hdr <- c(
        'Default name',
        'Synonyms',
        'UniProt ID',
        'LTD family',
        'Literature ligands',
        'Mammalian ligands main categories',
        'Non-mammalian ligands',
        'Tested in vitro',
        'Tested in vivo',
        'Number of ligands in final result in vivo (unique m/z)',
        'Ligands identified in vivo',
        'Number of ligands in final result in vitro (unique m/z)',
        'Ligands identified in vitro'
    )
    
    r <- get_results()
    p <- proteins_preprocess() %>%
        mutate(in_invivo = screen %in% c('A', 'AE'),
               in_invitro = screen %in% c('E', 'AE')) %>%
        select(protein, in_invitro, in_invivo)
    m1 <- suppressMessages(read_tsv(infile_master1))
    
    ra <- r$a %>%
        group_by(protein) %>%
        mutate(invivo_count = n(), invivo_binders = paste0(sort(unique(uhgroup)), collapse = ',')) %>%
        summarise_all(first) %>%
        select(protein, invivo_count, invivo_binders)
    
    re <- r$e %>%
        group_by(protein) %>%
        mutate(invitro_count = n(), invitro_binders = paste0(sort(unique(uhgroup)), collapse = ',')) %>%
        summarise_all(first) %>%
        select(protein, invitro_count, invitro_binders)
    
    m <- m1 %>%
        mutate(protein = default_name) %>%
        left_join(p, by = c('protein')) %>%
        mutate(in_invitro = ifelse(is.na(in_invitro), FALSE, in_invitro),
               in_invivo = ifelse(is.na(in_invivo), FALSE, in_invivo)) %>%
        left_join(ra, by = c('protein')) %>%
        left_join(re, by = c('protein')) %>%
        mutate(
            invitro_count   = ifelse(is.na(invitro_count), 0, invitro_count),
            invivo_count   = ifelse(is.na(invivo_count), 0, invivo_count)
        ) %>%
        arrange(desc(invivo_count), desc(invitro_count), desc(in_invitro), desc(in_invivo)) %>%
        select(default_name:non_mammalian_ligands, in_invivo, in_invitro, invivo_count:invitro_binders)
    
    if(output){
        
        write_tsv(as.data.frame(t(hdr)), 'master_part2.tsv', col_names = FALSE)
        write_tsv(m %>% filter(in_invivo | in_invitro), 'master_part2.tsv', append = TRUE)
        write('\nNon tested LTPs:', file = 'master_part2.tsv', append = TRUE)
        write_tsv(as.data.frame(t(hdr)), 'master_part2.tsv', append = TRUE)
        write_tsv(m %>% filter(!in_invivo & !in_invitro), 'master_part2.tsv', append = TRUE)
        
    }
    
    result <- list()
    result$a  <- r$a
    result$ra <- ra
    result$re <- re
    result$m  <- m
    result$p  <- p
    
    invisible(return(result))
    
}
