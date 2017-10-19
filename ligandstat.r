#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)

source('results.r')

infile_sa  <- 'stats-antonella.csv'
infile_se  <- 'stats-enric.csv'

domain_colors <- function(d, col){
    
    return(col[unique(as.character(d$protein))])
    
}

read_stats <- function(){
    
    coln <- c('protein', 'ionm', 'stat', 'val')
    
    sa <- suppressMessages(read_tsv(infile_sa, col_names = coln))
    se <- suppressMessages(read_tsv(infile_se, col_names = coln))
    
    return(list(sa = sa, se = se))
    
}


domains_assign_colors <- function(d){
    
    dcolors <- c('#3A7AB3', '#608784', '#03928C', '#CF5836', '#7575BE',
                '#D6708B', '#65B9B9', '#69B3D6', '#C441B3', '#9B998D')
    
    names(dcolors) <- sort(unique(na.omit(d$domain)))
    
    dprotein <- list()
    for(i in 1:dim(d)[1]){
        dprotein[[as.character(d$protein[i])]] <- as.character(d$domain[i])
    }
    
    pcolors <- as.character(dcolors[as.character(dprotein)])
    names(pcolors) <- names(dprotein)
    
    return(pcolors)
    
}


numof_hgs <- function(){
    
    m <- master_table(output = FALSE)
    d <- preprocess0()$c
    
    dcol <- domains_assign_colors(d)
    
    s <- d %>%
        group_by(protein) %>%
        mutate(
            cnt = n(),
            cnt_hg = n_distinct(hg0)
        ) %>%
        summarise_all(first) %>%
        ungroup() %>%
        arrange(desc(cnt_hg), domain, protein) %>%
        select(protein, domain, cnt, cnt_hg) %>%
        mutate(protein = factor(protein, levels = unique(protein), ordered = TRUE))
    
    p <- ggplot(s, aes(x = protein, y = cnt_hg)) +
        geom_bar(fill = 'black', stat = 'identity') +
        theme_linedraw() +
        xlab('Proteins') +
        ylab('Number of main lipid categories') +
        theme(
            text = element_text(family = 'DINPro'),
            axis.text.x = element_text(
                angle = 90, vjust = 0.5, size = 8, hjust = 1,
                color = domain_colors(s, dcol)
            )
        )
    
    ggsave('numof_different_hgs.pdf', device = cairo_pdf, width = 6, height = 4)
    
    s <- s %>%
        arrange(desc(cnt), domain, protein) %>%
        mutate(protein = factor(protein, levels = unique(protein), ordered = TRUE))
    
    p <- ggplot(s, aes(x = protein, y = cnt)) +
        geom_bar(fill = 'black', stat = 'identity') +
        theme_linedraw() +
        xlab('Proteins') +
        ylab('Number of identified features') +
        theme(
            text = element_text(family = 'DINPro'),
            axis.text.x = element_text(
                angle = 90, vjust = 0.5, size = 8, hjust = 1,
                color = domain_colors(s, dcol)
            ),
            panel.grid = element_blank()
        )
    
    ggsave('numof_features.pdf', device = cairo_pdf, width = 6, height = 4)
    
    p <- ggplot(s, aes(x = cnt, y = cnt_hg, label = protein)) +
    geom_text(color = domain_colors(s, dcol), angle = 30) +
        scale_x_log10() +
        scale_y_log10() +
        theme_linedraw() +
        xlab('Number of identified features (log)') +
        ylab('Number of different headgroups (log)') +
        theme(
            text = element_text(family = 'DINPro')
        )
    
    ggsave('numof_hgs_features.pdf', device = cairo_pdf, width = 6, height = 6)
    
}
