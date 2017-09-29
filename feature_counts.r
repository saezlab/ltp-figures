#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(tidyr)
require(lazyeval)
require(stringr)
require(ggrepel)

#
# constants
#

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_sa  <- 'stats-antonella.csv'
infile_se  <- 'stats-enric.csv'
stat_ord   <- c('all', 'valid', 'pp', 'c12')

#
# functions
#

preprocess <- function(){
    
    results_df <- function(fname){
        
        return(
            read.table(fname, sep = '\t', header = TRUE) %>%
            group_by(protein, ionm, screen) %>%
            summarise(cnt = n_distinct(id)) %>%
            ungroup() %>%
            group_by(protein, ionm) %>%
            summarise_all(first) %>%
            ungroup() %>%
            droplevels() %>%
            complete(protein, ionm, screen) %>%
            mutate(
                stat = 'pp',
                cnt  = ifelse(is.na(cnt), 0, cnt)
            ) %>%
            select(protein, ionm, stat, cnt, screen)
        )
        
    }
    
    class12_df <- function(fname, other){
        
        return(
            read.table(fname, sep = '\t', header = TRUE) %>%
            left_join(
                read.table(other, sep = '\t', header = TRUE) %>%
                filter(cls %in% c('I', 'II')) %>%
                mutate(otherscreen = ifelse(screen == 'A', 'E', 'A')) %>%
                select(protein, uhgroup, otherscreen) %>%
                group_by(protein, uhgroup, otherscreen) %>%
                summarise_all(first) %>%
                rename(screen = otherscreen) %>%
                mutate(other_confirmed = TRUE),
                by = c('protein', 'uhgroup', 'screen')
            ) %>%
            filter(cls == 'I' | (cls == 'II' & !is.na(other_confirmed))) %>%
            group_by(protein, ionm, screen) %>%
            summarise(cnt = n_distinct(id)) %>%
            ungroup() %>%
            group_by(protein, ionm) %>%
            summarise_all(first) %>%
            ungroup() %>%
            droplevels() %>%
            complete(protein, ionm, screen) %>%
            mutate(
                stat   = 'c12',
                cnt    = ifelse(is.na(cnt), 0, cnt)
            ) %>%
            select(protein, ionm, stat, cnt, screen)
        )
        
    }
    
    unified_df <- function(...){
        
        return(
            bind_rows(...) %>%
            complete(protein, stat, screen, ionm)
        )
        
    }
    
    sa  <- read.table(infile_sa, sep = '\t', header = FALSE)
    se  <- read.table(infile_se, sep = '\t', header = FALSE)
    names(sa) <- c('protein', 'ionm', 'stat', 'cnt')
    names(se) <- c('protein', 'ionm', 'stat', 'cnt')
    sa <- sa %>% mutate(screen = 'A')
    se <- se %>% mutate(screen = 'E')
    
    return(
        bind_rows(
            unified_df(
                sa,
                results_df(infile_a),
                class12_df(infile_a, infile_e)
            ),
            unified_df(
                se,
                results_df(infile_e),
                class12_df(infile_e, infile_a)
            )
        ) %>%
        mutate(scr_ionm = paste0(screen, ionm))
    )
    
}

preprocess_stage2 <- function(){
    
    s <- preprocess()
    
    return(
        s %>%
        filter(stat %in% stat_ord) %>%
        mutate(
            stat = factor(stat, levels = stat_ord, ordered = TRUE),
            grp  = paste0(screen, ionm, protein),
            cnt  = ifelse(is.na(cnt), 0.0, cnt)
        )
    )
    
}

counts_plot <- function(facets = TRUE){
    
    pdfwidth <- ifelse(facets, 4.5, 1.75)
    pdfname  <- ifelse(facets, 'ms_feature_counts.pdf', 'ms_feature_counts_1.pdf')
    ttlsize  <- ifelse(facets, 12, 7)
    title    <- ifelse(
        facets,
        'Number of features along the selection process',
        'Number of features\nalong the\nselection process'
    )
    
    s <- preprocess_stage2()
    
    p <- 
        
    
    if(facets){
        
        p <- ggplot(s, aes(x = stat, y = cnt + 1, group = protein)) +
            facet_grid(
                . ~ scr_ionm,
                labeller = as_labeller(c(
                    Aneg = 'In vivo, negative ion mode',
                    Apos = 'In vivo, positive ion mode',
                    Eneg = 'In vitro, negative ion mode',
                    Epos = 'In vitro, positive ion mode'
            )))
        
    }else{
        
        p <- ggplot(s, aes(x = stat, y = cnt + 1, group = grp))
        
    }
    
    p <- p +
        geom_line(size = .1, alpha = .2) +
        scale_y_log10(
            breaks = c(10, 100, 1000, 10000)
        ) +
        scale_x_discrete(
            labels = c(
                all   = 'All',
                valid = 'Quality\nfiltered',
                pp    = 'Fits protein\nprofile',
                c12   = 'Identified'
            ),
            expand = c(0.10, 0)
        ) +
        ggtitle(title) +
        xlab('Stage of selection') +
        ylab('Number of features (log)') +
        theme_linedraw() +
        theme(
            text = element_text(family = 'DINPro'),
            axis.text.x = element_text(angle = 90, vjust = .5, size = 7, hjust = 1),
            strip.text = element_text(size = 5),
            axis.title = element_text(size = 9),
            plot.title = element_text(size = ttlsize)
        )
    
    ggsave(pdfname, device = cairo_pdf, width = pdfwidth, height = 3)
    
}

where_greater <- function(d, stata, statb){
    
    return(
        d %>%
        filter(stat == stata) %>%
        inner_join(
            d %>% filter(stat == statb),
            by = c('protein', 'scr_ionm'),
            suffix = c(stata, statb)
        ) %>%
        filter_(paste0('cnt', stata, '>', 'cnt', statb))
    )
    
}
