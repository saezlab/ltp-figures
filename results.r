#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(readr)
require(dplyr)

source('clustering.r')

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_t   <- 'ltp_hptlc_static.tsv'

get_results <- function(with_hptlc = TRUE){
    
    by_t <- c('protein', 'uhgroup', 'screen', 'method')
    
    get_pairs <- function(d){
        
        return(
            d %>%
            filter(cls %in% c('I', 'II')) %>%
            select(protein, uhgroup) %>%
            group_by(protein, uhgroup) %>%
            summarise_all(first) %>%
            ungroup() %>%
            mutate(in_other = TRUE)
        )
        
    }
    
    a <- suppressMessages(read_tsv(infile_a)) %>% mutate(method = 'MS')
    e <- suppressMessages(read_tsv(infile_e)) %>% mutate(method = 'MS')
    
    if(with_hptlc){
        
        t <- suppressMessages(read_tsv(infile_t,
                                    col_names = c('protein', 'uhgroup', 'screen')))
        t <- t %>% mutate(method = 'HPTLC')
        ta <- t %>% filter(screen == 'A')
        te <- t %>% filter(screen == 'E')
        a <- a %>%
            full_join(ta, by = by_t) %>%
            mutate(cls = ifelse(method == 'HPTLC', 'I', cls))
        e <- e %>%
            full_join(te, by = by_t) %>%
            mutate(cls = ifelse(method == 'HPTLC', 'I', cls))
        
    }
    
    apairs <- get_pairs(a)
    epairs <- get_pairs(e)
    
    a <- a %>%
        mutate(lit = lit == 'True') %>%
        left_join(epairs, by = c('protein', 'uhgroup')) %>%
        filter((cls == 'I' | (cls == 'II' & in_other)) & uhgroup != 'P40')
    
    e <- e %>%
        mutate(lit = lit == 'True') %>%
        left_join(apairs, by = c('protein', 'uhgroup')) %>%
        filter((cls == 'I' | (cls == 'II' & in_other)) & uhgroup != 'P40')
    
    result <- list()
    result$a <- a
    result$e <- e
    if(with_hptlc){
        result$t <- t
    }
    
    return(result)
    
}


preprocess0 <- function(cluster_proteins = FALSE, method = 'ward.D2'){
    
    result <- list()
    r <- get_results()
    
    # preprocessing
    ae <- bind_rows(r$a, r$e) %>%
        group_by(protein, uhgroup) %>%
        mutate(screens = paste0(sort(unique(screen)), collapse = '')) %>%
        ungroup() %>%
        group_by(protein) %>%
        group_by(protein, ionm, screen, uhgroup, cls) %>%
        mutate(ihg = sum(as.numeric(intensity))) %>%
        ungroup() %>%
        group_by(protein, ionm, screen) %>%
        mutate(itotal = sum(as.numeric(ihg))) %>%
        ungroup() %>%
        group_by(protein, ionm, screen, uhgroup, cls) %>%
        mutate(irel = ihg / itotal) %>%
        mutate(lirel = log10(1 + irel)) %>%
        ungroup() %>%
        mutate(hg0 = gsub('-O$', '', gsub('^Lyso', '', as.character(uhgroup)))) %>%
        mutate(grp = ifelse(
            hg0 %in% gpl, 'GPL', ifelse(
            hg0 %in% gl, 'GL', ifelse(
            hg0 %in% fa, 'FA', ifelse(
            hg0 %in% ch, 'CH', ifelse(
            hg0 %in% sph, 'SPH', ifelse(
            hg0 %in% vit, 'VIT', NA
        ))))))) %>%
        mutate(hgcc0 = ifelse(is.na(hgcc), uhgroup, hgcc)) %>%
        group_by(grp) %>%
        mutate(cnt_grp = n()) %>%
        ungroup() %>%
        group_by(uhgroup) %>%
        mutate(cnt_hg = n()) %>%
        ungroup() %>%
        group_by(protein) %>%
        mutate(cnt_pro = n()) %>%
        ungroup()
    
    if(cluster_proteins){
        
        protein_ordr <- get_protein_ordr(ae, method = method)
        
    }
    
    l <- ae %>%
        select(uhgroup, hg0, grp, hgcc0, screen, cnt_grp, cnt_hg) %>%
        group_by(hg0) %>%
        mutate(screens = paste0(sort(unique(screen)), collapse = '')) %>%
        ungroup() %>%
        group_by(hgcc0) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(grp = factor(grp, levels = grp_ordr, ordered = TRUE)) %>%
        arrange(grp, hg0, uhgroup, hgcc0) %>%
        mutate(
            hg0 = factor(hg0, levels = unique(hg0), ordered = TRUE),
            uhgroup = factor(uhgroup, levels = unique(uhgroup), ordered = TRUE),
            hgcc0 = factor(hgcc0, levels = unique(hgcc0), ordered = TRUE)
        )
    
    p <- ae %>%
        select(protein, screen, cnt_pro) %>%
        group_by(protein) %>%
        mutate(screens = paste0(sort(unique(screen)), collapse = '')) %>%
        summarise_all(first) %>%
        ungroup()
    
    if(cluster_proteins){
        
        # ordering by proteins
        p <- p %>%
            mutate(
                protein = factor(
                    protein,
                    levels = protein_ordr,
                    ordered = TRUE
                )
            ) %>%
            arrange(protein)
        
    }else{
        
        # ordering by screens
        p <- p %>%
            mutate(
                screens = factor(
                    screens,
                    levels = scr_ordr,
                    ordered = TRUE
                )
            ) %>%
            arrange(screens, protein) %>%
            mutate(
                protein = factor(
                    protein,
                    levels = unique(protein),
                    ordered = TRUE
                )
            )
        
    }
    
    # full list of connections:
    result$c <- ae
    # all ligands:
    result$l <- l
    # all proteins:
    result$p <- p
    # proteins order:
    result$pordr <- switch(
        cluster_proteins,
        protein_ordr,
        unique(p$protein)
    )
    
    return(result)
    
}
