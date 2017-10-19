#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(readr)
require(dplyr)

source('clustering.r')
source('lipid_classes.r')
source('proteins.r')

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_t   <- 'ltp_hptlc_static.tsv'


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


pre_preprocess <- function(d, l, p){
    
    return(
        d %>%
        mutate(
            lit = lit == 'True',
            hg0 = gsub('-O$', '', gsub('^Lyso', '', as.character(uhgroup)))
        ) %>%
        mutate(
            hg0 = recode(
                hg0,
                CH = 'Sterols',
                PIP = 'PIPs'
            )
        ) %>%
        left_join(l, by = c('protein', 'hg0')) %>%
        mutate(lit0 = !is.na(lit0)) %>%
        left_join(p, by = c('protein', 'uhgroup')) %>%
        filter((cls == 'I' | (cls == 'II' & in_other)) & uhgroup != 'P40')
    )
    
}


get_results0 <- function(with_hptlc = TRUE){
    
    by_t <- c('protein', 'uhgroup', 'screen', 'method')
    
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
    
    result <- list()
    result$a <- a
    result$e <- e
    
    return(result)

}

get_results <- function(with_hptlc = TRUE){
    
    r <- get_results0(with_hptlc = with_hptlc)
    a <- r$a
    e <- r$e
    
    l <- literature_ligands()
    apairs <- get_pairs(a)
    epairs <- get_pairs(e)
    
    a <- pre_preprocess(a, l, epairs)
    e <- pre_preprocess(e, l, apairs)
    
    result <- list()
    result$a <- a
    result$e <- e
    if(with_hptlc){
        result$t <- t
    }
    
    return(result)
    
}


preprocess0 <- function(cluster_proteins = FALSE, method = 'ward.D2', with_hptlc = TRUE){
    
    result <- list()
    r <- get_results(with_hptlc = with_hptlc)
    
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
    
    r <- get_results0()
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

literature_ligands <- function(){
    
    m <- master_table(output = FALSE)$m %>%
        rename(
            protein = default_name,
            hg0 = mammalian_ligand_categories
        ) %>%
        filter(!is.na(hg0)) %>%
        unnest(hg0 = strsplit(hg0, ',')) %>%
        select(protein, hg0) %>%
        mutate(lit0 = TRUE)
    
}
