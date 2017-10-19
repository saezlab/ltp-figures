#!/usr/bin/Rscript

# (c) Denes Turei 2017 EMBL

source('results.r')

screens_lit_df <- function(with_hptlc = TRUE){
    
    r <- get_results(with_hptlc = with_hptlc)
    
    d <- rbind(r$a, r$e) %>% select(protein, hg0, screen)
    l <- literature_ligands() %>% mutate(screen = 'L')
    
    p <- d %>%
        select(protein) %>%
        group_by(protein) %>%
        summarise_all(first) %>%
        ungroup()
    
    d <- d %>%
        bind_rows(l) %>%
        group_by(protein, hg0) %>%
        select(protein, hg0, screen) %>%
        mutate(screen = paste0(unique(sort(screen)), collapse = '')) %>%
        mutate(in_lit = grepl('L', screen, fixed = TRUE)) %>%
        summarise_all(first) %>%
        ungroup() %>%
        inner_join(p, by = c('protein'))
    
    d %>% write_tsv('screens_venn.tsv')
    return(d)
    
}
