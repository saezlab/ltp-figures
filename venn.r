#!/usr/bin/Rscript

# (c) Denes Turei 2017 EMBL

infile_lit_ligands <- 'binding_properties_plain.csv'

source('results.r')

screens_lit_df <- function(with_hptlc = FALSE){
    
    r <- get_results()
    l <- suppressMessages(read_tsv(infile_lit_ligands,
                                   col_names = c('protein', 'uhgroup')))
    
    d <- rbind(r$a, r$e) %>% select(protein, uhgroup, screen)
    
    if(with_hptlc){
        
        t <- r$t %>% mutate(screen = 'T')
        d <- rbind(d, t)
        
    }
    
    d <- d %>%
        mutate(
            uhgroup = gsub('-O$', '',
                           gsub('^Lyso', '',
                                as.character(uhgroup)))) %>%
        bind_rows(
            l %>% mutate(screen = 'L')
        ) %>%
        mutate(
            uhgroup = recode(uhgroup,
                KCH = 'CH',
                HCH = 'CH',
                PCH = 'CH',
                `PAF, LPAF` = 'PAF',
                PHCH = 'CH',
                PUFA = 'FA',
                VLCFA = 'FA',
                LCFA = 'FA',
                GM2 = 'GM',
                GM3 = 'GM',
                GM1 = 'GM',
                TAIVA = 'IVA',
                LPS = 'PS',
                CCA = 'CA',
                CHS = 'CH'
            )
        ) %>%
        group_by(protein, uhgroup) %>%
        mutate(screen = paste0(unique(sort(screen)), collapse = '')) %>%
        mutate(in_lit = grepl('L', screen, fixed = TRUE)) %>%
        summarise_all(first) %>%
        ungroup()
    
    d %>% write_tsv('screens_venn.tsv')
    return(d)
    
}
