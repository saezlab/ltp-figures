#!/usr/bin/Rscript

# (c) Denes Turei 2017 EMBL

source('results.r')
source('proteins.r')

screens_lit_df <- function(with_hptlc = TRUE){
    
    r <- get_results(with_hptlc = with_hptlc)
    m <- master_table(output = FALSE)$m %>%
        rename(
            protein = default_name,
            uhgroup = mammalian_ligand_categories
        ) %>%
        filter(!is.na(uhgroup)) %>%
        unnest(uhgroup = strsplit(uhgroup, ',')) %>%
        select(protein, uhgroup)
    
    d <- rbind(r$a, r$e) %>% select(protein, uhgroup, screen)
    
    p <- d %>%
        select(protein) %>%
        group_by(protein) %>%
        summarise_all(first) %>%
        ungroup()
    
    d <- d %>%
        mutate(
            uhgroup = gsub('-O$', '',
                           gsub('^Lyso', '',
                                as.character(uhgroup)))) %>%
        bind_rows(
            m %>% mutate(screen = 'L')
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
        ungroup() %>%
        inner_join(p, by = c('protein'))
    
    d %>% write_tsv('screens_venn.tsv')
    return(d)
    
}
