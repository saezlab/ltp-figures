#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(tibble)
source('results.r')

gpl <- c('PE', 'PC', 'PG', 'PA', 'PI', 'PS', 'PG/BMP', 'BMP')
gl <- c('DAG', 'TAG', 'MAG')
fa <- c('FA')
ch <- c('CH', 'HCH')
sph <- c('SM', 'Sph', 'SphP', 'Cer', 'CerP', 'HexCer',
         'HexCerOH', 'Hex2Cer', 'Hex2CerOH', 'SHexCer')
vit <- c('VE', 'VA', 'VD')

grp_ordr <- c('GPL', 'GL', 'FA', 'SPH', 'CH', 'VIT')
scr_ordr <- c('A', 'AE', 'E')

karyotype_df <- function(){
    
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
        ungroup() %>%
        mutate(screens = factor(screens, levels = scr_ordr, ordered = TRUE)) %>%
        arrange(screens, protein) %>%
        mutate(protein = factor(protein, levels = unique(protein), ordered = TRUE))
    
    # full list of connections:
    result$c <- ae
    # all ligands:
    result$l <- l
    # all proteins:
    result$p <- p
    
    return(result)
    
}

karyotype_file <- function(){
    
    result <- list()
    
    d <- karyotype_df()
    
    ptotal <- sum(d$p$cnt_pro)
    ltotal <- sum(
        (
            d$l %>%
            group_by(grp) %>%
            summarise(cnt_grp = first(cnt_grp))
        )$cnt_grp
    )
    total <- ptotal + ltotal
    
    p <- d$p %>%
        add_column(
            'cr' = 'chr',
            'pr' = '-',
            protein2 = d$p$protein,
            .after = 1
        ) %>%
        mutate(
            lb = paste0('pro', screens),
            rn  = 1:nrow(d$p)
        ) %>%
        mutate(
            from = lag(cumsum(cnt_pro), 1),
            to   = cumsum(cnt_pro)
        ) %>%
        mutate(from = ifelse(is.na(from), 0, from)) %>%
        select(cr, pr, protein2, protein, protein, from, to, lb)
    
    l <- d$l %>%
        group_by(grp) %>%
        summarise_all(first)
    
    l <- l %>%
        add_column(
            'cr'  = 'chr',
            'pr' = '-',
            'protein'  = '',
            'protein2' = '',
            .after = 1
            ) %>%
        mutate(
            protein  = grp,
            protein2 = grp,
            lb = paste0('grp', grp),
            rn  = 1:nrow(l)
            ) %>%
        mutate(
            from = lag(cumsum(cnt_grp), 1) + max(p$to),
            to   = cumsum(cnt_grp) + max(p$to)
            ) %>%
        mutate(from = ifelse(is.na(from), max(p$to), from)) %>%
        select(cr, pr, protein2, protein, from, to, lb)
    
    all_chr <- suppressWarnings(bind_rows(p, l))
    
    suppressMessages(
        all_chr %>%
        write_delim('circos1/karyotype.txt', col_names = FALSE)
    )
    
    # karyotype ready
    
    # generating links
    
    coo <- d$c %>%
        #group_by(protein, hgcc0) %>%
        #summarise_all(first) %>%
        #ungroup() %>%
        left_join(
            p %>% select(protein, from, to),
            by = c('protein'),
            suffix = c('', '_pro')
        ) %>%
        left_join(
            l %>% select(protein, from, to),
            by = c('grp' = 'protein'),
            suffix = c('', '_grp')
        ) %>%
        arrange(grp, uhgroup, hgcc0) %>%
        mutate(
            grp = factor(grp, levels = grp_ordr, ordered = TRUE),
            uhgroup = factor(uhgroup, levels = unique(uhgroup), ordered = TRUE),
            hgcc0 = factor(hgcc0, levels = unique(hgcc0), ordered = TRUE)
        ) %>%
        group_by(protein) %>%
        mutate(pro_offset = from +  1:n()) %>%
        ungroup() %>%
        group_by(grp) %>%
        mutate(grp_offset = from_grp + 1:n()) %>%
        ungroup() %>%
        mutate(
            pro_offset2 = pro_offset - 1,
            grp_offset2 = grp_offset - 1
        ) %>%
        mutate(
            pro_offset = format(pro_offset, scientific = FALSE, justify = 'none', trim = TRUE),
            grp_offset = format(grp_offset, scientific = FALSE, justify = 'none', trim = TRUE),
            pro_offset2 = format(pro_offset2, scientific = FALSE, justify = 'none', trim = TRUE),
            grp_offset2 = format(grp_offset2, scientific = FALSE, justify = 'none', trim = TRUE)
        )
    
    li <- coo %>%
        select(protein, pro_offset2, pro_offset, grp, grp_offset2, grp_offset)
    
    suppressMessages(li %>% write_delim('circos1/links.txt',
                                        col_names = FALSE))
    
    # generating individual labels for lipids
    
    la <- coo %>%
        select(grp, grp_offset2, grp_offset, hgcc)
    
    suppressMessages(la %>% write_delim('circos1/labels.txt',
                                        col_names = FALSE))
    
    # returning data frames
    
    result$chr <- all_chr
    result$li  <- li
    
    invisible(return(result))
    
}
