#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

source('results.r')
source('clustering.r')

# apply a multiplier in the lipid section
# in order to have enough space for labels
lmul <- 2

preprocess0 <- function(cluster_proteins = FALSE){
    
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
        
        protein_ordr <- get_protein_ordr(ae)
        
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

gen_circos_inputs <- function(cluster_proteins = FALSE){
    
    result <- list()
    
    d <- circos_preprocess(cluster_proteins = cluster_proteins)
    
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
            from = lag(cumsum(cnt_grp), 1) * lmul + max(p$to),
            to   = cumsum(cnt_grp) * lmul + max(p$to)
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
        arrange(grp, hg0, hgcc0) %>%
        mutate(
            grp = factor(grp, levels = grp_ordr, ordered = TRUE),
            hg0 = factor(hg0, levels = unique(hg0), ordered = TRUE),
            uhgroup = factor(uhgroup, levels = unique(uhgroup), ordered = TRUE),
            hgcc0 = factor(hgcc0, levels = unique(hgcc0), ordered = TRUE)
        ) %>%
        group_by(protein) %>%
        mutate(pro_offset = from +  1:n()) %>%
        ungroup() %>%
        group_by(grp) %>%
        mutate(grp_offset = from_grp + 1:n() * lmul) %>%
        ungroup() %>%
        mutate(
            pro_offset2 = pro_offset - 1,
            grp_offset2 = grp_offset - lmul
        ) %>%
        mutate(
            pro_offset = format(pro_offset, scientific = FALSE, justify = 'none', trim = TRUE),
            grp_offset = format(grp_offset, scientific = FALSE, justify = 'none', trim = TRUE),
            pro_offset2 = format(pro_offset2, scientific = FALSE, justify = 'none', trim = TRUE),
            grp_offset2 = format(grp_offset2, scientific = FALSE, justify = 'none', trim = TRUE)
        )
    
    li <- coo %>%
        mutate(linkcolor = paste0('color=', tolower(hg0))) %>%
        select(protein, pro_offset2, pro_offset, grp, grp_offset2, grp_offset, linkcolor)
    
    suppressMessages(li %>% write_delim('circos1/links.txt',
                                        col_names = FALSE))
    
    # generating individual labels for lipids
    
    la <- coo %>%
        select(grp, grp_offset2, grp_offset, hgcc0)
    
    suppressMessages(la %>% write_delim('circos1/labels.txt',
                                        col_names = FALSE))
    
    # generating subclassing level #1
    
    sc1 <- coo %>%
        group_by(hg0) %>%
        mutate(
            sc11 = format(
                min(as.integer(grp_offset2)),
                scientific = FALSE, justify = 'none', trim = TRUE),
            sc12 = format(
                max(as.integer(grp_offset)),
                scientific = FALSE, justify = 'none', trim = TRUE)
        ) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(sc1color = paste0('color=', tolower(hg0)))
    
    sc1tiles  <- sc1 %>% select(grp, sc11, sc12, sc1color)
    sc1labels <- sc1 %>% select(grp, sc11, sc12, hg0)
    
    suppressMessages(sc1tiles  %>% write_delim('circos1/subclass1.txt',
                                               col_names = FALSE))
    suppressMessages(sc1labels %>% write_delim('circos1/subclass1lab.txt',
                                               col_names = FALSE))
    
    # generating subclassing level #2
    
    sc2 <- coo %>%
        group_by(uhgroup) %>%
        mutate(
            sc21 = format(
                min(as.integer(grp_offset2)),
                scientific = FALSE, justify = 'none', trim = TRUE),
            sc22 = format(
                max(as.integer(grp_offset)),
                scientific = FALSE, justify = 'none', trim = TRUE)
        ) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(sc2color = paste0('color=', tolower(uhgroup)))
    
    sc2p <- coo %>%
        group_by(protein, uhgroup) %>%
        mutate(
            sc21 = format(
                min(as.integer(pro_offset2)),
                scientific = FALSE, justify = 'none', trim = TRUE),
            sc22 = format(
                max(as.integer(pro_offset)),
                scientific = FALSE, justify = 'none', trim = TRUE)
        ) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(sc2color = paste0('color=', tolower(uhgroup)))
    
    sc2tiles  <- bind_rows(
        sc2 %>% select(grp, sc21, sc22, sc2color),
        sc2p %>% select(grp = protein, sc21, sc22, sc2color)
    )
    sc2labels <- sc2 %>% select(grp, sc21, sc22, uhgroup)
    
    suppressMessages(sc2tiles  %>% write_delim('circos1/subclass2.txt',
                                                col_names = FALSE))
    suppressMessages(sc2labels %>% write_delim('circos1/subclass2lab.txt',
                                                col_names = FALSE))
    
    # generating intensity plots
    
    ity <- coo %>%
        left_join(
            coo %>%
            filter(ionm == 'neg') %>%
            group_by(protein, screen, hgcc0) %>%
            mutate(nlirel = max(lirel)) %>%
            summarise_all(first) %>%
            select(protein, screen, hgcc0, nlirel) %>%
            mutate(ionm = 'pos'),
            by = c('ionm', 'protein', 'screen', 'hgcc0')
        ) %>%
        mutate(
            nlirel = ifelse(
                ionm == 'neg',
                lirel,
                ifelse(is.na(nlirel), 0.0, nlirel)
            )
        ) %>%
        left_join(
            coo %>%
            filter(ionm == 'pos') %>%
            group_by(protein, screen, hgcc0) %>%
            mutate(plirel = max(lirel)) %>%
            summarise_all(first) %>%
            select(protein, screen, hgcc0, plirel) %>%
            mutate(ionm = 'neg'),
            by = c('ionm', 'protein', 'screen', 'hgcc0')
        ) %>%
        mutate(
            plirel = ifelse(
                ionm == 'pos',
                lirel,
                ifelse(is.na(plirel), 0.0, plirel)
            )
        ) %>%
        mutate(
            nlirel = format(-1 * nlirel, scientific = FALSE, justify = 'none', trim = TRUE),
            plirel = format(plirel, scientific = FALSE, justify = 'none', trim = TRUE)
        )
    
    nity <- bind_rows(
        ity %>%
        select(
            protein = grp,
            pro_offset2 = grp_offset2,
            pro_offset = grp_offset,
            nlirel
        ),
        ity %>% select(protein, pro_offset2, pro_offset, nlirel)
    )
    pity <-bind_rows(
        ity %>%
        select(
            protein = grp,
            pro_offset2 = grp_offset2,
            pro_offset = grp_offset,
            plirel
        ),
        ity %>% select(protein, pro_offset2, pro_offset, plirel)
    )
    
    suppressMessages(pity %>% write_delim('circos1/int-histo-pos.txt',
                                          col_names = FALSE))
    suppressMessages(nity %>% write_delim('circos1/int-histo-neg.txt',
                                          col_names = FALSE))
    
    # returning data frames
    
    result$chr <- all_chr
    result$li  <- li
    
    invisible(return(result))
    
}

get_protein_ordr <- function(ae, return_cl = FALSE){
    
    adj <- ae %>%
        group_by(protein, uhgroup) %>%
        mutate(mlirel = max(lirel)) %>%
        summarise_all(first) %>%
        select(protein, uhgroup, mlirel) %>%
        spread(uhgroup, mlirel) %>%
        remove_rownames() %>%
        as.data.frame() %>%
        column_to_rownames('protein')
    
    adj[is.na(adj)] <- 0.0
    d <- dist(adj)
    cl <- hclust(d, method = 'ward.D2')
    
    if(return_cl){
        
        return(cl)
        
    }else{
        
        return(cl$labels[cl$order])
        
    }
    
}

summary_heatmap <- function(wide = FALSE){
    
    pdfname <- ifelse(wide, 'ms-summary-heatmap-w.pdf', 'ms-summary-heatmap.pdf')
    width   <- ifelse(wide, 8, 6.5)
    
    d  <- circos_preprocess(cluster_proteins = TRUE)
    #cl <- get_protein_ordr(d$c, return_cl = TRUE)
    
    if(wide){
        h <- d$c %>%
            group_by(protein, uhgroup)
    }else{
        h <- d$c %>%
            group_by(protein, hg0)
    }
    
    h <- h %>%
        mutate(
            hg0_lit = any(lit),
            screens = paste0(sort(unique(screen)), collapse = '')
        ) %>%
        mutate(screens = factor(
            screens,
            levels = c('A', 'E', 'AE'),
            ordered = TRUE)
        ) %>%
        summarise_all(first) %>%
        ungroup()
        #mutate(protein = factor(protein, levels = cl, ordered = TRUE))
    
    if(wide){
        p <- ggplot(h, aes(y = protein, x = uhgroup))
    }else{
        p <- ggplot(h, aes(y = protein, x = hg0))
    }
    
    p <- p +
        geom_tile(aes(fill = screens)) +
        geom_point(aes(alpha = lit), color = 'white') +
        scale_alpha_manual(
            guide = guide_legend(title = 'Novelty'),
            values = c(
                'FALSE' = 1.0,
                'TRUE'  = 0.0
            ),
            labels = c(
                'FALSE' = 'Novelty',
                'TRUE'  = 'Known from\nliterature'
            )
        ) +
            scale_fill_manual(
            guide = guide_legend(title = 'Screens'),
            values = c(
                'A' = '#6EA945',
                'E' = '#FCCC06',
                'AE'= '#007B7F'
            ),
            labels = c(
                'A'  = 'HEK cells',
                'E'  = 'E. coli & liposomes',
                'AE' = 'Both'
            )
        ) +
        ggtitle('MS/MS screening of LTP cargoes: summary') +
        xlab('Lipids: main groups') +
        ylab('Proteins') +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro-Medium'),
            axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1),
            legend.key = element_rect(fill = '#007B7F', color = 'white')
        )
    
    ggsave(pdfname, device = cairo_pdf, width = width, height = 9)
    
    return(h)
    
}
