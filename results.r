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
infile_sa  <- 'stats-antonella.csv'
infile_se  <- 'stats-enric.csv'

max_deltart <- 0.5
neg_min_int <-  30000
pos_min_int <- 150000

supp_1_header <- c(
    
    "Protein",
    "Ionmode",
    "Feature ID",
    "m/z",
    "Corrected m/z",
    "Intensity",
    "Result level",
    "Manual identification",
    "Mean RT",
    "RT of closest MS2 scan",
    "Min RT",
    "Max RT",
    "Lyso",
    "Prefix",
    "Fatty acyl total carbon count",
    "Fatty acyl total unsaturation",
    "Fatty acyl 1 prefix",
    "Fatty acyl 1 carbon count",
    "Fatty acyl 1 unsaturation",
    "Fatty acyl 2 prefix",
    "Fatty acyl 2 carbon count",
    "Fatty acyl 2 unsaturation",
    "Fatty acyl 3 prefix",
    "Fatty acyl 3 carbon count",
    "Fatty acyl 3 unsaturation",
    "Headgroup",
    "Headgroup with total carbon count and unsaturation",
    "Total carbon count and unsaturation",
    "Headgroup with fatty acyls",
    "Fatty acyls",
    "Screen",
    "Method (MS or HPTLC)",
    "Lipid transfer domain family",
    "Headgroup category",
    "Major lipid category",
    "Same category in literature",
    "Present in other screen",
    "Relative intensity of headgroup category",
    "Same category present in screens",
    "Headgroup category with total carbon count and unsaturation",
    "Count of features in major lipid category",
    "Count of features in headgroup category",
    "Count of features at protein",
    "Count of features at domain family",
    "Count of features in headgroup category at domain family",
    "Count of features in headgroup category at protein",
    "Difference of mean RT and closest MS2 RT"
    
)


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

#' Here we define what we consider final results.
#' Finally we agreed class I and class II if cross
#' confirmed between screens.
#' We also remove low intensity class II features.
#' Called from pre_preprocess()
select_results <- function(d){
   
    (
        d %>%
        filter((cls == 'I' | (cls == 'II' & in_other)) & uhgroup != 'P40') %>%
        remove_lowintensity_classII()
    )
   
}

#' Recodes some annotations and adds literature data.
#' Finally filters to the final results.
pre_preprocess <- function(d, l, p, result_fun = select_results){
   
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
        mutate(grp = ifelse(
            hg0 %in% gpl, 'GPL', ifelse(
            hg0 %in% gl, 'GL', ifelse(
            hg0 %in% fa, 'FA', ifelse(
            hg0 %in% ch, 'CH', ifelse(
            hg0 %in% sph, 'SPH', ifelse(
            hg0 %in% vit, 'VIT', NA
        ))))))) %>%
        left_join(l, by = c('protein', 'hg0')) %>%
        mutate(lit0 = !is.na(lit0)) %>%
        left_join(p, by = c('protein', 'uhgroup')) %>%
        result_fun()
    )
   
}

#' Removes features which are not class I and are
#' below intensity thresholds.
remove_lowintensity_classII <- function(d){
   
    (
        d %>%
        filter(cls == 'I' | (
                ionm == 'pos' & intensity >= pos_min_int
            ) | (
                ionm == 'neg' & intensity >= neg_min_int
            )
        )
    )
   
}

#' Reads MS results and optionally adds HPTLC results.
get_results0 <- function(with_hptlc = TRUE){
   
    by_t <- c('protein', 'uhgroup', 'screen', 'method')
   
    a <- suppressMessages(read_tsv(infile_a)) %>% mutate(method = 'MS')
    e <- suppressMessages(read_tsv(infile_e)) %>% mutate(method = 'MS')
    do <- get_domains()
   
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
   
    a <- a %>%
        select(-domain) %>%
        left_join(do, by = c('protein'))
   
    e <- e %>%
        select(-domain) %>%
        left_join(do, by = c('protein'))
   
    result <- list()
    result$a <- a
    result$e <- e
    result$t <- t
   
    return(result)

}

#' Reads literature ligands and calls pre_preprocess().
#' At the end we get a tibbles with HPTLC included and
#' final results selected. The two screens are still in
#' separate tibbles and there is one with the HPTLC al-
#' though that is also included in the other 2 tibbles.
get_results <- function(with_hptlc = TRUE, result_fun = select_results){
   
    r <- get_results0(with_hptlc = with_hptlc)
    a <- r$a
    e <- r$e
   
    l <- literature_ligands()
    apairs <- get_pairs(a)
    epairs <- get_pairs(e)
   
    a <- pre_preprocess(a, l, epairs, result_fun = result_fun)
    e <- pre_preprocess(e, l, apairs, result_fun = result_fun)
   
    result <- list()
    result$a <- a
    result$e <- e
    if(with_hptlc){
        result$t <- r$t
    }
   
    return(result)
   
}


preprocess0 <- function(cluster_proteins = FALSE, method = 'ward.D2', with_hptlc = TRUE){
   
    result <- list()
    r <- get_results(with_hptlc = with_hptlc)
   
    # preprocessing
    # here we calculate relative intensities
    # and log relative intensities
    # relative means proportion in total intensities
    # of all identified features
    # HPTLC results are not yet added, this
    # is fine becaues those does not have intensities
    ae <- bind_rows(r$a, r$e) %>%
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
        ungroup()
   
    if(with_hptlc){
       
        # here we add the rows from the HPTLC dataset
        ae <- ae %>%
            full_join(r$t, by = c('protein', 'uhgroup', 'screen', 'method'))
       
    }
   
    # below we introduce hg0 which is a less detailed grouping of headgroups
    # and we count number of connections within various groups
    ae <- ae %>%
        group_by(protein, uhgroup) %>%
        mutate(screens = paste0(sort(unique(screen)), collapse = '')) %>%
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
        ungroup() %>%
        group_by(domain) %>%
        mutate(cnt_dom = n()) %>%
        ungroup() %>%
        group_by(uhgroup, domain) %>%
        mutate(cnt_hg_dom = n()) %>%
        ungroup() %>%
        group_by(protein, uhgroup) %>%
        mutate(cnt_hg_pro = n()) %>%
        ungroup() %>%
        mutate(deltart = rtmean - rtms2) %>%
        select(-mhgroup, -fullhgroup, -headgroup)
   
    # we either cluster proteins by their binders or
    # order by their domains and then alphabetically
    if(cluster_proteins){
       
        protein_ordr <- get_protein_ordr(ae, method = method)
       
    }else{
       
        protein_ordr <- unique(
            (
                ae %>%
                arrange(domain, protein)
            )$protein
        )
       
    }
   
    # building a lipid oriented data frame
    l <- ae %>%
        select(protein, domain, uhgroup, hg0, grp, hgcc0, screen,
               cnt_grp, cnt_hg, cnt_hg_dom, cnt_hg_pro) %>%
        group_by(hg0) %>%
        mutate(screens = paste0(sort(unique(screen)), collapse = '')) %>%
        ungroup() %>%
        group_by(hgcc0) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(grp = factor(grp, levels = grp_ordr, ordered = TRUE)) %>%
        # here we order by lipid classification and within uhgroups
        # we order the connections by domain
        arrange(grp, hg0, uhgroup, domain, protein) %>%
        mutate(
            hg0 = factor(hg0, levels = unique(hg0), ordered = TRUE),
            uhgroup = factor(uhgroup, levels = unique(uhgroup), ordered = TRUE)
        )
   
    # building a protein oriented data frame
    # ordering is defined by protein_ordr
    p <- ae %>%
        select(protein, domain, screen, cnt_pro, grp) %>%
        group_by(protein) %>%
        mutate(screens = paste0(sort(unique(screen)), collapse = '')) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(
            protein = factor(
                protein,
                levels = protein_ordr,
                ordered = TRUE
            )
        ) %>%
        arrange(protein)
   
    # full list of connections:
    result$c <- ae
    # all ligands:
    result$l <- l
    # all proteins:
    result$p <- p
    # proteins order:
    result$pordr <- protein_ordr
   
    return(result)
   
}


supp_table_1 <- function(){
    
    d <- preprocess0()$c %>%
        select(-lit, -ihg, -itotal, -lirel)
    
    write_tsv(as.data.frame(t(supp_1_header)),
              'SupplementaryTable2.tsv',
              col_names = FALSE)
    write_tsv(d, 'SupplementaryTable2.tsv', append = TRUE)
    
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


get_domains <- function(){
   
    m <- suppressMessages(read_tsv(infile_master1))
   
    return(
        m %>%
        rename(
            protein = default_name,
            domain = ltd_family
        ) %>%
        group_by(protein, domain) %>%
        summarise_all(first) %>%
        ungroup() %>%
        select(protein, domain) %>%
        mutate(
            domain = recode(
                domain,
                scp2 = 'Scp2',
                'NCP1_like' = 'NPC1-like',
                'LBP_BPI_CETP' = 'LBP/BPI/CETP',
                'IP_trans' = 'IP-trans',
                lipocalin = 'Lipocalin'
            )
        )
    )
   
}


domain_colors <- function(d, col){
   
    return(col[as.character(levels(factor(d$protein)))])
   
}

lipid_colors <- function(d, var, col){
   
    return(col[unique(as.character(d[[var]]))])
   
}

read_stats <- function(){
   
    coln <- c('protein', 'ionm', 'stat', 'val')
   
    sa <- suppressMessages(read_tsv(infile_sa, col_names = coln))
    se <- suppressMessages(read_tsv(infile_se, col_names = coln))
   
    return(list(sa = sa, se = se))
   
}


assign_colors <- function(d, var, grp_var, cols = NULL){
   
    groups <- levels(factor(d[[grp_var]]))
   
    if(is.null(cols)){
        cols <- c('#3A7AB3', '#608784', '#03928C', '#CF5836', '#7575BE',
                  '#D6708B', '#65B9B9', '#69B3D6', '#C441B3', '#9B998D')
    }
   
    cols <- cols[1:length(groups)]
   
    gc <- bind_cols(list(g = groups, c = cols))
   
    print(d[[var]])
    print(var)
   
    v <- as_tibble(list(var = levels(factor(d[[var]]))))
    names(v) <- var
   
    vg <- d %>%
        select_(grp_var, var) %>%
        group_by_(grp_var, var) %>%
        summarise_all(first)
   
   
    vc <- v %>%
        inner_join(vg, by = var)
   
    vc <- vc %>%
        inner_join(gc, by = setNames('g', grp_var)) %>%
        select_(var, 'c')
   
    ord_col <- vc$c
    names(ord_col) <- vc[[var]]
   
    return(ord_col)
   
    names(cols) <- sort(unique(d[[grp_var]]))
   
    grp_col <- list()
    for(i in 1:dim(d)[1]){
        grp_col[[as.character(d[[var]][i])]] <- as.character(d[[grp_var]][i])
    }
   
    ord_col <- as.character(cols[as.character(grp_col)])
    names(ord_col) <- names(grp_col)
   
    return(ord_col)
   
}

domains_assign_colors <- function(d){
   
    return(
        assign_colors(d, 'protein', 'domain')
    )
   
}

lipid_groups_assign_colors <- function(d, var = 'hg0'){
   
    return(
        assign_colors(d, var, 'grp')
    )
   
}
