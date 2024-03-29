#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(lazyeval)
require(stringr)
require(ggrepel)
require(purrr)
require(readr)
require(tidyr)

source('results.r')

#
# constants
#

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_lit_ligands <- 'binding_properties_plain.csv'
lossfile <- 'fall_out.csv'

#
# functions
#

roc_preprocess_minimal <- function(){
    
    #a  <- read.table(infile_a, sep = '\t', header = TRUE)
    #e  <- read.table(infile_e, sep = '\t', header = TRUE)
    r <- get_results(result_fun = function(x){return(x)})
    ae <- rbind(r$a, r$e)
    
    return(ae)

}

roc_preprocess_by_hg <- function(){
    
    d <- roc_preprocess_minimal() %>%
        filter(cls %in% c('I', 'II') & hg0 != 'P40') %>%
        droplevels()
    
    d <- d %>%
        mutate(
            i90 = above_quantile(d, 0.90)$q,
            i75 = above_quantile(d, 0.75)$q,
            i50 = above_quantile(d, 0.50)$q,
            i25 = above_quantile(d, 0.25)$q
        )
    
    d <- d %>%
        group_by(protein, ionm, screen, hg0) %>%
        summarise(
            lit = first(lit),
            cls = ifelse('I' %in% cls, 'I', 'II'),
            i90 = any(i90),
            i75 = any(i75),
            i50 = any(i50),
            i25 = any(i25)) %>%
        ungroup() %>%
        mutate(detected = TRUE)
    
    return(d)
    
}

read_literature_ligands <- function(){
    
    d <- read.csv(infile_lit_ligands, header = FALSE, sep = '\t')
    names(d) <- c('carrier', 'ligand')
    
    litlig <- sapply(levels(d$carrier),
                     function(cr){as.character(d$ligand[d$carrier == cr])})
    names(litlig) <- levels(d$carrier)
    
    return(litlig)
    
}

all_ligands <- function(d, litlig, full_space = TRUE){
    
    alllig <- list()
    alllit <- unique(litlig$hg0)
    allneg <- unique(as.character((d %>% filter(ionm == 'neg'))$hg0))
    allpos <- unique(as.character((d %>% filter(ionm == 'pos'))$hg0))
    if(full_space){
        alllig$all <- union(alllit, union(allneg, allpos))
    }else{
        # this is the most unbiased version I think
        alllig$neg <- intersect(alllit, allneg)
        alllig$pos <- intersect(alllit, allpos)
        alllig$all <- intersect(alllit, union(allneg, allpos))
    }
    
    return(alllig)
    
}

by_hg_df <- function(full_space = TRUE){
    
    allpro <- unique(roc_preprocess_minimal()$protein)
    d      <- roc_preprocess_by_hg()
    
    d <- d %>% filter(!is.na(hg0))
    
    litlig <- literature_ligands()
    alllig <- all_ligands(d, litlig, full_space = full_space)
    # litlig <- litlig[sapply(litlig, function(x){as.logical(length(intersect(x, alllig$all)))})]
    
    if(full_space){
        
        allpro <- unique(union(allpro, litlig$protein))
        space <- expand.grid(
            protein = unique(d$protein),
            hg0     = unique(union(d$hg0, litlig$hg0)),
            ionm    = c('neg', 'pos'),
            screen  = c('A', 'E')
        )
        
    }else{
        
        allpro <- names(litlig)
        spacen <- expand.grid(
            protein = allpro,
            hg0     = alllig$neg,
            ionm    = c('neg'),
            screen  = c('A', 'E')
        )
        spacep <- expand.grid(
            protein = allpro,
            hg0     = alllig$pos,
            ionm    = c('pos'),
            screen  = c('A', 'E')
        )
        space  <- rbind(spacen, spacep)
    }
    bycols <- c('protein', 'ionm', 'hg0', 'screen')
    alllit <- litlig %>% mutate(in_lit = TRUE)
    
    print(sprintf('Size of space is %i.', dim(space)[1]))
    
    d0 <- d %>%
        filter(hg0 %in% alllig$all & protein %in% allpro) %>%
        select(protein, ionm, hg0, cls, screen,
               lit, detected, i90, i75, i50, i25) %>%
        right_join(space, by = bycols) %>%
        left_join(alllit, by = c('protein', 'hg0')) %>%
        mutate(
            detected = !is.na(detected),
            lit = !is.na(lit | in_lit)
        )
    
    return(d0)
    
}

sens_spec <- function(pos, neg){
    
    tp <- dim(pos %>% filter( lit))[1]
    fp <- dim(pos %>% filter(!lit))[1]
    
    fn <- dim(neg %>% filter( lit))[1]
    tn <- dim(neg %>% filter(!lit))[1]
    
    return(list(
        sens = tp / (tp + fn),
        spec = tn / (tn + fp),
        prec = tp / (tp + fp),
        neg  = tn + fn,
        pos  = tp + fp,
        fpr  = fp / (fp + tn),
        fnr  = fn / (tp + fn),
        fdr  = fp / (tp + fp),
        tp   = tp,
        fp   = fp,
        tn   = tn,
        fn   = fn,
        n    = tp + fp + tn + fn
    ))
    
}

roc_point <- function(name, pos, neg){
    
    result <- list(
        name = name,
        vals = sens_spec(pos, neg)
    )
    
    result$pos <- pos
    result$neg <- neg
    
    result$lost <- neg %>%
        filter(detected) %>%
        select(protein, hg0) %>%
        left_join(
            pos %>%
                filter(detected) %>%
                mutate(in_pos = TRUE) %>%
                select(protein, hg0, in_pos),
            by = c('protein', 'hg0')
        ) %>%
        filter(is.na(in_pos)) %>%
        select(protein, hg0) %>%
        group_by(protein, hg0) %>%
        summarise_all(first) %>%
        ungroup() %>%
        arrange(protein, hg0)
    
    write(sprintf('\n### %s', name), lossfile, append = TRUE)
    write.table(result$lost, lossfile, sep = '\t', quote = FALSE, append = TRUE, col.names = FALSE)
    
    return(result)
    
}

above_quantile <- function(d, q, op = `>`){
    
    aq <- d %>%
        mutate(q = (
            (
                ionm == 'neg' &
                screen == 'A' &
                op(intensity, quantile((d %>% filter(ionm == 'neg' & screen == 'A'))$intensity, q))
            ) |
            (
                ionm == 'pos' &
                screen == 'A' &
                op(intensity, quantile((d %>% filter(ionm == 'pos' & screen == 'A'))$intensity, q))
            ) |
            (
                ionm == 'neg' &
                screen == 'E' &
                op(intensity, quantile((d %>% filter(ionm == 'neg' & screen == 'E'))$intensity, q))
            ) |
            (
                ionm == 'pos' &
                screen == 'E' &
                op(intensity, quantile((d %>% filter(ionm == 'pos' & screen == 'E'))$intensity, q))
            )
        ))
    
    return(aq)
    
}

in_screen <- function(d, condition){
    
    condition <- enquo(condition)
    
    return(
        d %>%
        filter(!!condition) %>%
        select(protein, hg0) %>%
        group_by(protein, hg0) %>%
        summarise_all(first) %>%
        ungroup() %>%
        mutate(inn = TRUE)
    )
    
}

roc_conditions_by_hg_grouped <- function(full_space = FALSE){
    
    by_cols <- c('protein', 'uhgroup')
    
    d <- by_hg_df(full_space = full_space)
    
    inAI  <- in_screen(d, screen == 'A' & cls == 'I' ) %>% rename(in_ai  = inn)
    inAII <- in_screen(d, screen == 'A' & cls == 'II') %>% rename(in_aii = inn)
    inEI  <- in_screen(d, screen == 'E' & cls == 'I' ) %>% rename(in_ei  = inn)
    inEII <- in_screen(d, screen == 'E' & cls == 'II') %>% rename(in_eii = inn)
    inPos <- in_screen(d, ionm == 'pos') %>% rename(in_pos = inn)
    inNeg <- in_screen(d, ionm == 'neg') %>% rename(in_neg = inn)
    
    d <- d %>%
        left_join(inAI,  by = by_cols) %>%
        left_join(inAII, by = by_cols) %>%
        left_join(inEI,  by = by_cols) %>%
        left_join(inEII, by = by_cols) %>%
        left_join(inPos, by = by_cols) %>%
        left_join(inNeg, by = by_cols) %>%
        mutate(
            in_ai  = !is.na(in_ai ),
            in_aii = !is.na(in_aii),
            in_ei  = !is.na(in_ei ),
            in_eii = !is.na(in_eii),
            in_pos = !is.na(in_pos),
            in_neg = !is.na(in_neg)
        ) %>%
        group_by(protein, uhgroup) %>%
        mutate(cls = ifelse('I' %in% cls, 'I', 'II')) %>%
        summarise_if( is.logical, any) %>%
        ungroup()
    
    return(list(
        list(
            name = 'Class I & II (all)',
            vals = sens_spec(
                d %>% filter( detected),
                d %>% filter(!detected)
            )
        ),
        list(
            name = 'Only positive',
            vals = sens_spec(
                d %>% filter( detected &  in_pos),
                d %>% filter(!detected | !in_pos)
            )
        ),
        list(
            name = 'Only negative',
            vals = sens_spec(
                d %>% filter( detected &  in_neg),
                d %>% filter(!detected | !in_neg)
            )
        ),
        list(
            name = 'Only 25% highest intensity',
            vals = sens_spec(
                d %>% filter( detected &  i75),
                d %>% filter(!detected | !i75)
            )
        ),
        list(
            name = 'Only 10% highest intensity',
            vals = sens_spec(
                d %>% filter( detected &  i90),
                d %>% filter(!detected | !i90)
            )
        ),
        list(
            name = 'Only 50% highest intensity',
            vals = sens_spec(
                d %>% filter( detected &  i50),
                d %>% filter(!detected | !i50)
            )
        ),
        list(
            name = 'Only 75% highest intensity',
            vals = sens_spec(
                d %>% filter( detected &  i25),
                d %>% filter(!detected | !i25)
            )
        ),
        list(
            name = 'Only class I',
            vals = sens_spec(
                d %>% filter( detected & ( in_ai |  in_ei)),
                d %>% filter(!detected | (!in_ai & !in_ei))
            )
        ),
        list(
            name = 'Only class II',
            vals = sens_spec(
                d %>% filter( detected & ( in_aii |  in_eii)),
                d %>% filter(!detected | (!in_aii & !in_eii))
            )
        ),
        list(
            name = 'In vivo, class II',
            vals = sens_spec(
                d %>% filter( detected &  in_aii),
                d %>% filter(!detected | !in_aii)
            )
        ),
        list(
            name = 'In vitro, class II',
            vals = sens_spec(
                d %>% filter( detected &  in_eii),
                d %>% filter(!detected | !in_eii)
            )
        ),
        list(
            name = 'In vivo, class I',
            vals = sens_spec(
                d %>% filter( detected &  in_ai),
                d %>% filter(!detected | !in_ai)
            )
        ),
        list(
            name = 'In vitro, class I',
            vals = sens_spec(
                d %>% filter( detected &  in_ei),
                d %>% filter(!detected | !in_ei)
            )
        ),
        list(
            name = 'In vitro, class I & II',
            vals = sens_spec(
                d %>% filter( detected & ( in_ei |  in_eii)),
                d %>% filter(!detected | (!in_ei & !in_eii))
            )
        ),
        list(
            name = 'In vivo, class I & II',
            vals = sens_spec(
                d %>% filter( detected & ( in_ai |  in_aii)),
                d %>% filter(!detected | (!in_ai & !in_aii))
            )
        ),
        list(
            name = 'Class I and in vivo class II',
            vals = sens_spec(
                d %>% filter(             in_ai |  in_ei |  in_aii),
                d %>% filter(!detected & !in_ai & !in_ei & !in_aii)
            )
        ),
        list(
            name = 'Positive in vitro class I and in vivo class II',
            vals = sens_spec(
                d %>% filter(
                    detected &
                    in_pos & (
                        in_ei | in_aii
                    )
                ),
                d %>% filter(
                    !detected |
                    !in_pos | (
                        !in_ei & !in_aii
                    )
                )
            )
        ),
        list(
            name = 'Class I & in vivo class II confirmed by in vitro',
            vals = sens_spec(
                d %>% filter( detected &  (in_ai | in_ei | in_aii & (in_ei | in_eii))),
                d %>% filter(!detected | !(in_ai | in_ei | in_aii & (in_ei | in_eii)))
            )
        ),
        list(
            name = '50% highest intensity: class I & in vivo class II confirmed by in vitro',
            vals = sens_spec(
                d %>% filter( detected &  i50 &  (in_ai | in_ei | in_aii & (in_ei | in_eii))),
                d %>% filter(!detected | !i50 | !(in_ai | in_ei | in_aii & (in_ei | in_eii)))
            )
        ),
        list(
            name = 'Class I & class II confirmed by other screen',
            vals = sens_spec(
                d %>% filter(
                    detected & (
                        in_ai | in_ei |
                        (in_aii & (in_ei | in_eii)) |
                        (in_eii & (in_ai | in_aii))
                    )
                ),
                d %>% filter(
                    !detected | !(
                        in_ai | in_ei |
                        (in_aii & (in_ei | in_eii)) |
                        (in_eii & (in_ai | in_aii))
                    )
                )
            )
        ),
        list(
            name = '50% highest intensity: class I & class II confirmed by other screen',
            vals = sens_spec(
                d %>% filter(
                    detected &
                    i50 & (
                        in_ai | in_ei |
                        (in_aii & (in_ei | in_eii)) |
                        (in_eii & (in_ai | in_aii))
                    )
                ),
                d %>% filter(
                    !detected |
                    !i50 | !(
                        in_ai | in_ei |
                        (in_aii & (in_ei | in_eii)) |
                        (in_eii & (in_ai | in_aii))
                    )
                )
            )
        )
    ))
    
}

roc_conditions_by_hg <- function(group_sets = FALSE, full_space = FALSE, verbose = FALSE){
    
    write('', file = lossfile)
    
    by_cols <- c('protein', 'hg0')
    
    d <- by_hg_df(full_space = full_space)
    
    inAI  <- in_screen(d, screen == 'A' & cls == 'I' ) %>% rename(in_ai  = inn)
    inAII <- in_screen(d, screen == 'A' & cls == 'II') %>% rename(in_aii = inn)
    inEI  <- in_screen(d, screen == 'E' & cls == 'I' ) %>% rename(in_ei  = inn)
    inEII <- in_screen(d, screen == 'E' & cls == 'II') %>% rename(in_eii = inn)
    
    d <- d %>%
        left_join(inAI,  by = by_cols) %>%
        left_join(inAII, by = by_cols) %>%
        left_join(inEI,  by = by_cols) %>%
        left_join(inEII, by = by_cols) %>%
        mutate(
            in_ai  = !is.na(in_ai ),
            in_aii = !is.na(in_aii),
            in_ei  = !is.na(in_ei ),
            in_eii = !is.na(in_eii)
        )
    
    if(group_sets){
        
        d_any     <- d %>%
            group_by(protein, hg0) %>%
            mutate(cls = ifelse('I' %in% cls, 'I', 'II')) %>%
            summarise_if( is.logical, any) %>%
            summarise_if(is.character, first) %>%
            ungroup()
        
        d_ionm    <- d %>%
            group_by(protein, hg0, ionm) %>%
            summarise_if(is.logical, any) %>%
            ungroup()
        
        d_cls     <- d %>%
            group_by(protein, hg0) %>%
            summarise_if(is.logical, any) %>%
            ungroup()
        
        d_scr     <- d %>%
            group_by(protein, hg0, screen) %>%
            summarise_if(is.logical, any) %>%
            ungroup()
        
        d_scr_cls <- d %>%
            group_by(protein, hg0, screen, cls) %>%
            summarise_if(is.logical, any) %>%
            ungroup()
    
    }else{
        
        d_any     <- d
        d_ionm    <- d
        d_cls     <- d
        d_scr     <- d
        d_scr_cls <- d
        
    }
    
    return(list(
        roc_point(
            name = 'Class I | II',
            pos = d_any %>% filter( detected),
            neg = d_any %>% filter(!detected)
        ),
        roc_point(
            name = 'Positive',
            pos = d_ionm %>% filter( detected & ionm == 'pos'),
            neg = d_ionm %>% filter(!detected | ionm != 'pos')
        ),
        roc_point(
            name = 'Negative',
            pos = d_ionm %>% filter( detected & ionm == 'neg'),
            neg = d_ionm %>% filter(!detected | ionm != 'neg')
        ),
        roc_point(
            name = '25% highest intensity',
            pos = d_any %>% filter( detected &  i75),
            neg = d_any %>% filter(!detected | !i75)
        ),
        roc_point(
            name = '10% highest intensity',
            pos = d_any %>% filter( detected &  i90),
            neg = d_any %>% filter(!detected | !i90)
        ),
        roc_point(
            name = '50% highest intensity',
            pos = d_any %>% filter( detected &  i50),
            neg = d_any %>% filter(!detected | !i50)
        ),
        roc_point(
            name = '75% highest intensity',
            pos = d_any %>% filter( detected &  i25),
            neg = d_any %>% filter(!detected | !i25)
        ),
        roc_point(
            name = 'Class I',
            pos = d_cls %>% filter( detected & cls == 'I'),
            neg = d_cls %>% filter(!detected | cls != 'I')
        ),
        roc_point(
            name = 'Class II',
            pos = d_cls %>% filter( detected & cls == 'II'),
            neg = d_cls %>% filter(!detected | cls != 'II')
        ),
        roc_point(
            name = 'In vivo & class II',
            pos = d_scr_cls %>% filter(  detected & cls == 'II' & screen == 'A'),
            neg = d_scr_cls %>% filter(!(detected & cls == 'II' & screen == 'A'))
        ),
        roc_point(
            name = 'In vitro & class II',
            pos = d_scr_cls %>% filter( detected & cls == 'II' & screen == 'E'),
            neg = d %>% filter(!detected | cls != 'II' | screen != 'E')
        ),
        roc_point(
            name = 'In vivo & class I',
            pos = d_scr_cls %>% filter( detected & cls == 'I' & screen == 'A'),
            neg = d_scr_cls %>% filter(!detected | cls != 'I' | screen != 'A')
        ),
        roc_point(
            name = 'In vitro & class I',
            pos = d_scr_cls %>% filter( detected &  cls == 'I' & screen == 'E'),
            neg = d_scr_cls %>% filter(!detected |  cls != 'I' | screen != 'E')
        ),
        roc_point(
            name = 'In vitro & (class I | II)',
            pos = d_scr_cls %>% filter( detected & screen == 'E'),
            neg = d_scr_cls %>% filter(!detected | screen != 'E')
        ),
        roc_point(
            name = 'In vivo & (class I | II)',
            pos = d_scr_cls %>% filter( detected & screen == 'A'),
            neg = d_scr_cls %>% filter(!detected | screen != 'A')
        ),
        roc_point(
            name = 'Class I | (in vivo & class II)',
            pos = d_scr_cls %>% filter( cls == 'I' | (screen == 'A' & cls =='II')),
            neg = d_scr_cls %>% filter(!detected | (cls == 'II' & screen != 'A'))
        ),
        roc_point(
            name = 'Positive & ((in vitro & class I) | (in vivo & class II))',
            pos = d %>% filter(
                detected &
                ionm == 'pos' & (
                    (screen == 'E' & cls == 'I')  |
                    (screen == 'A' & cls =='II')
                )
            ),
            neg = d %>% filter(
                !detected |
                ionm == 'neg' | (
                    (screen == 'E' & cls == 'II') |
                    (screen == 'A' & cls =='I' )
                )
            )
        ),
        roc_point(
            name = 'Class I | (in vivo & class II & confirmed by in vitro)',
            pos = d_scr_cls %>% filter(detected & (cls == 'I'  | (screen == 'A' & (in_ei | in_eii)))),
            neg = d_scr_cls %>% filter(!detected | (cls == 'II' & (screen != 'A' | (!in_ei & !in_eii))))
        ),
        roc_point(
            name = '50% highest intensity & (class I | (in vivo & class II & confirmed by in vitro))',
            pos = d_scr_cls %>% filter( detected &  i50 & (cls == 'I'  | (screen == 'A' & (in_ei | in_eii)))),
            neg = d_scr_cls %>% filter(!detected | !i50 | (cls == 'II' & (screen != 'A' | (!in_ei & !in_eii))))
        ),
        roc_point(
            name = 'Class I | (class II & confirmed by other screen)',
            pos = d_scr_cls %>% filter(
                detected & (
                    cls == 'I'  |
                    (screen == 'A' & (in_ei | in_eii)) |
                    (screen == 'E' & (in_ai | in_aii))
                )
            ),
            neg = d_scr_cls %>% filter(
                !detected | (
                    cls != 'I'  & (
                        (screen == 'A' & (!in_ei & !in_eii)) |
                        (screen == 'E' & (!in_ai & !in_aii))
                    )
                )
            )
        ),
        roc_point(
            name = '50% highest intensity & (class I | (class II & confirmed by other screen))',
            pos = d_scr_cls %>% filter(
                detected &
                i50 & (
                    cls == 'I'  |
                    (screen == 'A' & (in_ei | in_eii)) |
                    (screen == 'E' & (in_ai | in_aii))
                )
            ),
            neg = d_scr_cls %>% filter(
                !detected |
                !i50 | (
                    cls != 'I'  & (
                        (screen == 'A' & (!in_ei & !in_eii)) |
                        (screen == 'E' & (!in_ai & !in_aii))
                    )
                )
            )
        )
    ))
    
}

roc_conditions <- function(){
    
    by_cols <- c('protein', 'uhgroup')
    
    d <- roc_preprocess_minimal() %>%
        filter(cls %in% c('I', 'II')) %>%
        group_by(protein, ionm, screen, id) %>%
        summarise_all(first) %>%
        ungroup()
    
    inAI  <- in_screen(d, screen == 'A' & cls == 'I' ) %>% rename(in_ai  = inn)
    inAII <- in_screen(d, screen == 'A' & cls == 'II') %>% rename(in_aii = inn)
    inEI  <- in_screen(d, screen == 'E' & cls == 'I' ) %>% rename(in_ei  = inn)
    inEII <- in_screen(d, screen == 'E' & cls == 'II') %>% rename(in_eii = inn)
    
    d <- d %>%
        left_join(inAI,  by = by_cols) %>%
        left_join(inAII, by = by_cols) %>%
        left_join(inEI,  by = by_cols) %>%
        left_join(inEII, by = by_cols) %>%
        mutate(
            in_ai  = !is.na(in_ai ),
            in_aii = !is.na(in_aii),
            in_ei  = !is.na(in_ei ),
            in_eii = !is.na(in_eii)
        )
    
    custom1pos <- d %>%
        filter(
            ionm == 'pos' & (
                (screen == 'E' & cls == 'I')  |
                (screen == 'A' & cls =='II')
            )
        ) %>% above_quantile(0.5, `>`) %>% filter(q)
    
    custom1neg <- d %>%
        left_join(custom1pos,
                  by = c('protein', 'ionm', 'screen', 'uhgroup', 'id'),
                  suffix = c('', 'y')
        ) %>%
        filter(is.na(q))
    
    return(list(
        list(
            name = 'Only positive',
            vals = sens_spec(
                d %>% filter(ionm == 'pos'),
                d %>% filter(ionm != 'pos')
            )
        ),
        list(
            name = 'Only negative',
            vals = sens_spec(
                d %>% filter(ionm == 'neg'),
                d %>% filter(ionm != 'neg')
            )
        ),
        list(
            name = 'Only 25% highest intensity',
            vals = sens_spec(
                above_quantile(d, 0.75, `>`) %>% filter(q),
                above_quantile(d, 0.75, `<=`) %>% filter(q)
            )
        ),
        list(
            name = 'Only 10% highest intensity',
            vals = sens_spec(
                above_quantile(d, 0.90, `>`) %>% filter(q),
                above_quantile(d, 0.90, `<=`) %>% filter(q)
            )
        ),
        list(
            name = 'Only 50% highest intensity',
            vals = sens_spec(
                above_quantile(d, 0.50, `>`) %>% filter(q),
                above_quantile(d, 0.50, `<=`) %>% filter(q)
            )
        ),
        list(
            name = 'Only 75% highest intensity',
            vals = sens_spec(
                above_quantile(d, 0.25, `>`) %>% filter(q),
                above_quantile(d, 0.25, `<=`) %>% filter(q)
            )
        ),
        list(
            name = 'Only class I',
            vals = sens_spec(
                d %>% filter(cls == 'I'),
                d %>% filter(cls != 'I')
            )
        ),
        list(
            name = 'Only class II',
            vals = sens_spec(
                d %>% filter(cls == 'II'),
                d %>% filter(cls != 'II')
            )
        ),
        list(
            name = 'In vivo, class II',
            vals = sens_spec(
                d %>% filter(cls == 'II' & screen == 'A'),
                d %>% filter(cls != 'II' | screen != 'A')
            )
        ),
        list(
            name = 'In vitro, class II',
            vals = sens_spec(
                d %>% filter(cls == 'II' & screen == 'E'),
                d %>% filter(cls != 'II' | screen != 'E')
            )
        ),
        list(
            name = 'In vivo, class I',
            vals = sens_spec(
                d %>% filter(cls == 'I' & screen == 'A'),
                d %>% filter(cls != 'I' | screen != 'A')
            )
        ),
        list(
            name = 'In vitro, class I',
            vals = sens_spec(
                d %>% filter(cls == 'I' & screen == 'E'),
                d %>% filter(cls != 'I' | screen != 'E')
            )
        ),
        list(
            name = 'In vitro, class I & II',
            vals = sens_spec(
                d %>% filter(screen == 'E'),
                d %>% filter(screen != 'E')
            )
        ),
        list(
            name = 'In vivo, class I & II',
            vals = sens_spec(
                d %>% filter(screen == 'A'),
                d %>% filter(screen != 'A')
            )
        ),
        list(
            name = 'Class I and in vivo class II',
            vals = sens_spec(
                d %>% filter(cls == 'I' | (screen == 'A' & cls =='II')),
                d %>% filter(cls == 'II' & screen != 'A')
            )
        ),
        list(
            name = 'Positive in vitro class I and in vivo class II',
            vals = sens_spec(
                d %>% filter(
                    ionm == 'pos' & (
                        (screen == 'E' & cls == 'I')  |
                        (screen == 'A' & cls =='II')
                    )
                ),
                d %>% filter(
                    ionm == 'neg' | (
                        (screen == 'E' & cls == 'II') |
                        (screen == 'A' & cls =='I' )
                    )
                )
            )
        ),
        list(
            name = 'Class I & in vivo class II confirmed by in vitro',
            vals = sens_spec(
                d %>% filter(cls == 'I'  | (screen == 'A' & (in_ei | in_eii))),
                d %>% filter(cls == 'II' & (screen != 'A' | (!in_ei & !in_eii)))
            )
        ),
        list(
            name = 'Class I & class II confirmed by other screen',
            vals = sens_spec(
                d %>% filter(
                    cls == 'I'  |
                    (screen == 'A' & (in_ei | in_eii)) |
                    (screen == 'E' & (in_ai | in_aii))
                ),
                d %>% filter(
                    cls != 'I'  & (
                        (screen == 'A' & (!in_ei & !in_eii)) |
                        (screen == 'E' & (!in_ai & !in_aii))
                    )
                )
            )
        )
    ))
    
}


roc_df <- function(by_hg = FALSE, group_sets = FALSE, full_space = FALSE, verbose = FALSE){
    
    if(group_sets){
        
        roc <- roc_conditions_by_hg_grouped(full_space = full_space, verbose = verbose)
        
    }else if(by_hg){
        
        roc <- roc_conditions_by_hg(full_space = full_space, verbose = verbose)
        
    }else{
        
        roc <- roc_conditions(verbose = verbose)
        
    }
    
    rocdf <- do.call(rbind, lapply(roc, function(r){data.frame(c(list(name = r$name), r$vals))}))
    
    return(as.data.frame(rocdf))
    
}


roc_plot <- function(by_hg = FALSE,
                     group_sets = FALSE,
                     full_space = FALSE,
                     lim = 1.0,
                     title = 'ROC over all identified features',
                     verbose = FALSE){
    
    rocdf <- roc_df(
        by_hg = by_hg,
        group_sets = group_sets,
        full_space = full_space,
        verbose = verbose
    )
    fname <- ifelse(by_hg,
        ifelse(full_space,
            'roc_classes_full',
            'roc_classes_by-hg'),
        'roc_classes'
    )
    
    #rocdf <- rocdf[1:6,]
    #rocdf <- rocdf %>% filter(sens > 0.0)
    
    p <- ggplot(rocdf, aes(y = sens, x = 1 - spec, label = name)) +
        geom_abline(intercept = 0, slope = 1, color = 'red') +
        geom_point(shape = 1) +
        geom_label_repel(
            family = 'DINPro',
            size = 3,
            max.iter = 40000,
            fill = 'black',
            color = 'white',
            segment.color = 'grey30',
            segment.size = 0.2,
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.35, "lines")
        ) +
        ggtitle(title) +
        xlab('1 - specificity') +
        ylab('Sensitivity') +
        xlim(0.0, lim) +
        ylim(0.0, lim) +
        theme_linedraw() +
        theme(
            text = element_text(family = 'DINPro'),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 21),
            plot.title = element_text(size = 24)
        )
    
    pdfname <- sprintf('%s.pdf', fname)
    
    cat(sprintf('Plotting into `%s`.', pdfname))
    
    ggsave(pdfname, device = cairo_pdf, width = 6, height = 6)
    rocdf %>% write_tsv(sprintf('%s.tsv', fname))
    
}

