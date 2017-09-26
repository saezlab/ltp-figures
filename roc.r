#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(ggplot2)
require(dplyr)
require(lazyeval)
require(stringr)
require(ggrepel)

#
# constants
#

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_lit_ligands <- 'binding_properties_plain.csv'

#
# functions
#

preprocess_minimal <- function(){
    
    a  <- read.table(infile_a, sep = '\t', header = TRUE)
    e  <- read.table(infile_e, sep = '\t', header = TRUE)
    ae <- rbind(a, e) %>% mutate(lit = lit == 'True')
    
    return(ae)

}

preprocess_by_hg <- function(){
    
    d <- preprocess_minimal() %>%
        filter(cls %in% c('I', 'II'))
    
    d <- d %>%
        mutate(
            i90 = above_quantile(d, 0.90)$q,
            i75 = above_quantile(d, 0.75)$q,
            i50 = above_quantile(d, 0.50)$q,
            i25 = above_quantile(d, 0.25)$q
        )
    
    #return(d)
    
    d <- d %>%
        group_by(protein, ionm, screen, uhgroup) %>%
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

all_ligands <- function(d, litlig){
    
    alllig <- list()
    alllit <- unique(reduce(litlig, c))
    allneg <- unique(as.character((d %>% filter(ionm == 'neg'))$uhgroup))
    allpos <- unique(as.character((d %>% filter(ionm == 'pos'))$uhgroup))
    # this is the most unbiased version I think
    alllig$neg <- intersect(alllit, allneg)
    alllig$pos <- intersect(alllit, allpos)
    alllig$all <- intersect(alllit, union(allneg, allpos))
    
    return(alllig)
    
}

by_hg_df <- function(){
    
    d      <- preprocess_by_hg()
    litlig <- read_literature_ligands()
    alllig <- all_ligands(d, litlig)
    litlig <- litlig[sapply(litlig, function(x){as.logical(length(intersect(x, alllig$all)))})]
    allpro <- names(litlig)
    spacen <- expand.grid(protein = allpro, uhgroup = alllig$neg, ionm = c('neg'), screen = c('A', 'E'))
    spacep <- expand.grid(protein = allpro, uhgroup = alllig$pos, ionm = c('pos'), screen = c('A', 'E'))
    space  <- rbind(spacen, spacep)
    bycols <- c('protein', 'ionm', 'uhgroup', 'screen')
    alllit <- unnest(
        data.frame(
            protein = names(litlig),
            uhgroup = I(litlig)
        )
    ) %>%
    mutate(in_lit = TRUE)
    
    print(sprintf('Size of space is %i.', dim(space)[1]))
    
    d0 <- d %>%
        filter(uhgroup %in% alllig$all & protein %in% allpro) %>%
        select(protein, ionm, uhgroup, cls, screen, lit, detected, i90, i75, i50, i25) %>%
        right_join(space, by = bycols) %>%
        left_join(alllit, by = c('protein', 'uhgroup')) %>%
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

roc_conditions <- function(){
    
    in_screen <- function(d, sc, cl){
        
        return(
            d %>%
            filter(cls == cl & screen == sc) %>%
            select(protein, screen, uhgroup) %>%
            group_by(protein, screen, uhgroup) %>%
            summarise_all(first) %>%
            ungroup() %>%
            mutate(inn = TRUE)
        )
        
    }
    
    by_cols <- c('protein', 'screen', 'uhgroup')
    
    d <- preprocess_minimal() %>%
        filter(cls %in% c('I', 'II')) %>%
        group_by(protein, ionm, screen, id) %>%
        summarise_all(first) %>%
        ungroup()
    
    inAI  <- in_screen(d, 'A', 'I')  %>% rename(in_ai  = inn)
    inAII <- in_screen(d, 'A', 'II') %>% rename(in_aii = inn)
    inEI  <- in_screen(d, 'E', 'I')  %>% rename(in_ei  = inn)
    inEII <- in_screen(d, 'E', 'II') %>% rename(in_eii = inn)
    
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

roc_df <- function(){
    
    roc <- roc_conditions()
    
    rocdf <- do.call(rbind, lapply(roc, function(r){data.frame(c(list(name = r$name), r$vals))}))
    
    return(as.data.frame(rocdf))
    
}

roc_plot <- function(){
    
    rocdf <- roc_df()
    
    #rocdf <- rocdf[1:6,]
    #rocdf <- rocdf %>% filter(sens > 0.0)
    
    p <- ggplot(rocdf, aes(x = sens, y = 1 - spec, label = name)) +
        geom_abline(intercept = 0, slope = 1, color = 'red') +
        geom_point() +
        geom_text_repel(family = 'DINPro') +
        ggtitle('Various conditions in ROC space') +
        xlab('1 - specificity') +
        ylab('Sensitivity') +
        xlim(0.0, 1.0) +
        ylim(0.0, 1.0) +
        theme_bw() +
        theme(
            text = element_text(family = 'DINPro')
        )
    
    ggsave('roc_classes.pdf', device = cairo_pdf, width = 6, height = 6)
    
}
