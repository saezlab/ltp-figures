#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(dplyr)
require(tibble)
require(tidyr)

clustering.get_protein_ordr <- function(ae, return_cl = FALSE, method = 'ward.D2'){
    
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
    cl <- hclust(d, method = method)
    
    if(return_cl){
        
        return(cl)
        
    }else{
        
        return(cl$labels[cl$order])
        
    }
    
}

clustering.do_clustering <- function(d, x, y, v, method = 'ward.D2'){
    
    res <- list()
    xs <- deparse(substitute(x))
    x <- enquo(x)
    y <- enquo(y)
    v <- enquo(v)
    
    adj <- d %>%
        group_by(!!x, !!y) %>%
        mutate(mv = max(!!v)) %>%
        summarise_all(first) %>%
        select(!!x, !!y, mv) %>%
        spread(!!y, mv) %>%
        remove_rownames() %>%
        as.data.frame() %>%
        column_to_rownames(xs)
    
    adj[is.na(adj)] <- 0.0
    dix <- dist(adj)
    diy <- dist(t(adj))
    clx <- hclust(dix, method = method)
    cly <- hclust(diy, method = method)
    
    res$x  <- clx
    res$y  <- cly
    res$dx <- dix
    res$dy <- diy
    res$a  <- adj
    
    return(res)
    
}
