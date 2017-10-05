#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(dplyr)
require(tibble)
require(tidyr)

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
