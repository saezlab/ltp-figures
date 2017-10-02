#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(readr)
require(dplyr)

infile_a   <- 'antonella_final.csv'
infile_e   <- 'enric_processed.csv'
infile_t   <- 'ltp_hptlc.tsv'

get_results <- function(){
    
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
    
    a <- suppressMessages(read_tsv(infile_a))
    e <- suppressMessages(read_tsv(infile_e))
    t <- suppressMessages(read_tsv(infile_t,
                                   col_names = c('protein', 'uhgroup')))
    
    apairs <- get_pairs(a)
    epairs <- get_pairs(e)
    
    a <- a %>%
        left_join(epairs, by = c('protein', 'uhgroup')) %>%
        filter((cls == 'I' | (cls == 'II' & in_other)) & uhgroup != 'P40')
    
    e <- e %>%
        left_join(apairs, by = c('protein', 'uhgroup')) %>%
        filter((cls == 'I' | (cls == 'II' & in_other)) & uhgroup != 'P40')
    
    result <- list()
    result$a <- a
    result$e <- e
    result$t <- t
    
    return(result)
    
}
