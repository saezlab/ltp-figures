#!/usr/bin/env Rscript

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

require(dplyr)
require(ggplot2)
require(grid)
require(gridExtra)
require(dendextend)
require(gtable)
require(ggdendro)

source('results.r')
source('clustering.r')

get_legend<-function(a.gplot){
    
    ## from
    ## http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
    leg <- tmp$grobs[[leg]]
    
    return(leg)
    
}

summary_heatmap <- function(wide = FALSE,
                            dendrogram = FALSE,
                            method = 'ward.D',
                            supervised = FALSE){
    
    dendrogram <- !supervised | dendrogram
    
    result <- list()
    
    pdfname <- sprintf(
        'ms-summary-heatmap%s%s%s.pdf',
        ifelse(wide, '-w', ''),
        ifelse(dendrogram, '-d', ifelse(supervised, '-s', '')),
        ifelse(dendrogram, sprintf('-%s', method), '')
    )
    width   <- ifelse(wide, 7.2, 5.7)
    height  <- 9
    main.width  <- .75
    main.height <- .80
    
    result$fname <- pdfname
    
    d  <- preprocess0(cluster_proteins = TRUE, method = method)
    #cl <- get_protein_ordr(d$c, return_cl = TRUE)
    
    if(wide){
        
        h <- d$c %>% group_by(protein, uhgroup)
        
        cl <- do_clustering(h, protein, uhgroup, lirel, method = method)
        
        h <- d$c %>%
            mutate(
                protein = factor(
                    protein,
                    levels = cl$x$labels[cl$x$order],
                    ordered = TRUE
                ),
                uhgroup = factor(
                    uhgroup,
                    levels = cl$y$labels[cl$y$order],
                    ordered = TRUE
                )
            ) %>%
            group_by(protein, uhgroup)
            
    }else{
        
        h <- d$c %>% group_by(protein, hg0)
        
        cl <- do_clustering(h, protein, hg0, lirel, method = method)
        
        h <- d$c %>%
            mutate(
                protein = factor(
                    protein,
                    levels = cl$x$labels[cl$x$order],
                    ordered = TRUE
                ),
                hg0 = factor(
                    hg0,
                    levels = cl$y$labels[cl$y$order],
                    ordered = TRUE
                )
            ) %>%
            group_by(protein, hg0)
        
    }
    
    h <- h %>%
        mutate(
            hg0_lit = any(lit0),
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
    
    if(supervised){
        
        h <- h %>%
            arrange(domain, protein) %>%
            mutate(
                protein = factor(protein, levels = unique(protein), ordered = TRUE)
            )
        
        if(wide){
            h <- h %>%
                mutate(grp = factor(grp, levels = grp_ordr, ordered = TRUE)) %>%
                arrange(grp, uhgroup) %>%
                mutate(
                    uhgroup = factor(uhgroup, levels = unique(uhgroup), ordered = TRUE)
                )
        }else{
            h <- h %>%
                mutate(grp = factor(grp, levels = grp_ordr, ordered = TRUE)) %>%
                arrange(grp, hg0) %>%
                mutate(
                    hg0 = factor(hg0, levels = unique(hg0), ordered = TRUE)
                )
        }
    }
    
    if(wide){
        p <- ggplot(h, aes(y = protein, x = uhgroup))
    }else{
        p <- ggplot(h, aes(y = protein, x = hg0))
    }
    
    p <- p +
        geom_tile(aes(fill = screens)) +
        geom_point(aes(alpha = hg0_lit), color = 'white') +
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
    
    if(!dendrogram){
        
        dcol <- domains_assign_colors(h)
        
        if(supervised){
            p <- p + theme(axis.text.y = element_text(color = domain_colors(h, dcol)))
            
            var <- ifelse(wide, 'uhgroup', 'hg0')
            
            lcol <- lipid_groups_assign_colors(h, var)
            p <- p +
                theme(
                    axis.text.x = element_text(
                        color = lipid_colors(h, var, lcol),
                        angle = 90, vjust = 0.5, size = 10, hjust = 1
                    )
                )
        }
        
        ggsave(pdfname, device = cairo_pdf, width = width, height = height)
        result <- list()
        result$data <- h
        result$heatmap <- p
        return(result)
    }
    
    # plotting dendrograms
    
#     dlp  <- as.dendrogram(cl$x) %>%
#         set('branches_lwd', c(.2, .2, .2)) %>%
#         as.ggdend() %>%
#         ggplot(labels = TRUE, horiz = TRUE) +
#         theme(
#             text = element_text(family = 'DINPro', size = 2),
#             plot.margin = unit(c(.3, .3, .3, .3), 'null')
#         )
#     dtp  <- dl  <- as.dendrogram(cl$x) %>%
#         set('branches_lwd', c(.2, .2, .2)) %>%
#         as.ggdend() %>%
#         ggplot(labels = TRUE) +
#         theme(text = element_text(family = 'DINPro'))
    
    ddatax <- dendro_data(cl$x, type = 'rectangle')
    dlp <- ggplot() +
        geom_segment(
            data = segment(ddatax),
            aes(x = x, y = y, xend = xend, yend = yend),
            lwd = .3
        ) +
        geom_text(
            data = label(ddatax),
            aes(x = x, y = y - 0.03, label = label, hjust = 0),
            family = 'DINPro-Medium',
            size = 3
        ) +
        coord_flip() +
        scale_y_reverse(
            expand = c(0, 0),
            limits = c(
                max(segment(ddatax)$y) * 1.2,
                min(segment(ddatax)$y) - 0.9
                
            )
        ) +
        scale_x_continuous(
            expand = c(0, .5)
        ) +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro'),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank()
        )
    
    ddatay <- dendro_data(cl$y, type = 'rectangle')
    dtp <- ggplot() +
        geom_segment(
            data = segment(ddatay),
            aes(x = x, y = y, xend = xend, yend = yend),
            lwd = .3
        ) +
        geom_text(
            data = label(ddatay),
            aes(x = x, y = y - 0.05, label = label, hjust = 1),
            family = 'DINPro-Medium',
            angle = 90,
            size = 3
        ) +
        scale_y_continuous(
            expand = c(0, 0),
            limits = c(
                min(segment(ddatay)$y) - 1.3,
                max(segment(ddatay)$y) * 1.2
            )
        ) +
        scale_x_continuous(
            expand = c(0, .5)
        ) +
        ggtitle(sprintf('Linkage method: %s', method)) +
        theme_minimal() +
        theme(
            text = element_text(family = 'DINPro'),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(size = 5)
        )
    
    cairo_pdf(pdfname, width = width, height = height)
        
        #grid.newpage()
        
#         top.layout <- grid.layout(
#             nrow = 2, ncol = 2,
#             widths  = unit(c(1 - main.width,  main.width ), 'null'),
#             heights = unit(c(1 - main.height, main.height), 'null')
#         )
        
        #pushViewport(viewport(layout = top.layout))
        
        p <- p +
            theme(
                legend.title = element_text(size = 11),
                legend.text  = element_text(size = 9),
                legend.spacing.y = unit(0.0, 'null')
            )
        
        leg <- get_legend(p)
        
        p <- p +
            theme(
                legend.position = 'none',
                plot.title = element_blank()
            ) +
            scale_y_discrete(position = 'right')
        
        dlp_gt <- ggplotGrob(dlp)
        dtp_gt <- ggplotGrob(dtp)
        p_gt   <- ggplotGrob(p)
        
        maxwidth <- unit.pmax(p_gt$widths, dtp_gt$widths)
        p_gt$widths   <- maxwidth
        dtp_gt$widths <- maxwidth
        
        maxheight <- unit.pmax(p_gt$heights, dlp_gt$heights)
        p_gt$heights   <- maxheight
        dlp_gt$heights <- maxheight
        
        grid.newpage()
        
        gm <- gtable_matrix(
            name = 'foobar',
            grobs = matrix(list(leg, dlp_gt, dtp_gt, p_gt),
            nrow = 2),
            widths = c(unit(c(0.25, 0.75), 'null')),
            heights = c(unit(c(0.2, 0.8), 'null'))
        )
        
        grid.draw(gm)
        
#         print(
#             p + theme(legend.position = 'none'),
#             vp = viewport(layout.pos.col = 2, layout.pos.row = 2)
#         )
#         
#         print(
#             dlp,
#             vp = viewport(layout.pos.col = 1, layout.pos.row = 2)
#         )
#         
#         print(
#             dtp,
#             vp = viewport(layout.pos.col = 2, layout.pos.row = 1)
#         )
        
#         leg <- get_legend(p)
#         pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
#         grid.draw(leg)
#         upViewport(0)
        
        # empty <- ggplot(data.frame(), aes(x = 1, y = 1)) + geom_blank()
        
    dev.off()
    
    result$data <- h
    result$grobs <- list()
    result$grobs$main <- p_gt
    result$grobs$left <- dlp_gt
    result$grobs$top  <- dtp_gt
    
    return(result)
    
}

test_linkage_methods <- function(){
    
    ws <- c(TRUE, FALSE)
    ms <- c('ward.D2', 'ward.D', 'single', 'complete', 'average', 'median', 'mcquitty', 'centroid')
    
    for(w in ws){
        
        for(m in ms){
            
            cat(sprintf('\tMethod: %s\n', m))
            
            nothing <- summary_heatmap(wide = w, dendrogram = TRUE, method = m)
            
        }
    }
    
}
