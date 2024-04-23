get_cpdb_df = function(cpdb_dir, celltypes){
    #read in the files
    cpdb_files_in = grep(".txt", list.files(cpdb_dir), value = T) 
    cpdb_in = list()
    for(i in cpdb_files_in){
        tab =  read.delim(file.path(cpdb_dir, i), check.names = FALSE)
        ispl = unlist(strsplit(i, "_"))
        i = paste0(ispl[[1]], "_", ispl[[2]], "_", ispl[[3]])    
        cpdb_in[[i]] = tab
    }
    #get a unified interaction index
    cpdb_index <- data.frame(cpdb_in$statistical_analysis_means[, 1:11],row.names = cpdb_in$statistical_analysis_means[, 1])
    #make a long interaction data frame
    means <- reshape::melt(data.frame(cpdb_in$statistical_analysis_means[, 12:ncol(cpdb_in$statistical_analysis_means)], row.names = cpdb_index$id_cp_interaction, "interaction" = cpdb_index$id_cp_interaction))
    p_vals <- reshape::melt(data.frame(cpdb_in$statistical_analysis_pvalues[, 12:ncol(cpdb_in$statistical_analysis_pvalues)], row.names = cpdb_index$id_cp_interaction, "interaction" = cpdb_index$id_cp_interaction))
    #make an overall interaction dataframe
    rownames(means) = paste0(means$interaction, "_", means$variable)
    rownames(p_vals) = paste0(p_vals$interaction, "_", p_vals$variable)
    interaction_df = data.frame('interaction' = means$interaction, 'cells' = means$variable,
                                'mean' = means$value, 'pval' = p_vals$value)

    #get the possible combinations of cell types
        combinations = list()
        p=1
        for(i in celltypes){
            for(j in celltypes){
                i = gsub("-", ".", i)
                j = gsub("-", ".", j)
            combinations[[p]] = data.frame("ct1" = i,  "ct2" = j)
                    p = p+1
        }}
    combinations = do.call(rbind, combinations)
    combinations$combination = paste0(combinations$ct1, ".", combinations$ct2)
    combinations$ct1 = gsub("[.]", "-", combinations$ct1 )
    combinations$ct2 = gsub("[.]", "-", combinations$ct2 )

    #remove entirely insignificant interactions 
    message("removing entirely insignificant interactions")
    agp = aggregate(interaction_df$pval, by = list('interaction' = interaction_df$interaction), sum)
    psum = agp$x
    names(psum) = agp$interaction
    psum = psum == unique(table(interaction_df$interaction))
    sig_interactions = names(psum)[!psum]
    interaction_df = interaction_df[interaction_df$interaction %in% sig_interactions , ]
    
        message("adding cell-cell info")
    #write in the combinations to interaction df
    idx = match(interaction_df$cells, combinations$combination)
    interaction_df[, 'sender'] = combinations[idx, 1]
    interaction_df[, 'receiver'] = combinations[idx, 2]
    interaction_df$cells = paste0(interaction_df$sender, "|", interaction_df$receiver)   

            #add the interactions 
        message("adding interactions/molecules")
        interactions <- cpdb_index[interaction_df$interaction, ]
        #then get the molecules
        interaction_df$molecule_a <- unlist(lapply(strsplit(as.character(interactions$interacting_pair), "[_]"), function(x){x[1]}))
        interaction_df$molecule_b <- unlist(lapply(strsplit(as.character(interactions$interacting_pair), "[_]"), function(x){x[2]}))
        cpdb_df = cbind(interaction_df, interactions)

        #now we want to get the ligands and receptors the right way round because inexplicably cpdb get's this all askew
      cpdb_df$receptor_cell <- NA
      cpdb_df$receptor_cell[cpdb_df$receptor_a %in% "True"] <- cpdb_df$sender[cpdb_df$receptor_a %in% "True"]
      cpdb_df$receptor_cell[cpdb_df$receptor_b %in% "True"] <- cpdb_df$receiver[cpdb_df$receptor_b %in% "True"]
      cpdb_df$ligand_cell <- NA
      cpdb_df$ligand_cell[cpdb_df$receptor_a %in% "False"] <- cpdb_df$sender[cpdb_df$receptor_a %in% "False"]
      cpdb_df$ligand_cell[cpdb_df$receptor_b %in% "False"] <- cpdb_df$receiver[cpdb_df$receptor_b %in% "False"]
      cpdb_df$receptor <- NA
      cpdb_df$receptor[cpdb_df$receptor_a %in% "True"] <- cpdb_df$molecule_a[cpdb_df$receptor_a %in% "True"]
      cpdb_df$receptor[cpdb_df$receptor_b %in% "True"] <- cpdb_df$molecule_b[cpdb_df$receptor_b %in% "True"]
      cpdb_df$ligand <- NA
      cpdb_df$ligand[cpdb_df$receptor_a %in% "False"] <- cpdb_df$molecule_a[cpdb_df$receptor_a %in% "False"]
      cpdb_df$ligand[cpdb_df$receptor_b %in% "False"] <- cpdb_df$molecule_b[cpdb_df$receptor_b %in% "False"]

      cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "ligand_cell"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "sender"]
      cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receptor_cell"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receiver"]
      cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "ligand"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "molecule_a"]
      cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receptor"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "molecule_b"]

      cpdb_df$LR <- paste0(cpdb_df$ligand, "->", cpdb_df$receptor)
      cpdb_df$SR <- paste0(cpdb_df$ligand_cell, "->", cpdb_df$receptor_cell)


      return(cpdb_df)  
}


#cellphone DB plot
min_max_normalisation <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
}

cellphone_db_plot <- function(cpdb_df, celltypes, senders, receivers, cluster_interactions = TRUE, 
                              filter_secreted  = FALSE, scaling = c("min_max", "scale"), remove_auto = TRUE,
                              filter_integrins = TRUE, 
                             gene_family = NULL, custom_genes = NULL){
    
      #get these ordered by celltype levels
  #egrid <- expand.grid(receivers, senders)
  #cpdb_df$SR <- factor(cpdb_df$SR, levels =  paste0(egrid$Var2, "->", egrid$Var1))
    egrid <- expand.grid(celltypes, celltypes)
   cpdb_df$SR <- factor(cpdb_df$SR, levels =  paste0(egrid$Var2, "->", egrid$Var1))

  #subset to just secreted interactions
  if(filter_secreted){
    cpdb_df <- cpdb_df[cpdb_df$secreted %in% "True", ]
  }
  
  #remove auto interactions
  if(remove_auto){
    cpdb_df <- cpdb_df[!cpdb_df$sender == cpdb_df$receiver, ]
    cpdb_df <- cpdb_df[!cpdb_df$molecule_a == cpdb_df$molecule_b, ]
  }
  
  if(filter_integrins){
    cpdb_df <- cpdb_df[cpdb_df$is_integrin %in% "False", ]
  }
  
  #subset to the cells we are interested in 
  cell_types <- union(senders, receivers)
  cpdb_df <- cpdb_df[cpdb_df$sender %in% cell_types + cpdb_df$receiver %in% cell_types == 2, ]
  
  #subset to senders and receivers
  cpdb_df <- cpdb_df[cpdb_df$ligand_cell %in% senders, ]
  cpdb_df <- cpdb_df[cpdb_df$receptor_cell %in% receivers, ]
  
#get the gene family we want to plot
if(is.null(gene_family)){
    message("no gene family selected")
}else{
    if(!is.null(custom_genes)){
           gene_family_list = list(
         'chemokines' =  grep("^CXC|^CCL|^CCR|^CX3|XCL|XCR|^IL|^AREG|^TGF|^CSF", cpdb_df$interacting_pair),
            'th1'=  grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4",
                cpdb_df$interacting_pair),
            'th2' = grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", cpdb_df$interacting_pair),
            'th17' = grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB",
                cpdb_df$interacting_pair),
            'treg' = grep("IL35|IL10|FOXP3|IL2RA|TGFB", cpdb_df$interacting_pair),
            'costimulatory' = grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]",
                cpdb_df$interacting_pair),
            'coinhibitory' = grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR",
                cpdb_df$interacting_pair),
            'niche' = grep("CSF", cpdb_df$interacting_pair),
               'custom_genes' =  grep(paste(custom_genes,collapse="|"), 
                        cpdb_df$interacting_pair)
            )

    interactions_sub = unique(unlist(gene_family_list[c(gene_family, 'custom_genes')]))  
    cpdb_df = cpdb_df[interactions_sub, ]
 
    }else{
         gene_family_list = list(
         'chemokines' =  grep("^CXC|^CCL|^CCR|^CX3|XCL|XCR|^IL|^AREG|^TGF|^CSF", cpdb_df$interacting_pair),
            'th1'=  grep("IL2|IL12|IL18|IL27|IFNG|IL10|TNF$|TNF |LTA|LTB|STAT1|CCR5|CXCR3|IL12RB1|IFNGR1|TBX21|STAT4",
                cpdb_df$interacting_pair),
            'th2' = grep("IL4|IL5|IL25|IL10|IL13|AREG|STAT6|GATA3|IL4R", cpdb_df$interacting_pair),
            'th17' = grep("IL21|IL22|IL24|IL26|IL17A|IL17A|IL17F|IL17RA|IL10|RORC|RORA|STAT3|CCR4|CCR6|IL23RA|TGFB",
                cpdb_df$interacting_pair),
            'treg' = grep("IL35|IL10|FOXP3|IL2RA|TGFB", cpdb_df$interacting_pair),
            'costimulatory' = grep("CD86|CD80|CD48|LILRB2|LILRB4|TNF|CD2|ICAM|SLAM|LT[AB]|NECTIN2|CD40|CD70|CD27|CD28|CD58|TSLP|PVR|CD44|CD55|CD[1-9]",
                cpdb_df$interacting_pair),
            'coinhibitory' = grep("SIRP|CD47|ICOS|TIGIT|CTLA4|PDCD1|CD274|LAG3|HAVCR|VSIR",
                cpdb_df$interacting_pair),
            'niche' = grep("CSF", cpdb_df$interacting_pair)
            )

    interactions_sub = unique(unlist(gene_family_list[gene_family]))  
    cpdb_df = cpdb_df[interactions_sub, ]
   
    }
}

    #remove entirely insignificant interactions 
    message("removing entirely insignificant interactions")
    agp = aggregate(cpdb_df$pval, by = list('interaction' = cpdb_df$interaction), sum)
    psum = agp$x
    names(psum) = agp$interaction
    psum = psum == unique(table(cpdb_df$interaction))
    sig_interactions = names(psum)[!psum]
    cpdb_df = cpdb_df[cpdb_df$interaction %in% sig_interactions , ]

  #and render nonsignificant p values null
  cpdb_df$pval[cpdb_df$pval == 1] = NA
    cpdb_df$is_significant = cpdb_df$pval < 0.05
    
  
  #make sender receiver a factor
  #cpdb_df$SR <- factor(cpdb_df$SR, levels = unique(cpdb_df$SR))
  cpdb_df$SR = droplevels(cpdb_df$SR)
    
#scale the data
  mat <- matrix(0, ncol = length(unique(cpdb_df$LR)), nrow = length(unique(cpdb_df$SR)))
  colnames(mat) <- unique(cpdb_df$LR)
  rownames(mat) <- unique(cpdb_df$SR)
  message("clustering and/or scaling step...")

    for(i in 1:nrow(cpdb_df)){
    lr = as.character(cpdb_df[i, 'LR'])
    sr = as.character(cpdb_df[i, 'SR'])
    mat[sr, lr] = cpdb_df[i, "mean"]
    }

    #cpdb_df = cpdb_df[cpdb_df$LR %in% colnames(mat), ]
    
#scaling the data    
  if(scaling == "min_max"){
          message("min max normalising")

    mat <- apply(mat, 2, min_max_normalisation)
  }
  if(scaling == "scale"){
    mat <- scale(mat)
  }
  
              message("applying normalised values")

  melted_mat <- reshape2::melt(mat)
  rownames(melted_mat) <- paste0(melted_mat$Var1, melted_mat$Var2)
  cpdb_df$scaled <- melted_mat[paste0(cpdb_df$SR, cpdb_df$LR), "value"]  
  
  
  #now cluster the interactions
  if(cluster_interactions){
    message("clustering interactions")
    cpdb_df$LR <- factor(cpdb_df$LR, levels = colnames(mat)[hclust(dist(t(mat)), method = "ward.D")$order])
  }else{
    cpdb_df$LR <- factor(cpdb_df$LR, levels = unique(cpdb_df$LR))
  }


    #plot
  message("plotting result")
  if(scaling == "scale"){
    scaled_colors <- c("palegreen4", "grey80", "plum3")
    max_scaled <- max(abs(cpdb_df$scaled))
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = scaled )) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradient2(low = scaled_colors[1], mid = scaled_colors[2], high = scaled_colors[3],
                                                                 limits = c(-max_scaled, max_scaled ))
    
  }
  if(scaling == "min_max"){
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, 
                              size = scaled,
                              #size = -log10(pval +1e-3 ), 
                              fill = scaled)) + #geom_point(pch = 21) + 
      #geom_point(pch = 21) +
      geom_tile() +
      theme_bw() +
      scale_size_continuous(limits = c(0, 1)) + 
      #scale_size_continuous(limits = c(0, 1)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, 1)) + coord_fixed()  +
       # geom_point(pch = 21, data = subset(cpdb_df, is_significant), aes(x = SR, y = LR, fill = scaled, size = scaled), stroke = 2, color = 'red')
    geom_tile(data = subset(cpdb_df, is_significant), aes(x = SR, y = LR, fill = scaled, size = scaled), lwd = 2, color = 'grey50')

  }
  if(scaling == FALSE){
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = log1p(mean))) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, max(log1p(cpdb_df$mean)))) 
    
  }
  message("done")
  return(pl)

}
