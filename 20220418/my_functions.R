
markers <- readRDS("E:/KTakahashi/R_project/iwabuchi_sensei/LiverAnnotation/markers.rds")
mouse_liver_list = markers$best %>% group_by(type) %>% select(Mgene) %>% filter(Mgene >0) %>% 
  nest() %>% mutate(data = map(data, ~pull(.x, Mgene))) %>% deframe()



run_cluster_signature <- function(object, gene_list) {
  score_mt <- run_signature_value(object = object, gene_list = gene_list)
  score_mt %>% gather(-label, key = "signature", value = "exp") %>%
    group_by(signature, label) %>%
    summarise(pct = sum(exp>0)/n(), mean = mean(exp), sd = sd(exp)) %>%
    group_by(signature) %>%
    mutate(max = max(mean),
           score = mean/max,
           class_score = mean*score)
}


#


run_auto_labeling <- function(object = object, gene_list, interaction = TRUE) {
  
  use_df <- run_cluster_signature(object, gene_list = gene_list)
  
  sig_df <- use_df %>% 
    group_by(label) %>%
    top_n(1, class_score) %>% 
    dplyr::select(signature, label)
  
  bb <<- group_by(use_df, label) %>%
    top_n(3, class_score)
  label <- object@ident %>% 
    plyr::mapvalues(from = as.character(sig_df$label), to = sig_df$signature)
  
  #label <- interaction(object$seurat_clusters, label, sep = "_", lex.order = T, drop = T)
  if(interaction) label <- interaction(label, object@ident, sep = "_", lex.order = T, drop = T)
  
  return(label)
}




run_dot_plot <- function(object = data, gene_list, text_size = 8, order = FALSE, ...) {
  
  stopifnot(class(object) == "seurat")
  gene_list <- purrr::map(gene_list, ~.[. %in% rownames(object@data)])
  
  feature <- unique(unlist(gene_list))
  label_df <- enframe(gene_list, name = "label",value = "gene") %>% 
    unnest %>% mutate(label = as.factor(label))
  
  use_df <- run_gene_value(object = object, feature = feature)
  
  
  use_df <- use_df %>% left_join(label_df, by = c("gene"))
  
  n <- length(gene_list)
  label_color <- gg_color_hue(n)
  
  #label_color_use <- label_color[as.numeric(as.factor(label_df$label))]
  use_df$label_color <- label_color[as.numeric(as.factor(use_df$label))]
  use_df$gene <- fct_relevel(use_df$gene, feature)
  
  p <- use_df %>% ggplot(aes(cluster, fct_rev(gene), size = pct,  colour = score)) + geom_point() +
    scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                           values = c(1.0,0.7,0.6,0.4,0.3,0))
  
  if(order){
    p <- use_df %>% ggplot(aes(cluster, fct_reorder2(gene, score, cluster), size = pct,  colour = score)) + geom_point() +
      scale_colour_gradientn(colours = c("red","yellow","white","lightblue","darkblue"),
                             values = c(1.0,0.7,0.6,0.4,0.3,0))
  }
  
  use_color <- use_df %>% distinct(gene, label_color) %>%
    arrange(desc(gene)) %>%
    pull(label_color)

  p + theme(axis.text.y = element_text(colour = use_color, size = text_size),
            axis.text.x = element_text(angle = 90, vjust = 0.5) 
  ) + ylab("label_gene")
}


# calculate gene mean expression per specific cluster
run_gene_value <- function(object, feature) {
  cat(paste0("removed_gene: ", feature[!feature %in% rownames(object)]))
  feature <- unique(feature[feature %in% rownames(object@data)])
  
  cluster_label <- object@ident 
  
  use_df <- object@data[feature,]
  
  if(length(feature) ==1){
    use_df <- use_df %>% as.tibble()
  }else{
    use_df <- t(as.matrix(use_df)) %>% as.tibble()
  }
  
  use_df<- use_df %>% add_column(cluster = cluster_label)
  use_df <- use_df %>% gather(-cluster, key = "gene", value = "logCPM") %>%
    group_by(cluster, gene) %>% 
    summarise(avg_logCPM = mean(logCPM), 
              pct = sum(logCPM>0)/n()) %>%
    group_by(gene) %>% 
    mutate(score = avg_logCPM/max(avg_logCPM), mean = mean(avg_logCPM)) %>% ungroup() %>% 
    mutate(gene = fct_relevel(gene, feature )) 
  return(use_df)
}


#return summarise expression
run_signature_value <- function(object = data, gene_list = gene_list, use_func = "mean") {
  
  use_func <- switch (use_func, "mean" = mean, "gm_mean" = gm_mean1)
  gm_mean1 = function(a){prod(a)^(1/length(a))}
  
  mt <- object@meta.data
  gene <- rownames(object@data)
  gene_list <- purrr::map(gene_list, ~.[. %in% gene])
  gene_list <- gene_list[map(gene_list, length)>1]
  count_mt <- object@data
  for(i in seq_along(gene_list)){
    sub_mt <- count_mt[gene_list[[i]],]
    value <- apply(sub_mt,2, use_func)
    mt[names(gene_list)[i]] <- value
  }
  mt <- mt[names(gene_list)]
  mt["label"] <-  object@ident
  return(mt)
  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
