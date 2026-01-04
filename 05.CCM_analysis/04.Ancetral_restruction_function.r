#' 更改文档格式
#'
#' @description
#' 本函数主要能将文件中乱序的物种排序成与进化树相同的顺序，便于后续进行分析
#'
#' @param tree  进化树格式文件（.nwk）
#' @param file  需要分析的文档（先转化到一个变量内）
#' @param plot_format  选择图片输出形式，默认为png文件
#' @param tree_type 决定进化树的形状，默认为"phylogram"，可选"cladogram","fan"等
#' @param width 画布宽度
#' @param height 画布高度
#' @param cex 调整进化树子节点字体大小
#' @export
#' @return 修改后的file矩阵
file.handle <- function(
    tree,
    file,
    plot_format="png",
    tree_type="phylogram",
    width,
    height,
    cex) {

  suppressPackageStartupMessages({
    library(ggtree)
    library(treeio)
    library(ape)
  })

  inner.nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
  if (plot_format == "png") {
    png(file.path(path, "tree.png"), width = width, height = height)
  } else if (plot_format == "pdf") {
    pdf(file.path(path, "tree.pdf"), width = width, height = height)
  } else {
    stop("Unsupported plot format. Use 'png' or 'pdf'.")
  }
  plot(tree,type = tree_type,use.edge.length = F,
       cex = cex,label.offset = 0.01, node.depth = 2)
  nodelabels(text=inner.nodes)
  dev.off()

  tree.un <- unroot(tree)
  inner.nodes.notheria <- (length(tree.un$tip.label) + 1):(length(tree.un$tip.label) + tree.un$Nnode)
  if (plot_format == "png") {
    png(file.path(path, "tree_un.png"), width = width, height = height)
  } else if (plot_format == "pdf") {
    pdf(file.path(path, "tree_un.pdf"), width = width, height = height)
  } else {
    stop("Unsupported plot format. Use 'png' or 'pdf'.")
  }
  plot(tree.un,type = tree_type,use.edge.length = F,
       cex = cex,label.offset = 0.01, node.depth = 2)
  nodelabels(text=inner.nodes.notheria)
  dev.off()

  last_colindex <- ncol(file)
  file <- file[, c(last_colindex, 1:(last_colindex - 1))]
  sort_order <- c("id",tree$tip.label)
  sort_order
  file <- file[,sort_order]

  return(file)
}

#' 计算基因转化率
#' @param fulldf  已被ASR类函数处理过的矩阵
#' @param tree  进化树
#' @export
#' @return 返回矩阵，包括gain;loss;changes;branch_length,changes_norm,gains_norm.loss_norm
calculateTransitionRates <- function(fulldf, tree) {
  n_branches <- nrow(tree$edge)
  results_df <- matrix(NA, nrow = n_branches, ncol = 7)
  colnames(results_df) <- c("gains", "losses", "changes", "branch_length", "changes_norm", "gains_norm", "losses_norm")

  branchlabels <- apply(tree$edge, 1, function(x) paste0("Branch ", x[1], "->", x[2]))
  rownames(results_df) <- branchlabels

  for (branch in 1:n_branches) {
    parent_node <- tree$edge[branch, 1]
    child_node <- tree$edge[branch, 2]

    deltas <- fulldf[, parent_node] - fulldf[, child_node]
    gains <- sum(deltas == -1)
    losses <- sum(deltas == 1)
    changes <- gains + losses

    branch_length <- tree$edge.length[branch]

    results_df[branch, 1] <- gains
    results_df[branch, 2] <- losses
    results_df[branch, 3] <- changes
    results_df[branch, 4] <- branch_length
    results_df[branch, 5] <- ifelse(!is.na(changes) && branch_length > 0, changes / branch_length, NA)
    results_df[branch, 6] <- ifelse(!is.na(gains) && branch_length > 0, gains / branch_length, NA)
    results_df[branch, 7] <- ifelse(!is.na(losses) && branch_length > 0, losses / branch_length, NA)
  }

  return(results_df)
}

#' 计算获得/损失性状以及节点表达性状
#' @param tree  已被ASR类函数处理过的矩阵
#' @param file  进化树
#' @export
#' @return 返回矩阵，包括gain;loss;changes;gain_gene;loss_gene;active_gene
calculateGenesetDynamics <- function(fulldf, tree) {

  tip <- as.matrix(tree$tip.label)
  mapping <- setNames(tip[1:10], 1:10)

  ntree <- sapply(tree$edge, function(x) {
    if (x %in% 1:10) mapping[as.character(x)] else x
  })
  ntree <- matrix(ntree, ncol = 2)

  result_df <- data.frame(
    Type = character(0),
    ID = character(0),
    Gains = integer(0),
    Losses = integer(0),
    Changes = integer(0),
    Gains_Genes = character(0),
    Losses_Genes = character(0),
    Active_Genes = character(0),
    stringsAsFactors = FALSE
  )

  for (branch in 1:nrow(tree$edge)) {
    parent_node <- tree$edge[branch, 1]
    child_node <- tree$edge[branch, 2]

    deltas <- fulldf[, parent_node] - fulldf[, child_node]
    gains <- sum(deltas == -1)
    losses <- sum(deltas == 1)
    changes <- gains + losses
    loss_genes <- rownames(fulldf)[deltas == 1]
    gain_genes <- rownames(fulldf)[deltas == -1]

    branch_id <- paste(
      ifelse(parent_node %in% 1:10, mapping[as.character(parent_node)], parent_node),
      "->",
      ifelse(child_node %in% 1:10, mapping[as.character(child_node)], child_node)
    )

    result_df <- rbind(result_df, data.frame(
      Type = "Branch",
      ID = branch_id,
      Gains = gains,
      Losses = losses,
      Changes = changes,
      Gains_Genes = if (gains > 0) paste(gain_genes, collapse = ", ") else "",
      Losses_Genes = if (losses > 0) paste(loss_genes, collapse = ", ") else "",
      Active_Genes = "",
      stringsAsFactors = FALSE
    ))
  }

  for (node in colnames(fulldf)) {

    active_genes <- rownames(fulldf)[fulldf[, node] == 1]
    active_genes_str <- if (length(active_genes) > 0)
      paste(sort(active_genes), collapse = ", ") else ""

    result_df <- rbind(result_df, data.frame(
      Type = "Node",
      ID = node,
      Gains = NA,
      Losses = NA,
      Changes = NA,
      Gains_Genes = "",
      Losses_Genes = "",
      Active_Genes = active_genes_str,
      stringsAsFactors = FALSE
    ))
  }

  return(result_df)
}

#' 替换字符
#'
#' @description
#' 本函数主要用于将数字内容进行对应字符转换
#'
#' @param number_string  数字矩阵
#' @param gene_list  基因列表
#' @export
#' @return 返回矩阵，用基因字符替换数字
convert_numbers_to_genes <- function(number_string, gene_list) {
  if (is.na(number_string) || number_string == "") {
    return("")
  }
  numbers <- as.numeric(unlist(strsplit(number_string, ",\\s*")))
  genes <- gene_list[numbers]
  paste(genes, collapse = ", ")
}



#' 最大简约法(等权重)
#' @param tree  进化树
#' @param data  符合要求的矩阵
#' @param output_path 输出路径
#' @param plot_format  选择图片输出形式，默认为png文件
#' @param tree_type 决定进化树的形状，默认为"phylogram"，可选"cladogram","fan"等
#' @param width 画布宽度
#' @param height 画布高度
#' @param cex 调整进化树子节点字体大小
#' @param cellname 细胞名称
#' @param convert_genes 默认TRUE,性状多的用来循环输出
#' @export
#' @return 两个.csv文件，一个为经主要函数处理后的文档,一个为计算获得/失去性状的文档,另保存一个.png文件,返回进化树各边变化情况
MPWeight11 <- function(
    tree,
    data,
    output_path,
    plot_format = "png",
    tree_type = "phylogram",
    width,
    height,
    cex,
    cellname,
    convert_genes = TRUE) {

  suppressPackageStartupMessages({
    library(castor)
    library(phytools)
    library(ape)
  })

  tip_labels <- tree$tip.label
  pre_data <- data[, tip_labels]
  transposed_data <- t(pre_data)

  fulldf <- data.frame()
  inner_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)

  for (igene in rownames(pre_data)) {
    x <- unlist(pre_data[igene, ])
    xmap <- replace(x, x == 1, 2)
    xmap <- replace(xmap, xmap == 0, 1)
    xres <- asr_max_parsimony(
      tree, xmap, Nstates = 2,
      transition_costs = "all_equal",
      edge_exponent = 0,
      weight_by_scenarios = TRUE,
      check_input = TRUE
    )
    xal <- xres$ancestral_likelihoods
    xwin <- as.integer(xal[, 2] > xal[, 1])
    names(xwin) <- inner_nodes
    xmerge <- c(x, xwin)
    xdf <- t(data.frame(xmerge))
    rownames(xdf) <- igene
    fulldf <- rbind(fulldf, xdf)
  }

  write.csv(
    fulldf,
    file.path(output_path, paste0(cellname, "_gene_node_forecast_MPweight11.csv")),
    row.names = TRUE
  )

  results_df <- calculateTransitionRates(fulldf, tree)

  plot_file <- file.path(
    output_path,
    paste0(cellname, "_ASR_MPweight11_tree.", plot_format)
  )

  if (plot_format == "png") {
    png(plot_file, width = width, height = height)
  } else if (plot_format == "pdf") {
    pdf(plot_file, width = width, height = height)
  } else {
    stop("Unsupported plot format. Use 'png' or 'pdf'.")
  }

  plot.phylo(tree, type = tree_type,edge.color = "black", edge.width = 2,
             use.edge.length = F,cex = cex,label.offset = 0.01, node.depth = 2)
  par(xpd = TRUE)
  edgelabels(
    text = paste0("+", results_df[, "gains"]),
    col = "green4", bg = "white", adj = c(0.5, -0.1)
  )
  edgelabels(
    text = paste0("-", results_df[, "losses"]),
    col = "red2", bg = "white", adj = c(0.5, 1.3)
  )

  dev.off()

  result <- calculateGenesetDynamics(fulldf, tree)
  result_df <- as.data.frame(result)
  gene_names <- data
  gene_names <- gene_names[[1]]

  convert_numbers_to_genes <- function(number_string, gene_list) {
    if (is.na(number_string) || number_string == "") {
      return("")
    }
    numbers <- as.numeric(unlist(strsplit(number_string, ",\\s*")))
    genes <- gene_list[numbers]
    paste(genes, collapse = ", ")
  }

  if (convert_genes) {
    columns_to_convert <- c("Gains_Genes", "Losses_Genes", "Active_Genes")
    for (col in columns_to_convert) {
      if (col %in% colnames(result_df)) {
        result_df[[col]] <- sapply(
          result_df[[col]],
          convert_numbers_to_genes,
          gene_list = gene_names
        )
      }
    }
  }

  write.csv(
    result_df,
    file.path(output_path, paste0(cellname, "_gene_expression_MPweight11.csv")),
    row.names = FALSE
  )

  return(print("success"))
}


#' 最大简约法(二倍权重)
#' @param tree  进化树
#' @param data  符合要求的矩阵
#' @param output_path 输出路径
#' @param plot_format 选择图片输出形式，默认为png文件
#' @param tree_type 决定进化树的形状，默认为"phylogram"，可选"cladogram","fan"等
#' @param width 画布宽度
#' @param height 画布高度
#' @param cex 调整进化树子节点字体大小
#' @param cellname 细胞名称
#' @param convert_genes 默认TRUE,性状多的用来循环输出
#' @export
#' @return 两个.csv文件，一个为经主要函数处理后的文档,一个为计算获得/失去性状的文档,另保存一个.png文件,返回进化树各边变化情况
MPWeight21 <- function(
    tree,
    data,
    output_path,
    plot_format = "png",
    tree_type = "phylogram",
    width,
    height,
    cex,
    cellname,
    gain=2,
    lost=1,
    convert_genes = TRUE) {

  suppressPackageStartupMessages({
    library(castor)
    library(phytools)
    library(ape)
  })

  tip_labels <- tree$tip.label
  pre_data <- data[, tip_labels]
  transposed_data <- t(pre_data)

  custom_costs <- matrix(c(0, gain,
                           lost, 0),
                         nrow = 2, byrow = TRUE)

  fulldf <- data.frame()
  inner_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)

  for (igene in rownames(pre_data)){
    x <- unlist(pre_data[igene,])
    xmap <- replace(x, x == 1, 2)
    xmap <- replace(xmap, xmap == 0, 1)
    xres <- asr_max_parsimony(tree, xmap, Nstates=2,
                              transition_costs = custom_costs,
                              edge_exponent=0, weight_by_scenarios=TRUE,
                              check_input=TRUE)
    xal <- xres$ancestral_likelihoods
    xwin <- xal[,2] > xal[,1]
    xwin <- replace(xwin, xwin == TRUE, 1)
    xwin <- replace(xwin, xwin == FALSE, 0)
    names(xwin) <- inner_nodes
    xmerge <- c(x, xwin)
    xdf <- t(data.frame(xmerge))
    rownames(xdf) <- igene
    fulldf <- rbind(fulldf, xdf)
  }

  write.csv(
    fulldf,
    file.path(output_path, paste0(cellname, "_gene_node_forecast_MPweight21.csv")),
    row.names = TRUE
  )

  results_df <- calculateTransitionRates(fulldf, tree)

  plot_file <- file.path(
    output_path,
    paste0(cellname, "_ASR_MPweight21_tree.", plot_format)
  )

  if (plot_format == "png") {
    png(plot_file, width = width, height = height)
  } else if (plot_format == "pdf") {
    pdf(plot_file, width = width, height = height)
  } else {
    stop("Unsupported plot format. Use 'png' or 'pdf'.")
  }

  plot.phylo(tree, type = tree_type,edge.color = "black", edge.width = 2,
             use.edge.length = F,cex = cex,label.offset = 0.01, node.depth = 2)
  par(xpd = TRUE)
  
  n_tips <- length(tree$tip.label)
  non_leaf_edges <- which(tree$edge[, 2] > n_tips) # tree$edge[,2]是子节点编号
  
  edgelabels(
    text = paste0("+", results_df[non_leaf_edges, "gains"]),
    col = "green4", bg = "white", adj = c(0.5, -0.1),
    edge = non_leaf_edges
  )
  edgelabels(
    text = paste0("-", results_df[non_leaf_edges, "losses"]),
    col = "red2", bg = "white", adj = c(0.5, 1.3),
    edge = non_leaf_edges
  )

  dev.off()

  result <- calculateGenesetDynamics(fulldf, tree)
  result_df <- as.data.frame(result)
  gene_names <- data
  gene_names <- gene_names[[1]]

  convert_numbers_to_genes <- function(number_string, gene_list) {
    if (is.na(number_string) || number_string == "") {
      return("")
    }
    numbers <- as.numeric(unlist(strsplit(number_string, ",\\s*")))
    genes <- gene_list[numbers]
    paste(genes, collapse = ", ")
  }

  if (convert_genes) {
    columns_to_convert <- c("Gains_Genes", "Losses_Genes", "Active_Genes")
    for (col in columns_to_convert) {
      if (col %in% colnames(result_df)) {
        result_df[[col]] <- sapply(
          result_df[[col]],
          convert_numbers_to_genes,
          gene_list = gene_names
        )
      }
    }
  }

  write.csv(
    result_df,
    file.path(output_path, paste0(cellname, "_gene_expression_MPweight21.csv")),
    row.names = FALSE
  )

  return(print("success"))
}


#' Wanger法
#' @param tree  进化树
#' @param data  符合要求的矩阵
#' @param out_group 输入外群物种
#' @param output_path 输出路径
#' @param plot_format 选择图片输出形式，默认为png文件
#' @param tree_type 决定进化树的形状，默认为"phylogram"，可选"cladogram","fan"等
#' @param width 画布宽度
#' @param height 画布高度
#' @param cex 调整进化树子节点字体大小
#' @param cellname 细胞名称
#' @param convert_genes 默认TRUE,性状多的用来循环输出
#' @export
#' @return 两个.csv文件，一个为经主要函数处理后的文档,一个为计算获得/失去性状的文档,另保存一个.png文件,返回进化树各边变化情况
Wanger <- function(
    tree,
    data,
    out_group,
    output_path,
    plot_format = "png",
    tree_type = "phylogram",
    width,
    height,
    cex,
    cellname,
    convert_genes = TRUE) {

  suppressPackageStartupMessages({
    library(castor)
    library(phytools)
    library(ape)
  })

  tree.un <- unroot(tree)

  tip_labels <- tree.un$tip.label
  pre_data <- data[, tip_labels]
  transposed_data <- t(pre_data)

  fulldf <- data.frame()
  inner_nodes <- (length(tree.un$tip.label) + 1):(length(tree.un$tip.label) + tree.un$Nnode)

  for (igene in rownames(pre_data)){
    x <- unlist(pre_data[igene,])
    xmap <- x
    xres <- MPR(xmap, tree.un, outgroup = out_group)
    xwin <- xres[,2]
    names(xwin) <- inner_nodes
    xmerge <- c(x, xwin)
    xdf <- t(data.frame(xmerge))
    rownames(xdf) <- igene
    fulldf <- rbind(fulldf, xdf)
  }

  write.csv(
    fulldf,
    file.path(output_path, paste0(cellname, "_gene_node_forecast_Wanger.csv")),
    row.names = TRUE
  )

  results_df <- calculateTransitionRates(fulldf, tree.un)

  plot_file <- file.path(
    output_path,
    paste0(cellname, "_Wanger_tree.", plot_format)
  )

  if (plot_format == "png") {
    png(plot_file, width = width, height = height)
  } else if (plot_format == "pdf") {
    pdf(plot_file, width = width, height = height)
  } else {
    stop("Unsupported plot format. Use 'png' or 'pdf'.")
  }

  plot.phylo(tree, type = tree_type,edge.color = "black", edge.width = 2,
             use.edge.length = F,cex = cex,label.offset = 0.01, node.depth = 2)
  par(xpd = TRUE)
  edgelabels(
    text = paste0("+", results_df[, "gains"]),
    col = "green4", bg = "white", adj = c(0.5, -0.1)
  )
  edgelabels(
    text = paste0("-", results_df[, "losses"]),
    col = "red2", bg = "white", adj = c(0.5, 1.3)
  )

  dev.off()

  result <- calculateGenesetDynamics(fulldf, tree.un)
  result_df <- as.data.frame(result)
  gene_names <- data
  gene_names <- gene_names[[1]]

  convert_numbers_to_genes <- function(number_string, gene_list) {
    if (is.na(number_string) || number_string == "") {
      return("")
    }
    numbers <- as.numeric(unlist(strsplit(number_string, ",\\s*")))
    genes <- gene_list[numbers]
    paste(genes, collapse = ", ")
  }

  if (convert_genes) {
    columns_to_convert <- c("Gains_Genes", "Losses_Genes", "Active_Genes")
    for (col in columns_to_convert) {
      if (col %in% colnames(result_df)) {
        result_df[[col]] <- sapply(
          result_df[[col]],
          convert_numbers_to_genes,
          gene_list = gene_names
        )
      }
    }
  }

  write.csv(
    result_df,
    file.path(output_path, paste0(cellname, "_gene_expression_Wanger.csv")),
    row.names = FALSE
  )

  return(print("success"))
}


#' REML法
#' @param tree  进化树
#' @param data  符合要求的矩阵
#' @param output_path 输出路径
#' @param plot_format 选择图片输出形式，默认为png文件
#' @param tree_type 决定进化树的形状，默认为"phylogram"，可选"cladogram","fan"等
#' @param width 画布宽度
#' @param height 画布高度
#' @param cex 调整进化树子节点字体大小
#' @param cellname 细胞名称
#' @param convert_genes 默认TRUE,性状多的用来循环输出
#' @export
#' @return 两个.csv文件，一个为经主要函数处理后的文档,一个为计算获得/失去性状的文档,另保存一个.png文件,返回进化树各边变化情况
REML <- function(
    tree,
    data,
    output_path,
    plot_format = "png",
    tree_type = "phylogram",
    width,
    height,
    cex,
    cellname,
    convert_genes = TRUE) {

  suppressPackageStartupMessages({
    library(castor)
    library(phytools)
    library(ape)
  })

  tip_labels <- tree$tip.label
  pre_data <- data[, tip_labels]
  transposed_data <- t(pre_data)

  fulldf <- data.frame()
  inner_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)

  for (igene in rownames(pre_data)){
    x <- unlist(pre_data[igene,])
    xres <- ace(x, tree, type='continuous', method="REML")
    xal <- xres$ace
    xwin <- as.integer(xres$ace > 0.5)
    names(xwin) <- inner_nodes
    xmerge <- c(x, xwin)
    xdf <- t(data.frame(xmerge))
    rownames(xdf) <- igene
    fulldf <- rbind(fulldf, xdf)
  }

  write.csv(
    fulldf,
    file.path(output_path, paste0(cellname, "_gene_node_forecast_REML.csv")),
    row.names = TRUE
  )

  results_df <- calculateTransitionRates(fulldf, tree)

  plot_file <- file.path(
    output_path,
    paste0(cellname, "_ASR_REML_tree.", plot_format)
  )

  if (plot_format == "png") {
    png(plot_file, width = width, height = height)
  } else if (plot_format == "pdf") {
    pdf(plot_file, width = width, height = height)
  } else {
    stop("Unsupported plot format. Use 'png' or 'pdf'.")
  }

  plot.phylo(tree, type = tree_type,edge.color = "black", edge.width = 2,
             use.edge.length = F,cex = cex,label.offset = 0.01, node.depth = 2)
  par(xpd = TRUE)
  edgelabels(
    text = paste0("+", results_df[, "gains"]),
    col = "green4", bg = "white", adj = c(0.5, -0.1)
  )
  edgelabels(
    text = paste0("-", results_df[, "losses"]),
    col = "red2", bg = "white", adj = c(0.5, 1.3)
  )

  dev.off()

  result <- calculateGenesetDynamics(fulldf, tree)
  result_df <- as.data.frame(result)
  gene_names <- data
  gene_names <- gene_names[[1]]

  convert_numbers_to_genes <- function(number_string, gene_list) {
    if (is.na(number_string) || number_string == "") {
      return("")
    }
    numbers <- as.numeric(unlist(strsplit(number_string, ",\\s*")))
    genes <- gene_list[numbers]
    paste(genes, collapse = ", ")
  }

  if (convert_genes) {
    columns_to_convert <- c("Gains_Genes", "Losses_Genes", "Active_Genes")
    for (col in columns_to_convert) {
      if (col %in% colnames(result_df)) {
        result_df[[col]] <- sapply(
          result_df[[col]],
          convert_numbers_to_genes,
          gene_list = gene_names
        )
      }
    }
  }

  write.csv(
    result_df,
    file.path(output_path, paste0(cellname, "_gene_expression_REML.csv")),
    row.names = FALSE
  )

  return(print("success"))
}
