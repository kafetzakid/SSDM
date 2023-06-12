#' @title Compute distance matrices for 1D lines
#'
#' @param data_lines 
#' @param p_mink 
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#' 
#' @author Danai Kafetzaki
#'
#' @examples
distance_matrix = function(data_lines, p_mink = 4) {
  
  distM_manh = matrix(0, nrow = dim(data_lines)[2], ncol = dim(data_lines)[2])
  distM_eucl = distM_ksts = distM_kend = distM_pear = distM_sper = distM_mink = distM_canb = distM_cheb = distM_manh
  
  for (i in 1:(dim(data_lines)[2]-1)) {
    for (j in (i+1):dim(data_lines)[2]) {
      
      distM_manh[i,j] = sum(abs(data_lines[,i] - data_lines[,j]))
      distM_eucl[i,j] = sum(abs(data_lines[,i] - data_lines[,j])^2)^(1/2)
      distM_cheb[i,j] = max(abs(data_lines[,i] - data_lines[,j]))
      distM_canb[i,j] = sum(abs(data_lines[,i] - data_lines[,j])/(abs(data_lines[,i]) + abs(data_lines[,j])), na.rm = TRUE)
      distM_mink[i,j] = sum(abs(data_lines[,i] - data_lines[,j])^p_mink)^(1/p_mink)
      distM_pear[i,j] = cor(cbind.data.frame(data_lines[,i],data_lines[,j]), method = "pearson")[1,2]
      distM_sper[i,j] = cor(cbind.data.frame(data_lines[,i],data_lines[,j]), method = "spearman")[1,2]
      distM_kend[i,j] = cor(cbind.data.frame(data_lines[,i],data_lines[,j]), method = "kendall")[1,2]
      distM_ksts[i,j] = ks.test(data_lines[,i], data_lines[,j])$statistic 
      
      distM_manh[j,i] = distM_manh[i,j]
      distM_eucl[j,i] = distM_eucl[i,j]
      distM_cheb[j,i] = distM_cheb[i,j]
      distM_canb[j,i] = distM_canb[i,j]
      distM_mink[j,i] = distM_mink[i,j]
      distM_pear[j,i] = distM_pear[i,j]
      distM_sper[j,i] = distM_sper[i,j]
      distM_kend[j,i] = distM_kend[i,j]
      distM_ksts[j,i] = distM_ksts[i,j]
    }
  }

  returns = list("distM_ksts" = distM_ksts, "distM_kend" = distM_kend, "distM_pear" = distM_pear, "distM_sper" = distM_sper, "distM_mink" = distM_mink, "distM_canb" = distM_canb, "distM_cheb" = distM_cheb, "distM_manh" = distM_manh, "distM_eucl" = distM_eucl)
  
  return(returns)
}



#' @title Compute distance matrices for GAM lines
#'
#' @param data_lines 
#' @param p_mink 
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#' 
#' @author Danai Kafetzaki
#'
#' @examples
distance_matrix_GAM = function(data_lines, p_mink = 4) {
  
  distM_manh = matrix(0, nrow = length(data_lines), ncol = length(data_lines))
  distM_eucl = distM_ksts = distM_kend = distM_pear = distM_sper = distM_mink = distM_canb = distM_cheb = distM_manh
  
  for (i in 1:(length(data_lines)-1)) {
    for (j in (i+1):length(data_lines)) {
      
      distM_manh[i,j] = sum(abs(data_lines[[i]][,2] - data_lines[[j]][,2]))
      distM_eucl[i,j] = sum(abs(data_lines[[i]][,2] - data_lines[[j]][,2])^2)^(1/2)
      distM_cheb[i,j] = max(abs(data_lines[[i]][,2] - data_lines[[j]][,2]))
      distM_canb[i,j] = sum(abs(data_lines[[i]][,2] - data_lines[[j]][,2])/(abs(data_lines[[i]][,2]) + abs(data_lines[[j]][,2])), na.rm = TRUE)
      distM_mink[i,j] = sum(abs(data_lines[[i]][,2] - data_lines[[j]][,2])^p_mink)^(1/p_mink)
      distM_pear[i,j] = cor(cbind.data.frame(data_lines[[i]][,2],data_lines[[j]][,2]), method = "pearson")[1,2]
      distM_sper[i,j] = cor(cbind.data.frame(data_lines[[i]][,2],data_lines[[j]][,2]), method = "spearman")[1,2]
      distM_kend[i,j] = cor(cbind.data.frame(data_lines[[i]][,2],data_lines[[j]][,2]), method = "kendall")[1,2]
      distM_ksts[i,j] = ks.test(data_lines[[i]][,2], data_lines[[j]][,2])$statistic 
      
      distM_manh[j,i] = distM_manh[i,j]
      distM_eucl[j,i] = distM_eucl[i,j]
      distM_cheb[j,i] = distM_cheb[i,j]
      distM_canb[j,i] = distM_canb[i,j]
      distM_mink[j,i] = distM_mink[i,j]
      distM_pear[j,i] = distM_pear[i,j]
      distM_sper[j,i] = distM_sper[i,j]
      distM_kend[j,i] = distM_kend[i,j]
      distM_ksts[j,i] = distM_ksts[i,j]
    }
  }
  
  returns = list("distM_ksts" = distM_ksts, "distM_kend" = distM_kend, "distM_pear" = distM_pear, "distM_sper" = distM_sper, "distM_mink" = distM_mink, "distM_canb" = distM_canb, "distM_cheb" = distM_cheb, "distM_manh" = distM_manh, "distM_eucl" = distM_eucl)
  
  return(returns)
}



#' @title Compute distance matrices for Metrics or Combo data 
#'
#' @param df 
#' @param p_mink 
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#' 
#' @author Danai Kafetzaki
#'
#' @examples
distance_matrix_df = function(df, p_mink = 4) {
  
  distM_cheb = matrix(0, nrow = dim(df)[1], ncol = dim(df)[1])
  distM_manh = distM_eucl = distM_ksts = distM_kend = distM_pear = distM_sper = distM_mink = distM_canb = distM_cheb
  
  for (i in 1:(dim(df)[1]-1)) {
    for (j in (i+1):dim(df)[1]) {
      
      distM_manh[i,j] = sum(abs(df[i,] - df[j,]))
      distM_eucl[i,j] = sum(abs(df[i,] - df[j,])^2)^(1/2)
      distM_cheb[i,j] = max(abs(df[i,] - df[j,]))
      distM_canb[i,j] = sum(abs(df[i,] - df[j,])/(abs(df[i,]) + abs(df[j,])), na.rm = T)
      distM_mink[i,j] = sum(abs(df[i,] - df[j,])^p_mink)^(1/p_mink)
      distM_pear[i,j] = cor(cbind.data.frame(df[i,], df[j,]), method = "pearson")[1,2]
      distM_sper[i,j] = cor(cbind.data.frame(df[i,], df[j,]), method = "spearman")[1,2]
      distM_kend[i,j] = cor(cbind.data.frame(df[i,], df[j,]), method = "kendall")[1,2]
      distM_ksts[i,j] = ks.test(as.numeric(df[i,]), as.numeric(df[j,]))$statistic 
      
      distM_manh[j,i] = distM_manh[i,j]
      distM_eucl[j,i] = distM_eucl[i,j]
      distM_cheb[j,i] = distM_cheb[i,j]
      distM_canb[j,i] = distM_canb[i,j]
      distM_mink[j,i] = distM_mink[i,j]
      distM_pear[j,i] = distM_pear[i,j]
      distM_sper[j,i] = distM_sper[i,j]
      distM_kend[j,i] = distM_kend[i,j]
      distM_ksts[j,i] = distM_ksts[i,j]
    }
  }
  returns = list("distM_ksts" = distM_ksts, "distM_kend" = distM_kend, "distM_pear" = distM_pear, "distM_sper" = distM_sper, "distM_mink" = distM_mink, "distM_canb" = distM_canb, "distM_cheb" = distM_cheb, "distM_manh" = distM_manh, "distM_eucl" = distM_eucl)
  
  return(returns)
}


  


#' @title Function for STAD 2D coordinates
#'
#' @importFrom stats as.dist
#' @importFrom stad stad
#' @importFrom stad plot_graph
#' @importFrom igraph layout_nicely
#' @importFrom igraph as_data_frame
#'
#' @param distM 
#' @param cols
#' @param seed
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#'
#' @author Danai Kafetzaki
#'
#' @examples
library(stad)
library(igraph)
stad_reconstruct = function(distM, cols, seed, show_plot = TRUE, vs = 3, ew = 0.2){

  distM = distM/max(distM)
  distM = stats::as.dist(distM)
  
  set.seed(seed)
  stad_list = stad::stad(distM)
  
  val_STAD_corr = max(stad_list$trace$correlation)
  
  if (show_plot) {
    par(mfrow = c(1,1))
    p1 = stad::plot_graph(stad_list, layout = igraph::layout_nicely,
                          vertex.color = cols,
                          vertex.frame.color = cols,
                          vertex.size = vs, edge.width = ew)
  }
  
  df_links = igraph::as_data_frame(stad_list$graph)
  names(df_links)[1] = "Source"; names(df_links)[2] = "Target"
  df_links$Source = as.numeric(df_links$Source)
  df_links$Target = as.numeric(df_links$Target)
  df_links$Value2 = ceiling(100*df_links$value)
  df_links = df_links[, colnames(df_links) %in% c("Source" , "Target", "Value2")]
  df_coord_nicely = igraph::layout_nicely(stad_list$graph)
  
  returns = list("df_coord_nicely" = df_coord_nicely, "df_links" = df_links, "graph" = stad_list$graph, "val_STAD_corr" = val_STAD_corr)
  
  return(returns)
  
}


#' @title Function for Shepard diagram analysis
#'
#' @importFrom MASS Shepard
#'
#' @param dfM1 
#' @param dfM2 
#' @param show_plot 
#' @param pchsel
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#' 
#' @author Danai Kafetzaki
#'
#' @examples
library(MASS)
shep_diagr_analysis = function(dfM1, dfM2, show_plot = TRUE, pchsel = ".") {
  
  plot.shepard = MASS::Shepard(dfM1, dfM2, p = 2)
  # Non-parametric fit
  val_Rsqr_np = round(cor(plot.shepard$yf, plot.shepard$x),3)
  # Regression model fit
  lmmod = lm(y~x-1, data = plot.shepard)
  lmmodsum = summary(lmmod);
  if (show_plot) {
    plot(plot.shepard, pch = pchsel,  col = "#00000080", main = "Shepard diagram for STAD (2D layout - naive)")
    lines(plot.shepard$x, plot.shepard$yf, type = "l", lwd = 3, col = "#991450")
    lines(plot.shepard$x, lmmod$fitted.values, type = "l", lwd = 3, col = "#ab7f07")
    text(x = 0.3*mean(plot.shepard$x), y = 0.9*max(plot.shepard$y), cex = 1.6, bquote(paste('R'^2,' = ',.(round(lmmodsum$r.squared,3)))), col = "#ab7f07")
  }
  val_Rsqr_lm = round(lmmodsum$r.squared,3)

  returns = list("val_Rsqr_lm" = val_Rsqr_lm, "val_Rsqr_np" = val_Rsqr_np)
  
  return(returns)
  
  
  }
  
  
# Function for STAD Dikstra's reconstruction 
# 
#' @title Adjacency matrix comparison 
#' 
#' @importFrom igraph shortest.paths
#'
#' @param g the graph object
#' @param distM the distance matrix used to create the graph
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#'
#' @examples
library(igraph)
compare_matrices = function(g, distM) {
  
  distM = distM/max(distM)
  
  pathsM = igraph::shortest.paths(g, algorithm = "dijkstra")
  pathsM = pathsM/max(pathsM)
  
  diffM = distM - pathsM

  dist_vect = round(abs(c(diffM[upper.tri(diffM)])),3)
  hist(dist_vect, breaks = 20, col = "#c1c29d", border = "#000000")
  diff_summary = summary(dist_vect)
  diff_mad = mad(dist_vect)
  diff_sd = sd(dist_vect)

  returns = list("diff_summary" = diff_summary, "diff_mad" = diff_mad, "diff_sd" = diff_sd)
  return(returns)
  
}


#' @title Distance metric selection
#'
#' @param df_results
#' @param iteration
#'
#' @return A list with the following items:
#' \itemize{
#'   \item itemname - item type. Description.
#'   \item itemname - item type. Description.
#' }
#' @export
#'
#' @examples
library(ggplot2)
library(stats)
library(car)
library(FSA)
distance_metric_selection = function(df_results, iteration = TRUE, criterion = "median_abs_difference", stattest = TRUE, show_plot = TRUE, ggp_title = "") {
  
  stat_result = NULL
  posthoc_result = NULL
  colnames(df_results)[which(colnames(df_results) == criterion)] = "criterion"
  
  if (iteration) {
    df_results_median = aggregate(criterion ~ metric, data = df_results, median)
    df_results_median = df_results_median[order(df_results_median$criterion, decreasing = FALSE), ]
    selected_metric = df_results_median[1,"metric"]
    metric_value = df_results_median[1,"criterion"]
    
    ggp = ggplot(data = df_results, aes(x = criterion, y = metric)) + 
            geom_boxplot() + 
            labs(x = "criterion",
                 y = 'distance metric',
                 title = as.character(ggp_title)) +
            theme(legend.position = "none",
                  legend.direction = "vertical",
                  legend.text=element_text(size = 28, face = "italic"),
                  legend.title = element_blank(),
                  legend.background = element_blank(),
                  legend.box.background = element_rect(colour = "grey10"),
                  legend.key.width = unit(1.8,"cm"),
                  panel.background = element_rect(fill = "white",
                                                  colour = "white",
                                                  size = 0.25, linetype = "solid"),
                  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                                  colour = "grey80"),
                  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                                  colour = "grey80"),
                  strip.text.x = element_text(size = 18, face = "plain"),
                  axis.text.x = element_text(color = "grey10", size = 24, face = "plain", angle = 0, hjust = 0.5),
                  axis.text.y = element_text(color = "grey10", size = 24, face = "plain"),
                  axis.title.x = element_text(color = "grey10", size = 26, face = "plain"),
                  axis.title.y = element_text(color = "grey10", size = 26, face = "plain"),
                  plot.title = element_text(colour = "black", size = 20, face = "italic", hjust = 0.5))
    if (show_plot) {print(ggp)}
    
    if (stattest) {
      data_test = data.frame("metric" = df_results$metric, "criterion" = df_results$criterion)
      data_test$metric = as.factor(data_test$metric) 
      stat_result = stats::aov(criterion ~ metric-1, data_test)
      tukey.t = stats::TukeyHSD(stat_result)
      posthoc_result = as.data.frame(tukey.t$metric)
      norm_result = stats::shapiro.test(residuals(stat_result))
      homo_result = car::leveneTest(stat_result)
      
      if (norm_result$p.value < 0.05 | homo_result$`Pr(>F)`[1] < 0.05) {
        stat_result = stats::kruskal.test(criterion ~ metric-1, data = data_test)
        dunnResult = FSA::dunnTest(criterion ~ metric-1, data = data_test, method = "bh")
        posthoc_result = as.data.frame(dunnResult$res)
      }
    }
    
  } else {
    if (stattest) {print("Statistically testing the difference between the performance of metrics is not possible without iterations.")}
    df_results_median = NULL
    df_results = df_results[order(df_results$criterion, decreasing = FALSE), ]
    selected_metric = df_results[1,"metric"]
    metric_value = df_results[1,"criterion"]
  }
  
  returns = list("df_results_median" = df_results_median, "selected_metric" = selected_metric, "metric_value" = metric_value, "comparison_plot" = ggp, "statistical_test_result" = list(stat_result), "posthoc_test_result" = list(posthoc_result))
  
  return(returns)
}



