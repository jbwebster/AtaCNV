#' Calculate CNV score from CNV profile
#'
#' @param copy_ratio a cell-by-bin copy ratio matrix.
#' @param cluster cell cluster identity of each cell
#' @param delXY logical, whether to exclude chromosome X and Y
#'
#' @return CNV burden score for each cluster
#' @export
#'
#' @examples
#' @noRd
CNV_score <- function(copy_ratio, cluster, delXY=T){
  S <- c(1:length(cluster))
  m <- 1
  if(delXY){
    bin <- str2bin(colnames(CNV))
    copy_ratio <- copy_ratio[,!bin$chr %in% c("chrX","chrY")]
  }
  for(i in unique(cluster)){
    f <- cluster==i
    temp <- copy_ratio[f,]
    s <- mean((apply(temp,2,mean)-m)^2)
    S[f] <- s
    print(paste0("cluster ",i,": ",sum(f),"cells ","score=",s))
  }
  return(S)
}

#' Plot CNV heatmap using ComplexHeatmap
#'
#' @param copy_ratio a cell-by-bin copy ratio matrix. Column names should be
#' formatted.
#' @param del_chrXY logical, whether to exclude chromosome X and Y?
#' @param cell_cluster a vector indicating cell cluster identity of each cell,
#' or "kmeans" if you want to use kmeans to assign cell cluster, or "none" if
#' there is no predefined cluster. When cluster identity is assigned,
#' [cluster_within_group()] will be used to decide row hierarchy of the heatmap.
#' @param K number of clusters.
#' @param cell_annotation a vector of cell annotations to be shown on the left
#' side of heat map, or "none".
#' @param save_hc logical, whether to save hierarchical clustering result?
#' @param add_noise logical. If TRUE, add random noise to copy ratio to
#' avoid clustering error.
#' @param output_dir output directory
#' @param output_name output image name, should be ".png" or ".pdf"
#'
#' @return
#'
#' @export
#'
#' @examples
plot_heatmap <- function(copy_ratio,
                         del_chrXY=T,
                         cell_cluster="none",
                         cluster_method="kmeans", K=5,
                         cell_annotation="none",
                         save_hc=F,
                         add_noise=T,
                         output_dir="./",
                         output_name="copy_ratio.png"
){
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  save_dir <- paste0(output_dir, output_name)
  bins <- str2bin(colnames(copy_ratio))
  if(del_chrXY){
    f <- !(bins$chr %in% c("chrX","chrY","X","Y"))
    copy_ratio <- copy_ratio[,f]
    bins <- bins[f,]
  }
  chr_bkp <- bin2chrbkp(bins)

  if(add_noise){
    noise <- matrix(rnorm(n = length(copy_ratio), mean = 0, sd = 0.005),
                    nrow = nrow(copy_ratio),
                    ncol = ncol(copy_ratio))
    copy_ratio <- copy_ratio + noise
  }

  ## cluster
  print("Step 1: clustering...")
  if(min(copy_ratio)<=0){
    temp <- copy_ratio
  }else{
    temp <- log(copy_ratio+0.01)
  }

  if(identical(cell_cluster,"none")){
    d <- parallelDist::parallelDist(temp)
    hc_re <- hclust(d, method="ward.D2")
  }else if(identical(cell_cluster,"kmeans")){
    kmeans_pc <- 50
    if(identical(kmeans_pc, "none")){
      temp_ <- temp
    }else{
      temp_ <- irlba::prcomp_irlba(temp, n=kmeans_pc)$x
    }
    group <- kmeans(temp_, centers=K)$cluster
    hc_re <- ComplexHeatmap::cluster_within_group(t(temp), group)
  }else{
    cell_cluster <- as.numeric(as.factor(cell_cluster))
    hc_re <- ComplexHeatmap::cluster_within_group(t(temp), cell_cluster)
    K <- length(unique(cell_cluster))
  }

  ## plot
  print("Step 2: plotting...")
  plot_data <- copy_ratio
  color_thre=c(0.2,1.8)
  c <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")))(100)
  colormap <-
    circlize::colorRamp2(
      c(1:100)*(color_thre[2]-color_thre[1])*0.01+color_thre[1], c)

  ## set row annotation
  if(identical(cell_annotation,"none")){
    a <- ComplexHeatmap::HeatmapAnnotation(
      cell_type=ComplexHeatmap::anno_empty(border=F),which="row")
  }else if(all(cell_annotation %in% c("tumor","normal"))){
    a <- ComplexHeatmap::HeatmapAnnotation(
      cell_type=cell_annotation,
      col=list(cell_type=c("normal"="gray",
                           "tumor"="purple")),
      which="row",show_annotation_name=F)
  }else{
    a <- ComplexHeatmap::HeatmapAnnotation(
      cell_type=cell_annotation, which="row",show_annotation_name=F)
  }

  ## set col annotation
  n_chr <- length(unique(bins$chr))
  temp <- rep("gray",n_chr)
  temp[2*c(1:floor(n_chr/2))] <- "white"
  names(temp) <- unique(bins$chr)
  b <- ComplexHeatmap::HeatmapAnnotation(
    chr=bins$chr, col=list(chr=temp), which="column", show_annotation_name=F)

  ## heatmap
  # ComplexHeatmap::ht_opt$message <- FALSE
  ht <- ComplexHeatmap::Heatmap(
                plot_data, name="copy ratio", col=colormap,
                cluster_rows=hc_re, cluster_columns=FALSE,
                show_row_names=FALSE, show_column_names=FALSE,
                split=K,
                left_annotation=a, top_annotation=b,
                use_raster=T)

  ## add vertical lines
  k <- K
  if(substr(save_dir, nchar(save_dir)-2, nchar(save_dir))=="png"){
    png(filename=save_dir, width=2000, height=1000)
  }else{
    pdf(file=save_dir, width=20, height=10)
  }
  print(ht)
  for(i in 1:k){
    for(j in chr_bkp){
      ComplexHeatmap::decorate_heatmap_body(
        "copy ratio",
        {grid::grid.lines(c(j/dim(plot_data)[2], j/dim(plot_data)[2]),
                    c(0,1), gp=grid::gpar(lty=1, lwd=2))},
        slice=i)
      }
  }
  dev.off()

  ## save cluster result
  if(save_hc){
    saveRDS(hc_re, paste0(output_dir, "heatmap_hc.rds"))
  }
}
