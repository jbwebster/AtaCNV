#' Determine CNV breakpoints by BICseq2 and output single cell copy ratio profile
#'
#' @param norm_count a normalized read count matrix. Obtained from previous step
#' \code{normalize()}.
#' @param baseline read count baseline of normal cells. Obtained from previous step
#' \code{normalize()}.
#' @param BICseq path to BICseq executable file.
#' @param lambda parameter of BICseq segmentation algorithm that controls the overall
#' resolution of segmentation result. In general, more breakpoints will be reported
#' when smaller \code{lambda} is used.
#' @param min_interval minimum bin number between breakpoints.
#' @param smooth logical. Whether to aggregate similar cells to smooth final CNV
#' results?
#' @param downsample an integer. If larger than 0, rows of input matrix will be
#' randomly down-sampled to this number to accelerate segmentation. 1000 is
#' recommended if the sample contains too many cells.
#' @param output_dir output directory.
#' The generated temporary files will also be output to this directory.
#' @param output_name output file name, should be ".rds"
#'
#' @return Return a list of the following objects:
#' \itemize{
#' \item "copy_ratio": a cell-by-bin copy ratio matrix.
#' \item "bkp": breakpoint positions.
#' \item "CNV_seg": a cell-by-segment copy ratio matrix.
#' @export
#'
#' @examples
#'
calculate_CNV <- function(norm_count,
                          baseline,
                          BICseq="default",
                          lambda=5, min_interval=2, smooth=F,
                          downsample=0,
                          output_dir="./",
                          output_name="CNV_result.rds"
){
  if(BICseq=="default"){
    if(Sys.info()["sysname"]=="Windows"){
      BICseq <- system.file("BICseq2", "mbic-seq.exe", package="AtaCNV")
    }else if(Sys.info()["sysname"]=="Linux"){
      BICseq <- system.file("BICseq2", "MBICseq", package="AtaCNV")
    }

  }
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  if(is.vector(baseline)){
    baseline <-
      matrix(data=baseline, nrow=nrow(norm_count), ncol=ncol(norm_count), byrow=T)
  }
  bins <- str2bin(colnames(norm_count))
  chr_bkp <- bin2chrbkp(bins)

  ## find breakpoints
  print("Running BICseq2 for segmentation...")
  bkp <- c()
  lambda <- lambda
  tmp_dir <- output_dir
  k <- 0
  if(downsample>0){
    idx <- sample(1:nrow(norm_count), downsample)
    norm_count_origin <- norm_count
    norm_count <- norm_count[idx,]
    baseline_origin <- baseline
    baseline <- baseline[idx,]
  }
  for(i in unique(bins$chr)){
    print(i)
    if(sum(bins$chr==i)<=2){
      next
    }
    norm_count_ <- norm_count[,bins$chr==i]
    bic_input <- matrix(ncol=2*nrow(norm_count_)+2,
                        nrow=ncol(norm_count_))
    bic_input[,1:2] <- as.matrix(str2bin(colnames(norm_count_))[,2:3])
    bic_input[,2*c(1:nrow(norm_count_))+1] <- t(norm_count_)
    bic_input[,2*c(1:nrow(norm_count_))+2] <- t(baseline[,bins$chr==i])
    s <- c(1:ncol(bic_input))
    s[1:2] <- c("start", "end")
    s[2*c(1:nrow(norm_count_))+1] <-
      rownames(norm_count_)
    s[2*c(1:nrow(norm_count_))+2] <-
      paste0(rownames(norm_count_),"_ref")
    colnames(bic_input) <- s

    write.table(x=bic_input, file=paste0(tmp_dir, i, ".bin"),
                quote=F, sep="\t",
                row.names=F, col.names=T, append=F)
    system(paste0(BICseq, " -i ", paste0(tmp_dir, i, ".bin"),
                  " -l ", lambda))
    seg_re <- read.table(paste0(tmp_dir, i, ".bin_seg"),
                         header=T, sep="\t")
    bkp_ <- cumsum(seg_re$binNum)

    bkp <- merge_bkp(bkp, bkp_+k, min_interval=2)
    k <- k+sum(bins$chr==i)
  }
  for(i in unique(bins$chr)){
    file.remove(paste0(tmp_dir, i, ".bin"))
    file.remove(paste0(tmp_dir, i, ".bin_seg"))
  }
  bkp <- merge_bkp(chr_bkp, bkp, min_interval=2)
  bkp <- sort(bkp)
  print(paste0("number of breakpoints: ", length(bkp)))

  ## calculate average copy ratio of segments
  if(downsample>0){
    norm_count <- norm_count_origin
    baseline <- baseline_origin
  }
  CNV_re <- matrix(0, nrow=nrow(norm_count),ncol=length(bkp)-1)
  for(i in 1:(length(bkp)-1)){
    if(bkp[i]+1==bkp[i+1]){
      CNV_re[,i] <- norm_count[,bkp[i+1]]/baseline[,bkp[i+1]]
      next
    }
    temp <- apply(norm_count[,(bkp[i]+1):bkp[i+1]], 1, mean)/
      apply(baseline[,(bkp[i]+1):bkp[i+1]], 1, mean)
    CNV_re[,i] <- temp
  }
  CNV <- norm_count
  for(i in 1:(length(bkp)-1)){
    CNV[,(bkp[i]+1):bkp[i+1]] <-
      matrix(CNV_re[,i], nrow=nrow(CNV), ncol=bkp[i+1]-bkp[i])
  }

  ## smooth using similar cells
  if(smooth){
    print("smoothing...")
    CNV <- smooth_CNV(CNV)
  }

  ## save results
  print("Saving results...")
  CNV <- round(CNV, digit = 3)
  result <- list(copy_ratio = CNV,
                 bkp = bkp,
                 CNV_seg = CNV_re)
  saveRDS(result, file = paste0(output_dir, output_name))
  return(result)
}

#' Arm-level CNV
#'
#' @param norm_count a normalized read count matrix. Obtained from previous step
#' \code{normalize()}.
#' @param baseline read count baseline of normal cells. Obtained from previous step
#' \code{normalize()}.
#' @param genome reference genome, "hg19", "hg38" or "mm10".
#' @param min_interval minimum bin number between breakpoints.
#' @param smooth logical. Whether to aggregate similar cells to smooth final CNV
#' results?
#' @param output_dir output directory
#' @param output_name output file name, should be ".rds"
#'
#' @return Return a list of the following objects:
#' \itemize{
#' \item "copy_ratio": a cell-by-bin copy ratio matrix.
#' \item "bkp": breakpoint positions.
#' @export
#'
#' @examples
#'
calculate_CNV_arm_level <- function(norm_count,
                                    baseline,
                                    genome="hg38",
                                    min_interval=5,
                                    smooth=F,
                                    output_dir="./",
                                    output_name="CNV_result.rds"
){
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  bins <- str2bin(colnames(norm_count))
  chr_list <- unique(bins$chr)
  if(genome=="hg19"){
    data("bin_info_hg19")
    bin <- bin_info_hg19
  }else if(genome=="hg38"){
    data("bin_info_hg38")
    bin <- bin_info_hg38
  }else if(genome=="mm10"){
    data("bin_info_mm10")
    bin <- bin_info_mm10
  }else{
    print("Genome should be one of {hg19, hg38, mm10}.")
  }
  f <- bin$arm %in% c("p","q")
  chr_arm <- bin[f, ]
  chr_arm$n <- 1:nrow(chr_arm)

  CNV_re <- matrix(1, nrow=nrow(norm_count), ncol=nrow(chr_arm))
  rownames(CNV_re) <- rownames(norm_count)
  colnames(CNV_re) <- bin2str(chr_arm)

  ## calculate CNV
  bkp <- c(0)
  for(j in chr_list){
    for(arm_ in c("p","q")){
      f <- chr_arm$chr==j & chr_arm$arm==arm_
      f1 <- bins$chr==j & bins$start %in% chr_arm$start[f]
      if(sum(f1) > 0){
        bkp <- c(bkp, bkp[length(bkp)] + sum(f1))
      }
      if(sum(f1) > min_interval){
        if(sum(f1)==1){
          CNV_re[, f] <- mean(norm_count[, f1])/mean(baseline[f1])
        }else{
          temp <- apply(norm_count[, f1], 1, mean)/mean(baseline[f1])
          CNV_re[, f] <-
            matrix(data=temp, nrow=nrow(CNV_re), ncol=sum(f), byrow=F)
        }
      }
    }
  }

  ## smooth using similar cells
  if(smooth){
    print("smoothing")
    CNV_re <- smooth_CNV(CNV_re)
  }

  ## save results
  result <- list(copy_ratio = CNV_re,
                 bkp = bkp)
  saveRDS(result, file = paste0(output_dir, output_name))
  return(result)
}

#' smooth CNV
#'
#' @param CNV
#'
#' @return
#' @export
#'
#' @examples
#' @noRd
smooth_CNV <- function(CNV){
  idx <- knn.index(CNV, k=10)
  idx <- cbind(idx, 1:nrow(idx))
  temp <- function(idx){
    return(apply(CNV[idx,],2,mean))
  }
  smooth_CNV <- apply(idx, 1, temp)
  smooth_CNV <- t(smooth_CNV)
  colnames(smooth_CNV) <- colnames(CNV)
  rownames(smooth_CNV) <- rownames(CNV)
  return(smooth_CNV)
}

