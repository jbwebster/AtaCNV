#' scATAC-seq count simulation using R package simATAC
#'
#' @param ref_count
#' @param ncell
#' @param gaussmean
#' @param gausssd
#' @param CN
#' @param coverage
#' @param random_seed
#'
#' @return
#' @export
#'
#' @examples
#' @noRd
sim_count <- function(ref_count, ncell, gaussmean,
                      gausssd, CN, coverage=2, random_seed=1234){
  # set.seed(random_seed)
  ATACobj = simATAC::newsimATACCount()
  ATACobj = simATAC::simATACEstimate(BiocGenerics::t(ref_count))
  ## Add Gaussian Noise
  ATACobj = simATAC::setParameters(ATACobj, noise.mean=gaussmean, noise.sd=gausssd)
  ## Set factor
  ATACobj = simATAC::setParameters(ATACobj, sparse.fac=coverage*CN/2)
  ## Set the cell number
  sim = simATAC::simATACSimulate(ATACobj, nCells=ncell)
  simcount = BiocGenerics::counts(sim)
  return(t(as.matrix(simcount)))
}

#' Random choose CNV positions
#'
#' @param nbin
#' @param ncell
#' @param ngroup
#' @param nCNV
#' @param chr_arm
#' @param CNV_candidates
#' @param total_group
#'
#' @return
#' @export
#'
#' @examples
#' @noRd
sim_true_CN <- function(nbin, ncell, ngroup, nCNV, chr_arm,
                        CNV_candidates, total_group=5){
  true_CN <- matrix(data=2, nrow=ncell, ncol=nbin)
  if(ngroup==0 | nCNV==0){
    return(list(true_CN=true_CN,CNV_list="none"))
  }
  k <- floor(ncell/total_group)
  CNV_list <- matrix(data=NA, nrow=ngroup, ncol=nCNV)
  for(i in 1:ngroup){
    CNV_arm <- sample(1:nrow(chr_arm), nCNV)
    for(j in 1:nCNV){
      temp1 <- sample(CNV_candidates,1)
      temp2 <- c(chr_arm$start[CNV_arm[j]]:chr_arm$end[CNV_arm[j]])
      true_CN[((i-1)*k+1):(i*k),temp2] <- temp1
      CNV_list[i,j] <- paste0(chr_arm$chr[CNV_arm[j]], " ", temp1)
    }
  }
  return(list(true_CN=true_CN,CNV_list=CNV_list))
}

#' Simulate data with given CNV
#'
#' @param ref_count
#' @param ncell
#' @param gaussmean
#' @param gausssd
#' @param CN
#' @param CNV_candidates
#' @param coverage
#' @param random_seed
#'
#' @return
#' @export
#'
#' @examples
#' @noRd
sim_count_given_CNV <- function(ref_count,
                                ncell,
                                gaussmean,
                                gausssd, CN, CNV_candidates,
                                coverage,
                                random_seed){
  simcount <- CN
  CN_states <- union(2, CNV_candidates)
  for(i in CN_states){
    if(sum(CN==i)==0){
      next
    }
    temp <-
      sim_count(ref_count=ref_count, ncell=ncell,
                gaussmean=gaussmean, gausssd=gausssd,
                CN=i, coverage=coverage, random_seed=random_seed)
    simcount[CN==i] <- temp[CN==i]
    print(dim(temp))
  }
  return(simcount)
}

#' Estimate copy ratio distribution for given copy number
#'
#' @param count
#' @param genome
#' @param bin
#' @param n_cnv
#' @param n_state
#'
#' @return
#' @export
#'
#' @examples
#' @noRd
estimate_cnv_state_paramater <- function(count, genome,
                                         bin,
                                         n_cnv=10, n_state=3,
                                         random_seed = 1234){
  # count: cell * bin
  set.seed(random_seed)
  count <- as(count, "dgCMatrix")
  f <- bin$arm %in% c("p", "q")
  chr_arm <- data.frame(chr=unique(paste0(bin$chr, bin$arm)[f]),
                        start=0, end=0)
  chr_arm$map_bin <- 0
  bin$n <- 1:nrow(bin)
  for(i in 1:nrow(chr_arm)){
    f <- paste0(bin$chr, bin$arm) == chr_arm$chr[i]
    bin_ <- bin[f,]
    chr_arm$start[i] <- min(bin_$n)
    chr_arm$end[i] <- max(bin_$n)
    chr_arm$map_bin[i] <- sum(bin_$map > 500000)
  }
  chr_arm$len <- chr_arm$end - chr_arm$start + 1
  chr_arm <- chr_arm[chr_arm$map_bin>10,]

  ncell <- 500
  # print(n_state)
  if(n_state == 3){
    CNV_candidates <- c(1,3)
  }else if(n_state == 5){
    CNV_candidates <- c(0,1,3,4)
  }else{
    print("Number of CNV state should be 3 or 5")
  }

  temp <- sim_true_CN(nbin=ncol(count), ncell=ncell,
                      ngroup=3, nCNV=n_cnv, chr_arm=chr_arm,
                      CNV_candidates=CNV_candidates)
  true_CN <- temp$true_CN
  simcount <-
    sim_count_given_CNV(ref_count=count, ncell=ncell,
                        gaussmean=0, gausssd=0,
                        CN=true_CN, CNV_candidates=CNV_candidates,
                        coverage=1, random_seed=1234)
  colnames(simcount) <- bin2str(bin)
  rownames(simcount) <- paste0("Cell", c(1:ncell))
  # cell_barcodes <- rownames(simcount)
  colnames(true_CN) <- colnames(simcount)
  rownames(true_CN) <- rownames(simcount)

  simcount <- t(simcount)
  temp <- colnames(simcount)[300:500]
  norm_re <- AtaCNV::normalize(count = t(simcount),
                               genome = genome,
                               mode = "normal cells",
                               normal_cells = temp,
                               cutoff = 3)

  copy_ratio <- norm_re$copy_ratio
  true_CN_filtered <- true_CN[ , bin2str(bin) %in% colnames(copy_ratio)]

  x <- matrix(nrow=5000, ncol=n_state)
  for(i in 1:n_state){
    if(n_state == 3){
      x[,i] <- sample(as.vector(copy_ratio[true_CN_filtered==i]), 5000)
    }
    if(n_state == 5){
      x[,i] <- sample(as.vector(copy_ratio[true_CN_filtered==(i-1)]), 5000)
    }
  }

  mu <- apply(x, 2, mean)
  sigma <- apply(x, 2, sd)
  return(list(mu, sigma, x))
}

#' Estimate copy number
#'
#' @param count the raw cell-by-bin read count matrix.
#' @param genome reference genome, "hg19", "hg38" or "mm10".
#' @param copy_ratio the estimated copy ratio matrix from \code{normalize()}.
#' @param bkp the detected breakpoint positions from \code{calculate_CNV()}.
#' @param label a vector indicating normal cell and tumor subclone. Its length
#' should be the same as the row number of count matrix, and normal cell
#' should be labeled as "normal".
#'
#' @return Return a list of the following objects:
#' \itemize{
#' \item "CN_state": a cell-by-bin copy number state matrix, where elements
#' are one of 0.5 (copy number loss), 1 (copy number neural) and 1.5 (copy
#' number gain).
#' \item "copy_ratio_adjusted": a cell-by-bin copy ratio matrix.
#' }
#' @export
#'
#' @examples
#'
estimate_cnv_state_cluster <- function(count,
                                       genome,
                                       copy_ratio,
                                       bkp,
                                       label,
                                       n_iter = 2){
  f <- label == "normal"
  if(sum(f) <= 1){
    print("No normal cells in label!")
  }
  if(genome=="hg19"){
    data("bin_info_hg19")
    bin <- bin_info_hg19
  }else if(genome=="hg38"){
    data("bin_info_hg38")
    bin <- bin_info_hg38
  }else if(genome=="mm10"){
    data("bin_info_mm10")
    bin <- bin_info_mm10
  }
  n_state <- 3
  print("Estimating copy ratio distribution...")
  tmp <- estimate_cnv_state_paramater(count = count[f,],
                                      genome = genome,
                                      bin = bin,
                                      n_state = n_state)
  mu <- tmp[[1]]
  sigma <- tmp[[2]]

  ## initial estimate
  print("Assigning initial copy number state...")
  CNV_state <- copy_ratio
  CNV_state[,] <- 1
  prob_mat <- CNV_state

  for(k in unique(label)){
    if(k == "normal"){
      next
    }
    f <- label == k
    for(i in 1:(length(bkp) - 1)){
      seg <- (bkp[i]+1):bkp[i+1]
      x <- as.vector(copy_ratio[f, seg])
      prob <- rep(0, length(mu))
      for(i in 1:length(prob)){
        tmp <- log(dnorm(x, mean=mu[i], sd=sigma[i]) + 1e-6)
        prob[i] <- mean(tmp)
      }
      prob <- exp(prob)/sum(exp(prob))
      CNV_state[f, seg] <- which.max(prob) / 2
      prob_mat[f, seg] <- max(prob)
    }
  }
  CNV_state_init <- CNV_state

  ## adjust and re-estimate
  for(iter in 1:n_iter){
    f <- colMeans(CNV_state) == 1
    print(paste0("Number of diploid region: ", sum(f)))
    print("Adjusting result...")
    adjust_ratio <- rowMeans(copy_ratio[,f])
    copy_ratio_adjusted <- copy_ratio / adjust_ratio

    CNV_state <- copy_ratio_adjusted
    CNV_state[,] <- 1
    prob_mat <- CNV_state


    for(k in unique(label)){
      if(k == "normal"){
        next
      }
      f <- label == k
      for(i in 1:(length(bkp) - 1)){
        seg <- (bkp[i]+1):bkp[i+1]
        x <- as.vector(copy_ratio_adjusted[f, seg])
        prob <- rep(0, length(mu))
        for(i in 1:length(prob)){
          tmp <- log(dnorm(x, mean=mu[i], sd=sigma[i]) + 1e-6)
          prob[i] <- mean(tmp)
        }
        prob <- exp(prob)/sum(exp(prob))
        CNV_state[f, seg] <- which.max(prob) / 2
        prob_mat[f, seg] <- max(prob)
      }
    }
  }

  CNV_state_final <- CNV_state
  return(list(CN_state = CNV_state_final,
              copy_ratio_adjusted = copy_ratio_adjusted))
}
