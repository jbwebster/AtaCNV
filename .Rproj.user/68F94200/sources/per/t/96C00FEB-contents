#' Merge two vector of breakpoints together
#'
#' @param bkp1 breakpoint vector 1
#' @param bkp2 breakpoint vector 2
#' @param min_interval minimum interval between two breakpoint
#'
#' @return merged breakpoint vector
#' @export
#'
#' @examples
#' @noRd
merge_bkp <- function(bkp1, bkp2, min_interval=3){
  # bkp2 all reserved
  for(i in bkp1){
    if(min(abs(bkp2-i)>min_interval)){
      bkp2[length(bkp2)+1] <- i
    }
  }
  return(bkp2)
}

#' Strings to bin data frame
#'
#' @param bin_name a string vector. Elements are like "chr1_0_1000000".
#'
#' @return a data frame containing "chr", "start", and "end".
#' @export
#'
#' @examples
#' @noRd
str2bin <- function(bin_name){
  temp <- strsplit(bin_name, split="_")
  temp2 <- unlist(temp)
  dim(temp2) <- c(3,length(temp2)/3)
  bins <- data.frame(chr=temp2[1,],start=as.numeric(temp2[2,]),
                     end=as.numeric(temp2[3,]))
  bins$n <- c(1:nrow(bins))
  return(bins)
}

#' Bin data frame to string
#'
#' @param bin a data frame containing "chr", "start", and "end".
#'
#' @return a string vector,
#' @export
#'
#' @examples
#' @noRd
bin2str <- function(bin){
  return(paste0(bin$chr,"_",bin$start,"_",bin$end))
}

#' Obtain breakpoints between chromosomes from bin data frame
#'
#' @param bins a data frame containing "chr", "start", and "end".
#'
#' @return breakpoints between chromosomes
#' @export
#'
#' @examples
#' @noRd
bin2chrbkp <- function(bins){
  chr_bkp <- c(0)
  temp <- unique(bins$chr)
  for(i in 1:length(temp)){
    chr_bkp[i+1] <- max(bins$n[bins$chr==temp[i]])
  }
  return(chr_bkp)
}

#' Delete specific chromosomes
#'
#' @param count a cell-by-bin matrix with formatted column names
#' @param chr_list a list of chromosomes to remove
#'
#' @return a cell-by-bin matrix
#' @export
#'
#' @examples
del_chr <- function(count, chr_list=c("chrX","chrY")){
  bin <- str2bin(colnames(count))
  return(count[,!bin$chr %in% chr_list])
}
