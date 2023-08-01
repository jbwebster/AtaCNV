#' Generate input matrix from fragment file for AtaCNV pipeline
#'
#' @param fragment_file path of fragment file.
#' @param cell_barcodes a string vector giving barcodes that represent cells in
#' scATAC-seq expriment.
#' @param genome reference genome, "hg19" or "hg38".
#' @param output_dir output directory
#'
#' @return a cell-by-bin read count matrix
#' @export
#'
#' @examples
#'
generate_input <- function(fragment_file, cell_barcodes, genome="hg38",
                           output_dir="./"){
  if(genome=="hg19"){
    data("chr_length_hg19")
    chr_info <- chr_length_hg19
  }else if(genome=="hg38"){
    data("chr_length_hg38")
    chr_info <- chr_length_hg38
  }
  print(paste0("Reference genome: ", genome))
  chr_list <- chr_info$chrom
  chr_length <- chr_info$size
  barcode <- cell_barcodes
  K <- 1000000

  start <- c()
  end <- c()
  chr <- c()
  for(i in 1:length(chr_list)){
    n <- floor(chr_length[i]/K)
    temp <- K*c(0:n)
    temp <- c(temp, chr_length[i])
    start <- c(start, temp[1:(n+1)])
    end <- c(end, temp[2:(n+2)])
    chr <- c(chr, rep(chr_list[i], n+1))
  }
  bin <- data.frame(chr, start, end)

  count <- matrix(data=0, nrow=length(barcode), ncol=length(start))
  row.names(count) <- barcode

  t1 <- proc.time()
  # con <- file(fragment_file, "r")
  con <- gzcon(file(fragment_file,open="r"))
  flag <- 0
  nchr <- 1
  nchr_start <- 0
  while(TRUE){
    flag <- flag+1
    line <- readLines(con, n=1)
    # print(line)
    if(length(line)==0){
      break
      print("all finished")
    }
    temp <- strsplit(line, "\t")[[1]]
    cell_barcode <- temp[4]
    if(!(cell_barcode %in% barcode) | !(temp[1]) %in% chr_list){
      next
    }
    if(temp[1]!=chr_list[nchr]){
      print(paste(chr_list[nchr], "finished"))
      nchr_start <- nchr_start+floor(chr_length[nchr]/K)+1
      nchr <- nchr+1
      if(nchr > length(chr_list)){
        break
      }
    }
    pos <- as.numeric(temp[2])
    # count[cell_barcode, nchr_start+floor(pos/K)+1] <-
    #   count[cell_barcode, nchr_start+floor(pos/K)+1]+1
    count[cell_barcode, nchr_start+floor(pos/K)+1] <-
      count[cell_barcode, nchr_start+floor(pos/K)+1]+as.numeric(temp[5])
  }
  close(con)
  t2 <- proc.time()
  print(t2-t1)

  print(paste0("Saving to ", output_dir))
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(count,file=paste0(output_dir,"input.rds"))
}
