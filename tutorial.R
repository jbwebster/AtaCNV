library("AtaCNV")

cell_info <- readRDS("./example/cell_info.rds")
count <- readRDS("./example/count.rds")
count_paired <- readRDS("./example/count_paired.rds")

## normalization
# option1: mode="matched normal sample"
norm_re1 <- AtaCNV::normalize(count,
                              genome="hg19", mode="matched normal sample",
                              cell_cluster=cell_info$cluster,
                              count_paired=count_paired,
                              output_dir="./",
                              output_name="norm_re2.rds"
)

# option2: mode="normal cells"
norm_re2 <- AtaCNV::normalize(count,
                              genome="hg19", mode="normal cells",
                              normal_cells=(cell_info$cell_type=="normal"),
                              output_dir="./",
                              output_name="norm_re2.rds"
)

# option3: mode="all cells"
norm_re3 <- AtaCNV::normalize(count,
                              genome="hg19", mode="all cells",
                              cell_cluster=cell_info$cluster,
                              output_dir="./",
                              output_name="norm_re3.rds"
)

## option4: mode="none"
norm_re4 <- AtaCNV::normalize(count,
                              genome="hg19", mode="none",
                              cell_cluster=cell_info$cluster,
                              output_dir="./",
                              output_name="norm_re4.rds"
)

# visualize copy ratio after normalization
plot_heatmap(copy_ratio=norm_re1$copy_ratio, 
             cell_cluster=cell_info$cluster,
             output_dir="./",
             output_name="copy_ratio1.png")

## segmentation
CNV_re <- calculate_CNV(norm_count=norm_re1$norm_count,
                        baseline=norm_re1$baseline,
                        output_dir="./",
                        output_name="CNV_re1.rds")

# visulize copy ratio after segmentation
plot_heatmap(copy_ratio=CNV_re$CNV, 
             cell_cluster=cell_info$cluster,
             output_dir="./",
             output_name="CNV1.png")
