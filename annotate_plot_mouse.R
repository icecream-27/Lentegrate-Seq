# 安装必要的包（如果尚未安装）
install.packages("ggplot2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")

# 加载包
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(scales)
library(tools)

# 加载TxDb数据库
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene

# 定义处理单个BED文件的函数
process_bed_file <- function(bed_file, output_dir) {
  file_base_name <- file_path_sans_ext(basename(bed_file))
  binding_sites <- import(bed_file, format = "BED", extraCols = c(reads = "integer"))
  
  # 将导入的GRanges对象转换为数据框
  binding_sites_df <- as.data.frame(binding_sites)
  
  # 确保reads列被正确读取
  print(head(binding_sites_df))
  
  # 使用annotatePeak函数进行注释, 从TSS向上游3000bp到TSS向下游3000bp的范围内进行注释
  annotated_sites <- annotatePeak(binding_sites, TxDb = txdb, tssRegion=c(-3000, 3000), verbose=TRUE)
  
  # 将注释结果转换为数据框
  annotated_sites_df <- as.data.frame(annotated_sites)
  
  # 合并reads数信息到注释结果中
  annotated_sites_df$reads_count <- binding_sites_df$reads
  
  # 打印结果，查看是否包含reads数信息
  print(head(annotated_sites_df))
  
  # 提取注释类别，并使用reads数作为计数
  annotation_counts <- tapply(annotated_sites_df$reads_count, annotated_sites_df$annotation, sum)
  
  # 将注释类别计数转换为数据框
  annotation_df <- as.data.frame(annotation_counts)
  annotation_df$Annotation <- rownames(annotation_df)
  colnames(annotation_df) <- c("Count", "Annotation")
  
  # 汇总注释类别
  summarized_df <- annotation_df %>%
    mutate(Annotation = case_when(
      grepl("Exon", Annotation) ~ "Exon",
      grepl("Intron", Annotation) ~ "Intron",
      grepl("Intergenic", Annotation) ~ "Intergenic",
      grepl("Promoter", Annotation) ~ "Promoter",
      grepl("UTR", Annotation) ~ "UTR",
      TRUE ~ "Others"
    )) %>%
    group_by(Annotation) %>%
    summarise(Count = sum(Count)) %>%
    mutate(Percentage = floor(Count / sum(Count) * 1000) / 10) %>%  # 直接截断小数位
    mutate(Label = ifelse(Percentage < 1, "<1%", paste0(Percentage, "%")))
  
  # 设置科研配色
  colors <- c("Exon" = "#D55E00", "Intron" = "#542788", "Intergenic" = "#0072B2", 
              "Promoter" = "#CCB974", "UTR" = "#008837", "Others" = "#262626")
  
  # 计算总数
  total_count <- sum(summarized_df$Count)
  
  # 自定义标签
  summarized_df$Label <- paste0(summarized_df$Annotation, " (", summarized_df$Label, ")")
  
  # 绘制饼图并添加总数标题
  p <- ggplot(summarized_df, aes(x = "", y = Count, fill = Annotation)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors, labels = summarized_df$Label) +
    labs(title = paste("Annotation Distribution", "\nTotal =", total_count), x = "", y = "") +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10))
  
  # 创建输出文件夹
  png_dir <- file.path(output_dir, "png")
  annotated_dir <- file.path(output_dir, "annotated")
  if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)
  if (!dir.exists(annotated_dir)) dir.create(annotated_dir, recursive = TRUE)
  
  # 保存图片
  output_file <- file.path(png_dir, paste0(file_base_name, ".png"))
  ggsave(output_file, plot = p, width = 8, height = 6)
  
  # 将 annotated_sites_df 写入文件
  write.table(annotated_sites_df, file = file.path(annotated_dir, paste0(file_base_name, "-annotated.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 将 annotation_counts 转换为数据框并写入文件
  write.table(annotation_df, file = file.path(annotated_dir, paste0(file_base_name,"-annotation_counts.txt")), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 打印注释类别的计数
  print(annotation_df)
}

# 定义主函数
process_all_bed_files <- function(input_dir) {
  bed_files <- list.files(input_dir, pattern = "\\.bed$", full.names = TRUE)
  output_dir <- input_dir
  
  for (bed_file in bed_files) {
    process_bed_file(bed_file, output_dir)
  }
}

# 执行主函数，输入存储BED文件的文件夹路径
input_directory <- "P428-bed"  # 修改为实际的文件夹路径
process_all_bed_files(input_directory)

