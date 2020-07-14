#!/usr/bin/env Rscript

##########################################################################################################################
#    Part of the FASAS project
#    File name: Tools-SampleHeatmap.r
#    Auther: Ke Zhang
#    Version: 1.0.0
#    Data: 2019.10.10
##########################################################################################################################

##########################################################################################################################
#    Script of plotting sample heatmap
#    Copyright (C) 2019 CapitalBio Corporation
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
##########################################################################################################################

#Data Example
#Tax	Ctrl_1	Ctrl_2	Ctrl_3
#Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Blautia;Blautia wexlerae	23	54	67

library(optparse)
library(pheatmap)
library(grDevices)

my.options = list(
    make_option(c("-i", "--input_file"), type = "character", default = NULL,
                help = "Input, this is input file [default = %default]", metavar = "character"),

    make_option(c("-o", "--output_file"), type = "character", default = NULL,
                help = "Output, this is output file [default = %default]", metavar = "character"),

    make_option(c("-a", "--anno_file"), type = "character", default = NULL,
                help = "Output, this is annotation file [default = %default]", metavar = "character"),

    make_option(c("-n", "--normalization"), type = "logical", default = TRUE,
                help = "Normalization [default= %default]", metavar = "logical"),

    make_option(c("-t", "--top_number"), type = "integer", default = 100,
                help = "按平均丰度排序，取前top_number个用于heatmap [default = %default]", metavar = "integer"),

    make_option(c("-s", "--scale_type"), type = "character", default = 'none',
                help = "Heatmap的scale模式 [default = %default]", metavar = "character"),

    make_option(c("-r", "--anno_rank"), type = "integer", default = 2,
                help = "需要注释的Rank级别 [default = %default]", metavar = "integer")
)
opt <- parse_args(OptionParser(option_list = my.options))

#检查参数
if (length(opt) <= 5){
    stop("Usage: Rscript Tools-SampleHeatmap.r -h/--help")
}
if (opt$anno_rank > 7){
    error_info <- paste(sep = " ", opt$anno_rank, ": Unsupported value, supported range is 1 to 7")
    stop(error_info)
}
if (! grepl('row|column|none', opt$scale_type, perl = T)){
    error_info = paste(sep = " ", opt$scale_type, ": Must be one of none, row, or column")
}

read.table(opt$input_file, sep = "\t", header = T, row.names = 1, check.names = F) -> df

#归一化
if (opt$normalization){
    df.sum <- apply(df, 2, sum)
    df.nor <- df / df.sum
    df <- df.nor
} else {
    warning("no normalization")
}
#计算均值并取前top_number个
apply(df, 1, mean) -> df.mean
df[order(df.mean,decreasing = T),] -> df.sort
if (opt$top_number > dim(df)[1]){
    error_info <- paste(sep = ' ', inputfile, ": This file doesn't have that much classification information")
    stop(error_info)
}
df.top <- df.sort[1:opt$top_number,]
#创建注释信息
taxon <- data.frame(taxon = row.names(df.top), stringsAsFactors = F)
#创建用于split的function
my.split <- function (x) { strsplit(x, split = ';') }
#对每一行进行split，每一行都返回一个list
taxon.split <- apply(taxon, 1, my.split)
#顺序追加信息
anno <- vector()
for (i in 1:length(taxon.split)) {
    anno <- append(anno, taxon.split[[i]]$taxon[opt$anno_rank])
}
anno <- data.frame(row.names = row.names(df.top), Anno_rank = anno)
#获取注释颜色
newCols <- colorRampPalette(rainbow(length(unique(anno[,1]))))
mycolors <- newCols(length(unique(anno[,1])))
names(mycolors) <- unique(anno[,1])
mycolors <- list(Anno_rank = mycolors)
#重新选择注释级别的名称
rank_level <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
names(mycolors) <- rank_level[opt$anno_rank]
colnames(anno) <- rank_level[opt$anno_rank]
#作图
p <- pheatmap(df.top, cellheight = 5, cellwidth = 12,
         cluster_rows = T, show_rownames = F,
         file = opt$output_file,
         scale = opt$scale_type,
         annotation_row = anno, annotation_colors = mycolors)
#导出注释表格
#从p中直接得到顺序，然后df.top[row.names]得到注释表格
df.anno <- df.top[p$tree_row$labels[p$tree_row$order],]
write.table(df.anno, sep = "\t", file = opt$anno_file, quote = F, row.names = T, col.names = NA)
