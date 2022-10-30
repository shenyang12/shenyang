
rm(list=ls())
setwd("C:\\Users\\admin\\Desktop\\reAnalysisTP53\\08.mutation_maf")

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")


library(maftools)

#TCGA-pancancer MAF file (gz)----download https://gdc.cancer.gov/about-data/publications/pancanatlas

maf= read.maf(maf = "mc3.v0.2.8.PUBLIC.maf.gz")


#Shows sample summry.
getSampleSummary(maf)
#Shows gene summary.
getGeneSummary(maf)
#Shows all fields in MAF
getFields(maf)

#查看maf基本信息
slotNames(maf)

table(maf@data$Variant_Type)

table(maf@data$Variant_Classification)

table(maf@data$Mutation_Status)

table(maf@data$Variant_Classification)
#Writes maf summary to an output file with basename maf.
write.mafSummary(maf = maf, basename = 'maf')

##提取maf子集
maf_TP53=subsetMaf(maf = maf, genes = c('TP53','TP63','TP73'), mafObj = TRUE)

slotNames(maf_TP53)
variant_class=table(maf_TP53@data$Variant_Classification)
write.table(variant_class,"variant_class.txt",sep="\t")
variant_Type=table(maf_TP53@data$Variant_Type)
write.table(variant_Type,"variant_Type.txt",sep="\t")

#pl#plotmafSummary
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = maf, top = 20)
oncoplot(maf = maf, genes = c("TP53",'TP63', 'TP73'))

pdf(file="oncostrip.pdf", width=20, height=10)
oncostrip(maf = maf, genes = c("TP53",'TP63', 'TP73'))
dev.off()


# 修正TCGA名称
rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # 数据来自PanCanAtlas
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
write.table(samAnno,"simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)

unique(maf@data$Variant_Classification)
# 设置coding区域的突变以及非coding区域的突变
## coding区域（个人认为是非沉默突变）
vc_nonSyn <- c("Frame_Shift_Del",
               "Frame_Shift_Ins",
               "Splice_Site",
               "Translation_Start_Site",
               "Nonsense_Mutation",
               "Nonstop_Mutation",
               "In_Frame_Del",
               "In_Frame_Ins",
               "Missense_Mutation")
## 非coding区域（个人认为是沉默突变）
vc_Syn <- c("3'UTR",
            "5'UTR",
            "3'Flank",
            "5'Flank",
            "IGR",
            "Intron",
            "Silent")

# 获取不同癌种的样本数
mafsam <- data.frame(samID = maf@data$Tumor_Sample_Barcode,
                     simple_barcode = substr(maf@data$Tumor_Sample_Barcode,1,15),
                     stringsAsFactors = F)
mafsam <- mafsam[!duplicated(mafsam$simple_barcode),]
mafsam <- merge(mafsam,samAnno, by = "simple_barcode", all.x = TRUE)
mafsam <- as.data.frame(na.omit(mafsam))

txt <- data.frame(label = paste0(names(table(mafsam$`cancer type`))," (n=", as.numeric(table(mafsam$`cancer type`)),")"),
                  "cancer type" = names(table(mafsam$`cancer type`)),
                  check.names = F) # 将肿瘤和样本数标签合并
mafsam <- merge(mafsam, txt, by = "cancer type", all.x = TRUE)
rownames(mafsam) <- mafsam$simple_barcode



maf@data$Tumor_Sample_Barcode <- substr(maf@data$Tumor_Sample_Barcode,1,15)
mafAnno <- mafsam[which(mafsam$simple_barcode %in% maf@data$Tumor_Sample_Barcode),]
colnames(mafAnno)[1:2] <- c("Cancer","Tumor_Sample_Barcode")


mutcol <- brewer.pal(n = 10, name = 'Paired')
names(mutcol) <- c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')

cancercolors <- NMF:::ccRamp(brewer.pal(n = 12,name = 'Paired'), length(unique(mafAnno$Cancer)))
names(cancercolors) <- unique(mafAnno$Cancer)

annocolors = list(Cancer = cancercolors)
maf@data$Tumor_Sample_Barcode <- substr(maf@data$Tumor_Sample_Barcode,1,15)


pdf(file = "oncoprint.pdf", width = 10, height = 5)
oncoplot(maf = maf, # MAF文件
         colors = mutcol, # 突变类型的颜色
         genes = c("TP53","TP63","TP73"),
         bgCol = "grey95", # 瀑布图背景色
         borderCol = NA, # 突变的边框（由于样本太多所以不设置边框）
         annotationDat = mafAnno, # 样本注释文件
         annotationColor = annocolors,
         clinicalFeatures = "Cancer", # 从注释文件中显示的信息
         sortByAnnotation = T) # 按照mafAnno的第一列排序，即按照Cancer排序
invisible(dev.off())



