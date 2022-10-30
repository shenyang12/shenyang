
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

library(limma)
library(ggpubr)
rm(list=ls())
setwd("E:\\桌面临时文件\\TP53_Family _2\\06_P53_clicor")                     #修改工作目录
file="expTime.txt"                                                      #输入文件
rt=read.table(file,sep="\t",header=T,check.names=F) 
rt=rt[,c(1,3,5,6)]
write.table(rt,file="p63Exp.txt",sep="\t",quote=F,row.names = F)

#读取表达数据文件
cli=read.table("clinical.txt",sep="\t",header=T,check.names=F) 
dim(cli)
cli[duplicated(cli),]  #查看cli重复行
# Remove duplicates based on sample rows
cli1 <- cli[!duplicated(cli$sample), ]
dim(cli1)

write.table(cli1,file="cli1.txt",sep="\t",quote=F,row.names = F)


#读取临床数据文件

file="p53Exp.txt"                                                      #输入文件
rt=read.table(file,sep="\t",header=T,check.names=F,row.names=1)              #读取表达数据文件
cli=read.table("cli1.txt",sep="\t",header=T,check.names=F,row.names=1)    #读取临床数据文件
gene=colnames(rt)[1]
clinical=colnames(cli)[1]


#对肿瘤类型循环，进行临床相关性分析
outTab=data.frame()
for(i in levels(factor(rt[,"CancerType"]))){
  rt1=rt[(rt[,"CancerType"]==i),]
  #交集样品
  data=cbind(rt1,gene=rt1[,gene])
  data=as.matrix(data[,c(gene,"gene")])
  if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
    row.names(data)=gsub(".$","",row.names(data))}
  data=avereps(data)
  sameSample=intersect(row.names(data),row.names(cli))
  sameData=data[sameSample,]
  sameClinical=cli[sameSample,]
  library(stringr)
   cliExpData=cbind(as.data.frame(sameClinical),sameData)
   cliExpData=cliExpData[order(cliExpData$sameClinical),]
   cliExpData=as.data.frame(cliExpData)
  #sameClinical=str_sort(sameClinical,decreasing = FALSE)
  #sameData=str_sort(sameData,decreasing=FALSE)
  if(nrow(cliExpData)==0){next}
  #设置比较组
  group=levels(factor(cliExpData$sameClinical))
  comp=combn(group,2)
  my_comparisons=list()
  for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
  #绘制boxplot
  boxplot=ggboxplot(cliExpData, x="sameClinical", y="gene", color="sameClinical",
                    xlab=clinical,
                    ylab=paste(gene,"expression"),
                    legend.title=clinical,
                    title=paste0("Cancer: ",i),
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  pdf(file=paste0(clinical,".",i,".pdf"),width=5.5,height=5)
  print(boxplot)
  dev.off()

}

