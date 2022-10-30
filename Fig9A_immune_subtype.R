
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)
library(reshape2)

expFile="panGeneExp.txt"                         #表达输入文件
subtypeFile="Subtype_Immune_Model_Based.txt"     #免疫分组文件
setwd("E:\\desktop\\TP53 family\\11_immune Type")     #设置工作目录

#读取表达文件，并对输入文件整理
exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]
exp=as.matrix(exp[,1:(ncol(exp)-2)])
row.names(exp)=gsub(".$","",row.names(exp))
exp=avereps(exp)

#读取免疫分组文件
subtype=read.table(subtypeFile, header=T,sep="\t",row.names=1,check.names=F)

#样品取交集
sameSample=intersect(row.names(subtype),row.names(exp))
subtype=subtype[sameSample,]
subtype=gsub(".+Immune |\\)","",subtype)
exp=exp[sameSample,]
exp=cbind(as.data.frame(exp),subtype)

#差异分析
outTab=data.frame()
geneSig=c()
for(gene in colnames(exp)[1:(ncol(exp)-1)]){
	rt1=exp[,c(gene,"subtype")]
	colnames(rt1)=c("expression","subtype")
	ks=kruskal.test(expression ~ subtype, data = rt1)
	p=ks$p.value
	outTab=rbind(outTab,cbind(gene,pvalue=p))
	Sig=ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","")))
	geneSig=c(geneSig,Sig)
}
geneSig=c(geneSig,"")
colnames(exp)=paste0(colnames(exp),geneSig)
write.table(outTab,file="immuneType.xls",sep="\t",row.names=F,quote=F)

#把数据转换成ggplot2输入文件
data=melt(exp)
colnames(data)=c("Subtype","Gene","Expression")

#绘制图形
p1=ggplot(data,aes(x=Subtype,
                y=Expression,
                fill=Subtype))+
    guides(fill=guide_legend(title="Immune Subtype"))+
    labs(x = "Immune Subtype", y = "Gene expression")+
	geom_boxplot()+ facet_wrap(~Gene,nrow =1)+ theme_bw()
pdf(file="immuneType.pdf",width=9,height=5)
print(p1)
dev.off()
