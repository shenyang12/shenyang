
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#引用包
library(limma)

gene="TP53"            #基因名称，改成自己研究的基因
setwd("C:\\Users\\admin\\Desktop\\TP53_Family \\16_geneCor")      #设置工作目录
geneRT=read.table("gene.txt",sep="\t",header=F)                 #读取基因表达文件
files=dir()
files=grep("^symbol.",files,value=T) ##来源于UCSC xena下载

#按肿瘤类型循环进行相关性检验
outTab=data.frame()
corTab=data.frame()
sameGenes=c()
for(i in files){
	#读取表达文件
	CancerType=gsub("^symbol\\.|\\.txt$","",i)
	rt=read.table(i, header=T,sep="\t",check.names=F)
	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	
	#去除正常样品
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
	#提取gene.txt里面基因的表达量
	sameGenes=intersect(as.vector(geneRT[,1]),row.names(data))
	data=data[c(sameGenes,gene),]
	
	#对基因进行循环
	x=as.numeric(data[gene,])
    outVector=data.frame(CancerType)
    corVector=data.frame(CancerType)
	for(j in sameGenes){
		y=as.numeric(data[j,])
		if(sd(y)>0.01){
			corT=cor.test(x,y)
			cor=corT$estimate
			pValue=corT$p.value
			outVector=cbind(outVector,pValue)
			corVector=cbind(corVector,cor)
		}
		else{
			outVector=cbind(outVector,pValue=1)
			corVector=cbind(corVector,cor=0)
		}
	}
	outTab=rbind(outTab,outVector)
	corTab=rbind(corTab,corVector)
}
#输出相关性检验p值结果
colnames(outTab)=c("CancerType",sameGenes)
write.table(outTab,file="geneCor.pvalue.txt",sep="\t",row.names=F,quote=F)
#输出相关性结果
colnames(corTab)=c("CancerType",sameGenes)
write.table(corTab,file="geneCor.cor.txt",sep="\t",row.names=F,quote=F)
