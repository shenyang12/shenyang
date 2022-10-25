
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)             #引用包
geneFile="gene.txt"        #基因列表文件
setwd("E:\\desktop\\TP53 family\\03_panGeneExp/")    #设置工作目录

#读取基因文件
gene=read.table("gene.txt",sep="\t",header=F,check.names=F)
genelist=as.vector(gene[,1])
genelist=gsub(" ","",genelist)

#读取目录下的文件
files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
	#读取文件
	CancerType=gsub("symbol\\.|\\.txt","",i)
	rt=read.table(i,sep="\t",header=T,check.names=F)

	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)

	#得到样品的Type
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	Type=ifelse(group==0,"Tumor","Normal")
	geneExp=t(data[genelist,])
	#geneExp=geneExp/data["TBP",]       #是否按照持家基因TBP对数据矫正，需要矫正去掉前面#号
	outTab=rbind(outTab,cbind(geneExp,Type,CancerType))
}

#输出结果表格
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="panGeneExp.txt",sep="\t",quote=F,row.names=F)

