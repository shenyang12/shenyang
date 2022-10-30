
#install.packages("corrplot")


library(corrplot)                    #引用包
expFile="panGeneExp.txt"             #表达输入文件
scoreFile="estimateScores.txt"       #免疫微环境结果文件
scoreType="StromalScore"             #打分类型
setwd("E:\\desktop\\TP53 family\\13_estimateCor")     #设置工作目录

#读取表达文件
exp=read.table(expFile, header=T,sep="\t",row.names=1,check.names=F)
exp=exp[(exp[,"Type"]=="Tumor"),]

#读取肿瘤微环境文件
TME=read.table(scoreFile, header=T,sep="\t",row.names=1,check.names=F)

#样品取交集
sameSample=intersect(row.names(TME),row.names(exp))
TME=TME[sameSample,]
exp=exp[sameSample,]

#相关性检验
outTab=data.frame()
pTab=data.frame()
#按肿瘤类型循环
for(i in levels(factor(exp[,"CancerType"]))){
    exp1=exp[(exp[,"CancerType"]==i),]
    TME1=TME[(TME[,"CancerType"]==i),]
    x=as.numeric(TME1[,scoreType])
    pVector=data.frame(i)
    outVector=data.frame(i)
	#按基因循环
	genes=colnames(exp1)[1:(ncol(exp1)-2)]
	for(j in genes){
		y=as.numeric(exp1[,j])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pValue=corT$p.value
		pVector=cbind(pVector,pValue)
		outVector=cbind(outVector,cor)
	}
	pTab=rbind(pTab,pVector)
	outTab=rbind(outTab,outVector)
}
#输出相关性的表格
colNames=c("CancerType",colnames(exp1)[1:(ncol(exp1)-2)])
colnames(outTab)=colNames
write.table(outTab,file="estimateCor.cor.txt",sep="\t",row.names=F,quote=F)
#输出相关性检验p值的表格
colnames(pTab)=colNames
write.table(pTab,file="estimateCor.pval.txt",sep="\t",row.names=F,quote=F)

#绘制微环境相关性图形
pdf(file=paste0(scoreType,".pdf"),height=3.5,width=20)
row.names(outTab)=outTab[,1]
outTab=outTab[,-1]
corrplot(corr=as.matrix(t(outTab)),
		 title=paste0("\n\n\n\n",scoreType),
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
