
#install.packages("ggpubr")

library(ggpubr)               #引用包
expFile="panGeneExp.txt"      #表达输入文件
setwd("E:\\desktop\\TP53 family\\05_diff")                    #设置工作目录
data=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)     #读取箱线图输入文件

#统计每个肿瘤正常样品数目,挑选正常样品数目>=5的肿瘤
Normal=data[data[,"Type"]=="Normal",]
NormalNum=table(Normal[,"CancerType"])
NormalNum=NormalNum[NormalNum>=5]
NormalCacner=names(NormalNum)
data=data[which(data[,"CancerType"] %in% NormalCacner),]

#对基因进行循环，并做差异分析
genelist=colnames(data[,(1:(ncol(data)-2))])
for(gene in genelist){
	rt1=data[,c(gene,"Type","CancerType")]
	colnames(rt1)[1]="expression"
	#绘制箱型图
	p=ggboxplot(rt1, x="CancerType", y="expression", color = "Type",
	     ylab=paste0(gene," expression"),
	     xlab="Cancer type",
	     palette = c("blue","red") )
	p=p+rotate_x_text(60)
	p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")
	pdf(file=paste0(gene,".diff.pdf"),width=8,height=5)    #输出图片文件
	print(p1)
	dev.off()
}

