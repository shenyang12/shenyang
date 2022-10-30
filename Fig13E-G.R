
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

setwd("C:\\Users\\Administrator\\Desktop\\TP53_family_3\\02_methyDrug")
#引用包
library(limma)
rm(list=ls())
options(stringsAsFactors = FALSE)

##########################计算相关性
rm(list=ls())

methy=read.table("panMethyExp.txt", header=T, sep="\t", check.names=F)
methy=unique(methy)
row.names(methy)=methy[,1]
drug=read.table("panDrugExp.txt", header=T, sep="\t", check.names=F,row.names = 1)
sameGene=intersect(rownames(methy),rownames(drug))
sameMethy=methy[sameGene,]
samedrug=drug[sameGene,]
methySelect=sameMethy[,1:4]
drugSelect=samedrug[,4:ncol(samedrug)]
corTest=cbind(methySelect,drugSelect)
write.table(corTest,"CorTest.txt",sep="\t", quote=F,row.names = F)


rt <- read.delim("MethyDrug.txt",row.names = 1)

rt=na.omit(rt)#删除缺失值所在的行

drug=read.table("gene.txt",sep="\t",header=F)  
druglist=as.vector(drug[,1])
fix(rt) 

#x=as.numeric(rt[druglist,])
gene="TP53"
#x=as.numeric(t(rt)[gene,])

#rm(outTab)
#rm(corTab)
outTab=data.frame()
corTab=data.frame()

for (CancerType in levels (factor(rt[,"CancerType"]))){
  rt1=rt[(rt[,"CancerType"]==CancerType),]
  outVector=data.frame(CancerType)
  corVector=data.frame(CancerType)
  for(j in druglist){
    
    rt2=t(rt1)
    y=as.numeric(rt2[j,])
    x=as.numeric(rt2[gene,])
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
colnames(outTab)=c("CancerType",druglist)
write.table(outTab,file="geneCorTP53.pvalue.txt",sep="\t",row.names=F,quote=F)
#输出相关性结果
colnames(corTab)=c("CancerType",druglist)
write.table(corTab,file="geneCorTP53.cor.txt",sep="\t",row.names=F,quote=F)

##heatmap------------------------------------------------------------------------------
rm(list=ls())
library(reshape2)
library(RColorBrewer)


up <- read.table("geneCorTP53.pvalue.txt",sep = "\t",check.names = F,header = T,row.names=1)  #读取左上角的数据
dn <- read.table("geneCorTP53.cor.txt",sep = "\t",check.names = F,header = T,row.names=1)     #读取右下角数据
dn=t(dn)
up=t(up)
#设置颜色
colVector=c("#DDDDDD","#00FFFF","#FFFFFF","#FF0000")

#行名和列名
gene.level <- as.character(rownames(dn)) 
cancer.level <- as.character(colnames(dn))

#把行转为列
dn.long <- setNames(melt(dn), c('Gene', 'Cancer', 'Frequency'))
dn.long$Categrory <- "DN"
up.long <- setNames(melt(up), c('Gene', 'Cancer', 'Frequency'))
up.long$Categrory <- "UP"

#右下角颜色
dn.long$range <- cut(dn.long$Frequency, 
                     breaks = seq(floor(min(dn.long$Frequency)),
                                  ceiling(max(dn.long$Frequency)),0.01))
rangeMat1 <- levels(dn.long$range) # 提出分割区间
rbPal1 <- colorRampPalette(colors = c(colVector[2],"white",colVector[4]))
col.vec1 <- rbPal1(length(rangeMat1)); names(col.vec1) <- rangeMat1
dn.long$color <- col.vec1[as.character(dn.long$range)]

#左上角颜色
up.long$range <- cut(up.long$Frequency, breaks = seq(floor(min(up.long$Frequency)),ceiling(max(up.long$Frequency)),0.01)) 
rangeMat2 <- levels(up.long$range)
rbPal2 <- colorRampPalette(colors = c(colVector[1],colVector[3]))
col.vec2 <- rbPal2(length(rangeMat2)); names(col.vec2) <- rangeMat2
up.long$color <- col.vec2[as.character(up.long$range)]

dn.long=cbind(dn.long,dn.long$Frequency)
up.long=cbind(up.long,dn.long$Frequency)
#合并右下角和左上角
heatmat <- rbind.data.frame(dn.long,up.long) #汇总热图矩阵
pdf("heatmap1.TP53.pdf",width =7,height =12)
layout(mat=matrix(c(1,0,1,2,1,0,1,3,1,0),5,2,byrow=T),widths=c(length(cancer.level),2))

#热图绘制区域
par(bty="n", mgp = c(2,0.5,0), mar = c(5.1, 5.5, 3, 3),tcl=-.25,xpd = T)
x=as.numeric(factor(heatmat$Cancer,levels = cancer.level))
y=as.numeric(factor(heatmat$Gene,levels = gene.level))
#创建空白画布
plot(1,xlim=c(1,length(unique(x))+1),ylim=c(1,length(unique(y))+1),
     xaxs="i", yaxs="i",xaxt="n",yaxt="n",
     type="n",bty="n",xlab="",ylab="",
     main = "TP53",cex.main=2)
#填充颜色

for(i in 1:nrow(heatmat)) {
  if(heatmat$Categrory[i]=="DN") polygon(x[i]+c(0,1,1),y[i]+c(0,0,1),col=heatmat$color[i]) #填充左上角
  if(heatmat$Categrory[i]=="UP") {
    polygon(x[i]+c(0,1,0),y[i]+c(0,1,1),col=heatmat$color[i]) #填充右下角
    if(abs(heatmat$`dn.long$Frequency`[i])>=0.2 & heatmat$Frequency[i]<0.001){
      text(x[i]+0.5,y[i]+0.8,"***",cex=0.8)
    }else if(abs(heatmat$`dn.long$Frequency`[i])>=0.2 & heatmat$Frequency[i]<0.01){
      text(x[i]+0.5,y[i]+0.8,"**",cex=0.8)
    }else if(abs(heatmat$`dn.long$Frequency`[i])>=0.2 & heatmat$Frequency[i]<0.05){
      text(x[i]+0.5,y[i]+0.8,"*",cex=0.8)
    }
  }
}



#基因名和癌症名
axis(1,at = sort(unique(x)) + 0.5,labels = cancer.level,lty = 0,las = 2)  #添加x轴坐标
axis(2,at = sort(unique(y)) + 0.5,labels = gene.level,lty = 0,las = 1)    #添加y轴坐标
#mtext("Cancer types",side = 1,line = 3.5,cex=1.2)    #x轴名称

#绘制图例
par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
barplot(rep(1,length(col.vec2)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec2)),horiz=TRUE,
        axes = F, col=col.vec2)  # Loss
axis(4,at=c(1,ceiling(length(col.vec2)/2),length(col.vec2)),c(round(min(up),1),'Pvalue',round(max(up),1)),tick=FALSE)
par(mar=c(0,0,0,2),xpd = T,cex.axis=1.6)
barplot(rep(1,length(col.vec1)),border = NA, space = 0,ylab="",xlab="",ylim=c(1,length(col.vec1)),horiz=TRUE,
        axes = F, col=col.vec1)  # Gain
#axis(4,at=c(1,ceiling(length(col.vec1)/2),length(col.vec1)),c(round(min(dn),1),'Cor',round(max(dn),1)),tick=FALSE)
axis(4,at=c(1,ceiling(length(col.vec1)/2),length(col.vec1)),c(round(-1,1),'Cor',round(1,1)),tick=FALSE)

dev.off()


