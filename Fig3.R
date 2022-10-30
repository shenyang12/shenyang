setwd("C:\\Users\\admin\\Desktop\\reAnalysisTP73\\cox")



library(survival)
                    
rt=read.table("TPMexpTime.txt",header=T,sep="\t",check.names=F,row.names=1)    #????????????

exprSet=rt[,(3:5)]

exprSet=log2(exprSet+1)

a=cbind(rt[,1:2],exprSet)
rt=cbind(a,rt[,6:7])
rt$futime=rt$futime/365

allGene=colnames(rt)[3:(ncol(rt)-2)]
allType=levels(factor(rt[,"CancerType"]))
pdf(file="forest.pdf", width = 15,height =10)
geneNum=length(allGene)
n=length(allType)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(1:(geneNum+1),nc=geneNum+1))
forestCol=rainbow(geneNum)


xlim = c(0,1)
par(mar=c(4.5,1,2,2))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text(0.5,n:1,allType,adj=0,cex=1.8)


i=0
for(gene in colnames(rt)[3:(ncol(rt)-2)]){
  i=i+1
  
  outTab=data.frame()
  for(CancerType in levels(factor(rt[,"CancerType"]))){
    rt1=rt[(rt[,"CancerType"]==CancerType),]
    cox=coxph(Surv(futime, fustat) ~ rt1[,gene], data = rt1)
    coxSummary = summary(cox)
    coxCoef=coxSummary$coefficients
    coxP=coxCoef[,"Pr(>|z|)"]
    outTab=rbind(outTab,
                 cbind(cancer=CancerType,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP) )
  }
  outFile=paste0(gene,".cox.txt")
  write.table(outTab,file=outFile,sep="\t",row.names=F,quote=F)  
  
  
  coxRT <- read.table(outFile,header=T,sep="\t",row.names=1,check.names=F)
  hr <- sprintf("%.3f",coxRT$"HR")
  hr[as.numeric(hr)<0.001]=0.001
  hrLow  <- sprintf("%.3f",coxRT$"HR.95L")
  hrLow[as.numeric(hrLow)<0.0001]=0.0001
  hrHigh <- sprintf("%.3f",coxRT$"HR.95H")
  hrHigh[as.numeric(hrHigh)>1000]=1000
  
  
  par(mar=c(4.5,1,2,2),mgp=c(2,0.5,0))
  LOGindex =5
  hrLow = log(as.numeric(hrLow),LOGindex)
  hrHigh = log(as.numeric(hrHigh),LOGindex)
  hr = log(as.numeric(hr),LOGindex)
  xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.01,col=forestCol[i],lwd=3)
  segments(log(1,LOGindex), 0, x1 = log(1,LOGindex), y1 = n)
  points(as.numeric(hr), n:1, pch = 15, col = forestCol[i], cex=2)
  a1 = axis(1,labels=F,tick=F)
  axis(1,a1,5^a1)
  text(log(1,LOGindex),n+1,gene,cex=2,xpd=T)
  if(i==round(geneNum/2,0)){
    text(log(1,LOGindex),-1.6,"Hazard Ratio",cex=2.5,xpd=T)
  }
}
dev.off()

########针对每个基因绘制forest plot

rm(list=ls())


library(forestplot)
library(survival)

  coxFile="TP73.cox.txt"
  #读取输入文件
  rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  data=as.matrix(rt)
  coxSummary=data[,1:3]
  hr=sprintf("%.3f",coxSummary[,"HR"])
  hrLow=sprintf("%.3f",coxSummary[,"HR.95L"])
  hrHigh=sprintf("%.3f",coxSummary[,"HR.95H"])
  pVal=data[,"pvalue"]
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))

  
  forest <- data.frame(Cancer = rownames(coxSummary),
                       HR = paste0(round(coxSummary[,"HR"],3),"(",
                                   round(coxSummary[,"HR.95L"],3),"-",round(coxSummary[,"HR.95H"],3),")"),
                       P_value = round(data[,"pvalue"],5))
  labeltext <- as.matrix(forest) 
  forest$HR_95L <- coxSummary[,"HR.95L"]
  forest$HR_95H <- coxSummary[,"HR.95H"]
  colnames(labeltext)

  singleCox <- forestplot(labeltext,  mean = coxSummary[,"HR"],  # 图形文本与HR部分 
                          graph.pos=3, #为箱线图所在的位置
                          lower =  forest$HR_95L, upper = forest$HR_95H, # 95%CI
                          zero = 1.0, #箱线图中基准线的位置
                          title="OS for TP73",
                          lwd.zero = 2, # 设置无效线的横坐标和线条宽度
                          xlab = "<Hazard ratio>", xlog = TRUE,
                          xticks = c(0.03125,0.0625,0.125,0.25,0.5,1,2,4,8,16,32),  # 设置x轴刻度ticks
                          xticks.digits = 8,clip = c(0.05,100),cex=0.9, lineheight = "auto",
                          lwd.xaxis = 2, boxsize = 0.5, lty.ci = 1, lwd.ci = 2, # 森林图置信区间的线条类型和宽度
                          ci.vertices = TRUE,  # 森林图置信区间两端添加小竖线，默认FALSE
                          ci.vertices.height = 0.2, # 设置森林图置信区间两端小竖线的高度，默认0.1
                          align = "l",  # 文字对齐方式，"l"、"r"和"c"分别为左对齐，右对齐和居中对齐
                          col=fpColors(box="blue",lines = 'black',zero = '#7AC5CD'),# 设置坐标轴标题大小
                          txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1.1)))
  singleCox
  
  
  pdf("TP73.OS.pdf",width = 10)
  singleCox
  dev.off()
  
  
