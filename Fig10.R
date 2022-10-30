setwd("C:\\Users\\Administrator\\Desktop\\reAnalysisTP53\\TMB_MSI_surv\\TMB/TP63/")
rm(list=ls())

library(survival)
library(survminer)

df<-read.csv("easy_input.csv",row.names = 1)
head(df)
df$futime=df$futime/365     #生存单位改成年
rt=df
#对基因进行循环
outTab=data.frame()
for(gene in colnames(rt)[3:(ncol(rt)-2)]){
  #对肿瘤类型进行循环
  for(CancerType in levels(factor(rt[,"CancerType"]))){
    rt1=rt[(rt[,"CancerType"]==CancerType),]
    #按照中位数划分表达量的高与低
    rt1$medA <- cut(rt1[,3],breaks=c(-Inf,median(rt1[,3]), Inf),
                    labels=c("0","1"))
    rt1$medB <- cut(rt1[,4],breaks=c(-Inf,median(rt1[,4]), Inf),
                    labels=c("0","1"))
    rt1$medGroup <- paste0(rt1$medA,rt1$medB)
    #按照两个基因表达量高低的排列组合分成4组
    rt1$medGroup <- ifelse(rt1$medGroup=='10','1',ifelse(rt1$medGroup=='00','2',ifelse(rt1$medGroup=="01","3","4")))
    
    rt1$group<-rt1$medGroup
    #定义结局变量，生存时间和结局
    y <- Surv(rt1$futime, rt1$fustat==1)
    
    
    comp1 <- survdiff(Surv(futime,fustat) ~ group, data = rt1)
    length=length(levels(factor(rt1[,"group"])))
    P=1-pchisq(comp1$chisq, df=length-1)
    if(P<0.001){
      P="P<0.001"
    }else{
      P=paste0("P=",sprintf("%.03f",P))
    }
    
    
    
    #logrank检验两两比较的结果
    comp <- pairwise_survdiff(Surv(futime,fustat) ~ group, data = rt1)
    pvalue<-as.vector(unlist(comp$p.value))
    
    #将p value和两两比较的组别生成两列的数据集
    name <- array(dim=c(3,3))
    for(i in 1:3){
      for(j in 2:4){
        name[i,j-1] = print(paste(i,j,sep = " vs "))
      }
    }
    
    pvalue_name <- as.vector(t(name))
    logrank <- data.frame(pvalue_name,pvalue)
    logrank
    #挑选p value小于0.05的记录
    logrank_sig <- subset(logrank, pvalue<0.05)
    #如果p值太小，就写“<0.0001”，否则保留小数点后4位
    logrank_sig$pvalue <- lapply(logrank_sig$pvalue,function(i)
      ifelse (i<0.0001,"<0.0001",round(i, 4)))
    logrank_sig
    
    #进行COX回归，导出每组对应的HR值，并将1赋值给对照组
    coxph.fit <- coxph(y ~ as.factor(group), data = rt1) 
    hr <- round(coef(summary(coxph.fit))[,2],3)
    HR <- c(1,as.vector(unlist(hr)))
    
    
    #KM曲线的绘制
    kmfit <- survfit(y~rt1$group,)
    
    #写legend
    A = c("Low; ", "High;", "Low; ", "High;")
    B = c("Low; ", "Low; ", "High;", "High;")
    Group = c(1,2,3,4)
    text.legend1 <- paste0(colnames(rt1)[3], " = ", A, colnames(rt1)[4], " = ", B, " Group = ", Group, ", HR = ", HR)
    text.legend2 <- paste0(logrank_sig$pvalue_name," : ",logrank_sig$pvalue)
    text.legend3 <- P
    #自定义足够你用的颜色
    mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
    
    #画图，并保存到pdf文件
    pdf(paste0(gene,"_",CancerType,".pdf"),width = 7,height = 6)
    par(xpd = T, 
        mar = par()$mar + c(0,0,0,5)#,#在右侧给图例留出空间
        #cex.axis=0.8, cex.lab=1, cex.main=1, cex.sub=1 #修改字体大小
    )
    plot(kmfit, 
         col= mycol, 
         lwd = 1.4,#线的粗细
         #pch = "o", #如果你也想让观测点的形状是o，就运行这行
         xlab="years from diagnosis", 
         ylab="Survival Distribution Function", 
         main=paste0("Cancer: ",CancerType),
         mark.time=TRUE)
    legend("bottomleft", lty=c("solid","solid","solid","solid"),
           col=mycol, 
           legend=text.legend1, 
           bty="n", 
           lwd = 2,
           cex=0.8)
    legend("topright", 
           inset=c(-0.3,0), #图例画到图外面
           legend=c("Pairwise comparison",text.legend2), 
           bty="n", 
           cex=0.8)
    legend("topright", 
           inset=c(0.1,0), 
           legend=c(text.legend3), 
           bty="n", 
           cex=1.2)
    
    dev.off()
    
  }
}

