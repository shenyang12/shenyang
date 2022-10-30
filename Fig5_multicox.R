setwd("C:\\Users\\Administrator\\Desktop\\reAnalysisTP53/results")
install.packages('forestplot')
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)
getStudies(cbio)
##查看所有癌症的缩写
getStudies(cbio)[["studyId"]]


####获取UCEC信息进行COX分析----------------------------------------


samps <- sampleLists(cbio, "ucec_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "ucec_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","CLINICAL_STAGE","GRADE"
                           , "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId"
)] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]
  
sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID
rownames(multicox)=multicox$ID
multicox=same[,c(2:7,10:12)]


write.table(multicox,"UCEC.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("UCEC_COX.txt",header=T,sep="\t",check.names=F,row.names=1)
#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="UCEC.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="UCEC",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()



####获取SKCM信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "skcm_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "skcm_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_PATHOLOGIC_TUMOR_STAGE","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_TUMOR_PATHOLOGIC_PT","AJCC_NODES_PATHOLOGIC_PN"
                               , "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:9,12:14)]


write.table(multicox,"SKCM.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("SKCM_COX2.txt",header=T,sep="\t",check.names=F,row.names = 1)



#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="SKCM.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="SKCM",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()




####获取PRAD信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "prad_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "prad_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","CLIN_M_STAGE","PATH_N_STAGE","PATH_T_STAGE", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:8,11:13)]


write.table(multicox,"PRAD.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("PRAD_COX.txt",header=T,sep="\t",check.names=F,row.names = 1)



#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="PRAD.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="PRAD",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()


####获取READ信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "coadread_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "coadread_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_TUMOR_PATHOLOGIC_PT","AJCC_PATHOLOGIC_TUMOR_STAGE", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:9,12:16)]
multicox=multicox[multicox$CancerType=="READ",]
multicox=multicox[,c(1:11)]
write.table(multicox,"READ.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("READ_COX2.txt",header=T,sep="\t",check.names=F,row.names = 1)



#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="READ2.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="READ2",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()



####获取PCPG信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "pcpg_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "pcpg_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:5,8:10)]

write.table(multicox,"PCPG.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("PCPG.txt",header=T,sep="\t",check.names=F,row.names = 1)

# 全部样品os_status=0
#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="PCPG.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="PCPG",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()





####获取LUSC信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "lusc_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "lusc_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_PATHOLOGIC_TUMOR_STAGE","AJCC_TUMOR_PATHOLOGIC_PT", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:9,12:14)]

write.table(multicox,"LUSC.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("LUSC_COX2.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="LUSC2.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="LUSC",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()

####获取LGG信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "lgg_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "lgg_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","GRADE", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:6,8:11)]

write.table(multicox,"LGG.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("LGG_COX.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="LGG.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="LGG",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()

####获取KIRC信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "kirc_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "kirc_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_PATHOLOGIC_TUMOR_STAGE","GRADE","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_TUMOR_PATHOLOGIC_PT", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:10,13:15)]

write.table(multicox,"KIRC.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("KIRC_COX3.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="KIRC.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="KIRC",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()


####获取HNSC信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "hnsc_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "hnsc_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_PATHOLOGIC_TUMOR_STAGE","GRADE","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_TUMOR_PATHOLOGIC_PT", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:10,13:15)]

write.table(multicox,"HNSC.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("HNSC_COX2.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="HNSC2.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="HNSC",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()




####获取CESC信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "cesc_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "cesc_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","CLINICAL_STAGE","GRADE","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_TUMOR_PATHOLOGIC_PT", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:10,13:15)]

write.table(multicox,"CESC.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("CESC_COX2.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="CESC.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="CESC",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()


####获取BRCA信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "brca_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "brca_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_PATHOLOGIC_TUMOR_STAGE","AJCC_METASTASIS_PATHOLOGIC_PM","AJCC_NODES_PATHOLOGIC_PN","AJCC_TUMOR_PATHOLOGIC_PT", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:9,12:14)]

write.table(multicox,"BRCA.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("BRCA_COX3.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS, OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="BRCA3.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="BRCA",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()




####获取acc信息进行COX分析----------------------------------------
rm(list=ls())
library(cBioPortalData)
library(DT)
cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs"
)



samps <- sampleLists(cbio, "acc_tcga")
samps[, c("category", "name", "sampleListId")]



#Obtaining Clinical Data

clinicalData=clinicalData(cbio, "acc_tcga")

DT::datatable(clinicalData)
write.table(clinicalData,"clinicaldata.txt",sep="\t",quote=F)


targetAnno <- clinicalData[, c("patientId","AJCC_PATHOLOGIC_TUMOR_STAGE","CLIN_M_STAGE","PATH_N_STAGE","PATH_T_STAGE", "OS_MONTHS", "OS_STATUS", "SEX","AGE","sampleId")] 

RPKM=read.table("panGeneExp.txt", header=T, sep="\t", check.names=F)
RPKM=as.matrix(RPKM)
rownames(RPKM)=RPKM[,1]
library(stringr)
ID=str_sub(rownames(RPKM),start=1,end=12)
rownames(RPKM)=ID

sameID=intersect(targetAnno$patientId, rownames(RPKM))

targetAnno=as.matrix(targetAnno)
rownames(targetAnno)=targetAnno[,1]

sametarget=targetAnno[sameID,]
sameexp=RPKM[sameID,]
same=cbind(sametarget,sameexp)

same=as.data.frame(same)
same$OS_STATUS=sapply(strsplit(same$OS_STATUS,":"),"[",1)
rownames(same)=same$ID

multicox=same[,c(2:9,12:14)]

write.table(multicox,"ACC.txt",sep="\t",quote=F)

#multicox分析-------------------------------------------

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.table("ACC_COX4.txt",header=T,sep="\t",check.names=F,row.names = 1)

#删掉正常样品
group=sapply(strsplit(rownames(rt),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
rt=rt[group==0,]

multiCox=coxph(Surv(OS_MONTHS,OS_STATUS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="ACC.pdf",
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           title="ACC",
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()






