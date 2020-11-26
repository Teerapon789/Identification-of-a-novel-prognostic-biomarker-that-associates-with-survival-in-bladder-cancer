## "TCGA Analysis For Biomarker in BLCA
## Auther: Teerapon Sahwangarrom
## date: 26/11/2019

### Importing RPPA data into R
rppa.blca <- read.table("/Users/teeraponsahwangarrom/RPPA", sep="\t", head=T,row.names = 1)
dim(rppa.blca)
rppa.blca <- as.matrix(rppa.blca)
class(rppa.blca)


### Importing phenotype information (clinical data)
clin.data <- read.table("/Users/teeraponsahwangarrom/BLCA_clinicalMatrix", sep="\t", head=T, row.names = 1)
class(clin.data)
rownames(clin.data) <- gsub(rownames(clin.data),pattern="-",replace=".")
dim(clin.data)

sample.normal <- clin.data$sample_type=="Solid Tissue Normal"
sample.ID.normal <- rownames(clin.data[sample.normal,])
sample.ID.normal
colnames(rppa.blca) %in% sample.ID.normal

### Performing survival analysis
library(survival)
library(survminer)
os.time <- clin.data[colnames(rppa.blca),"OS.time"]
os.event <- as.numeric(clin.data[colnames(rppa.blca), "vital_status"]=="DECEASED")
blca.os <- Surv(os.time,os.event)
coxph(blca.os ~ rppa.blca[1])

all.hrs <- rep(NA,nrow(rppa.blca))
all.pvals <- rep(NA,nrow(rppa.blca))
all.lci <- rep(NA,nrow(rppa.blca))
all.uci <- rep(NA,nrow(rppa.blca))
for(i in 1:nrow(rppa.blca)){
  coxphmodel <- coxph(blca.os ~ rppa.blca[i,])
  all.hrs[i] <- summary(coxphmodel)$coef[1,2]
  all.pvals[i] <- summary(coxphmodel)$coef[1,5]
  all.lci[i] <- summary(coxphmodel)$conf.int[1,3]
  all.uci[i] <- summary(coxphmodel)$conf.int[1,4]
  
}

rppa.coxph.df <- data.frame(Protein=rownames(rppa.blca), HR=all.hrs, LCI = all.lci, UCI = all.uci, p.value=all.pvals)
rppa.coxph.df$adj.p <- p.adjust(all.pvals)
rppa.coxph.df <- rppa.coxph.df[order(rppa.coxph.df$p.value, decreasing=F),]
rppa.coxph.df[1:10,]

Annexin.high <- as.numeric(rppa.blca["Annexin-1-M-E",]>median(rppa.blca["Annexin-1-M-E",]))
surv_pvalue(survfit(blca.os ~ Annexin.high),data=as.data.frame(rppa.blca), method="survdiff")
plot(survfit(blca.os ~ Annexin.high), col=c("black","red"),lwd=2,mark.time=T,main="Annexin-1-M-E")
legend("topright",legend=c("Annexin-high", "Annexin-low"), col=c("red","black"),lwd=2)
boxplot(rppa.blca["Annexin-1-M-E",]~Annexin.high, main="Annexin.High vs Annexin.Low for RPPA")

### Correlation between RNA and RPPA

rna.blca <- read.table("/Users/teeraponsahwangarrom/HiSeqV2",sep="\t",head=T,row.names = 1)
dim(rna.blca)
class(rna.blca)

 ### Subtracting patient samples between RPPA and RNAseq
wh <- which(is.element(colnames(rppa.blca),colnames(rna.blca)))
rppa.blca1 <- rppa.blca[,wh]
rna.blca1 <- rna.blca[,colnames(rppa.blca1)]
rppa.blca1 <- as.matrix(rppa.blca1)
rna.blca1 <- as.matrix(rna.blca1)
dim(rna.blca1)
cor(rppa.blca1["Annexin-1-M-E",],rna.blca1["ANXA1",])
cor.test(rppa.blca1["Annexin-1-M-E",],rna.blca1["ANXA1",])
plot(rppa.blca1["Annexin-1-M-E",]~ rna.blca1["ANXA1",],main="The correlation plot between RNAseq and RPPA")
boxplot(rppa.blca1["Annexin-1-M-E",],rna.blca1["ANXA1",])