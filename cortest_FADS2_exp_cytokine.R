cyto<-read.csv("~/Desktop/PhD data/500fg/pheno_91cytokines_4Raul.csv",header=T,row.names=1)
exp<-read.table("~/Desktop/PhD data/500fg/500FG_RNAseq_TMM_normalized_read_counts.txt",header=T,row.names=1)

FADS2<-exp$ENSG00000134824
FADS2<-as.data.frame(FADS2)
rownames(FADS2)<-rownames(exp)
names(FADS2)<-"FADS2"
#remove dup
FADS2<-FADS2[-c(10,20,30,40,50,100),]
FADS2<-as.data.frame(FADS2)
rownames(FADS2)<-rownames(exp)[-c(10,20,30,40,50,100)]
rm(exp)

cyto<-t(cyto[,-1])

ol<-intersect(rownames(FADS2),rownames(cyto))

FADS2<-FADS2[ol,]
cyto<-cyto[ol,]
cyto<-as.data.frame(cyto)

#cor.test

p=rep(1,ncol(cyto))
r=rep(0,ncol(cyto))

for (i in 1:ncol(cyto)){
  t<-cor.test(FADS2,cyto[,i],method = "spearman")
  p[i]<-t$p.value  
  r[i]<-t$estimate[[1]]
}

res<-data.frame(p=p,r=r)
rownames(res)<-names(cyto)


#draw
plot(FADS2)
abline(h=350)
#FADS2 has 2 group


#IL22_MTB_PBMC_7days
cyto_2<-cyto[which(FADS2<350),"IL22_MTB_PBMC_7days"]
cyto_1<-cyto[which(FADS2>350),"IL22_MTB_PBMC_7days"]

#get the p value
t.test(cyto_1,cyto_2)

#fill the p value
boxplot(cyto_2,cyto_1,names = c("FADS2 low","FADS2 high"),main="IL22_MTB_PBMC_7days")
text(1.5,13,"p=0.042(t.test)")
