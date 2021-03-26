setwd("~/Desktop/metabolites/")

dat<-read.table("AA_snp.txt",header=F)
names(dat)<-c("SNP","metabolites","beta","t-stat","p","FDR")
dat<-transform(dat,se=sqrt(beta^2/(qchisq(p,1,lower.tail = FALSE))))

library(devtools)        
#install_github("MRCIEU/TwoSampleMR")       
library(TwoSampleMR)

#devtools::install_github("MRCIEU/MRInstruments")

#library(MRInstruments)
#data(gwas_catalog)
#head(gwas_catalog)
#table(gwas_catalog$Phenotype)

#Crohn<-subset(gwas_catalog,Phenotype=="Crohn's disease")
#Crohn<-format_data(Crohn)

#Crohn<-transform(Crohn,id.exposure=rep("Crohn's_disease",nrow(Crohn)))
#Crohn<-Crohn[,-13]
#Crohn<-Crohn[c(1,)]

dat<-subset(dat,p<1e-5)

names(dat)<-c("SNP","id.exposure","beta.exposure","t-stat","pval.exposure","FDR","se.exposure")
list<-clump_data(dat,clump_r2 = 0.1)


effect_allele<-read.table("snp_effect_allele.txt",header=T)

exp<-merge(list,effect_allele,by.x="SNP")



out<-read.table("Crohns_disease_EUR_2015_26192919_hg19.txt",header=T)
names(out)<-c("SNP","MinorAllele","beta.outcome","pval.outcome")
out<-transform(out,se.outcome=sqrt(beta.outcome^2/(qchisq(pval.outcome,1,lower.tail = FALSE))))
 
out<-read.table("Neutrophiles_EUR_2016_27346689_hg19.txt",header=T)
names(out)<-c("SNP","MinorAllele","beta.outcome","pval.outcome")
out<-transform(out,se.outcome=sqrt(beta.outcome^2/(qchisq(pval.outcome,1,lower.tail = FALSE))))


out<-read.table("Eczema_2015_26482879_hg19.txt",header=T)
names(out)<-c("SNP","MinorAllele","beta.outcome","pval.outcome")
out<-transform(out,se.outcome=sqrt(beta.outcome^2/(qchisq(pval.outcome,1,lower.tail = FALSE))))


merged<-merge(exp,out,by.x="SNP")              
merged<-transform(merged,id.outcome=rep("Crohn's_disease",nrow(merged)))
merged<-transform(merged,mr_keep.exposure=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep.outcome=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep=rep(TRUE,nrow(merged)))

tmp<-merged[,-c(4,6)]

for (i in 1:nrow(tmp)){
  if(tmp$A1[i]!=tmp$MinorAllele[i]){tmp$beta.outcome[i]=0-tmp$beta.outcome[i]}else{tmp$beta.outcome[i]=tmp$beta.outcome[i]}
}

tmp<-transform(tmp,exposure=rep("g_303.2331",nrow(tmp)))
tmp<-transform(tmp,outcome=rep("Crohn's_disease",nrow(tmp)))

mr(tmp)

names(tmp)
names(tmp)<-c("SNP","id.outcome","beta.outcome","pval.outcome","se.outcome","A1","MinorAllele","beta.exposure","pval.exposure","se.exposure","id.exposure","mr_keep.exposure","mr_keep.outcome","mr_keep","outcome","exposure")
mr(tmp)




####FADS2
out<-read.table("~/Desktop/metabolites/FADS2.txt",header=T)
names(out)<-c("SNP","MinorAllele","beta.outcome","pval.outcome")



###### AA
exp<-read.table("N6_ArA_2014_24823311_rs.txt",header=T)
names(exp)<-c("SNP","A1","A2","EA","beta.exposure","pval.exposure")
exp<-exp[,-c(2,3)]
exp<-subset(exp,pval.exposure<1e-5)
exp<-transform(exp,se.exposure=sqrt(beta.exposure^2/(qchisq(pval.exposure,1,lower.tail = FALSE))))
exp<-clump_data(exp,clump_r2 = 0.1)


merged<-merge(exp,out,by.x="SNP")              
merged<-transform(merged,id.outcome=rep("Crohn's_disease",nrow(merged)))
merged<-transform(merged,mr_keep.exposure=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep.outcome=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep=rep(TRUE,nrow(merged)))

tmp<-merged

for (i in 1:nrow(tmp)){
  if(tmp$EA[i]!=tmp$MinorAllele[i]){tmp$beta.outcome[i]=0-tmp$beta.outcome[i]}else{tmp$beta.outcome[i]=tmp$beta.outcome[i]}
}

tmp<-transform(tmp,exposure=rep("g_303.2331",nrow(tmp)))
tmp<-transform(tmp,outcome=rep("Crohn's_disease",nrow(tmp)))

mr(tmp)


names(tmp)
names(tmp)<-c("SNP","EA","beta.outcome","pval.outcome","se.outcome","id.outcome","MinorAllele","beta.exposure","pval.exposure","se.exposure","id.exposure","mr_keep.exposure","mr_keep.outcome","mr_keep","outcome","exposure")
mr(tmp)


######w6
exp<-read.table("FAw6_2016_27005778_hg19.txt",header=T)
names(exp)<-c("SNP","EA","beta.exposure","pval.exposure")
exp<-subset(exp,pval.exposure<1e-5)
exp<-transform(exp,se.exposure=sqrt(beta.exposure^2/(qchisq(pval.exposure,1,lower.tail = FALSE))))
exp<-clump_data(exp,clump_r2 = 0.1)


merged<-merge(exp,out,by.x="SNP")              
merged<-transform(merged,id.outcome=rep("Crohn's_disease",nrow(merged)))
merged<-transform(merged,mr_keep.exposure=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep.outcome=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep=rep(TRUE,nrow(merged)))

tmp<-merged

for (i in 1:nrow(tmp)){
  if(tmp$EA[i]!=tmp$MinorAllele[i]){tmp$beta.outcome[i]=0-tmp$beta.outcome[i]}else{tmp$beta.outcome[i]=tmp$beta.outcome[i]}
}

tmp<-transform(tmp,exposure=rep("g_303.2331",nrow(tmp)))
tmp<-transform(tmp,outcome=rep("Crohn's_disease",nrow(tmp)))

mr(tmp)


names(tmp)
names(tmp)<-c("SNP","EA","beta.outcome","pval.outcome","se.outcome","id.outcome","MinorAllele","beta.exposure","pval.exposure","se.exposure","id.exposure","mr_keep.exposure","mr_keep.outcome","mr_keep","outcome","exposure")
mr(tmp)


####LA/EPA/DHA/DPA

exp<-read.table("N3_DPA_2011_21829377_rs.txt",header=T)
names(exp)<-c("SNP","A1","A2","EA","beta.exposure","pval.exposure")
exp<-exp[,-c(2,3)]
exp<-subset(exp,pval.exposure<1e-5)
exp<-transform(exp,se.exposure=sqrt(beta.exposure^2/(qchisq(pval.exposure,1,lower.tail = FALSE))))
exp<-clump_data(exp,clump_r2 = 0.1)


merged<-merge(exp,out,by.x="SNP")              
merged<-transform(merged,id.outcome=rep("Crohn's_disease",nrow(merged)))
merged<-transform(merged,mr_keep.exposure=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep.outcome=rep(TRUE,nrow(merged)))
merged<-transform(merged,mr_keep=rep(TRUE,nrow(merged)))

tmp<-merged

for (i in 1:nrow(tmp)){
  if(tmp$EA[i]!=tmp$MinorAllele[i]){tmp$beta.outcome[i]=0-tmp$beta.outcome[i]}else{tmp$beta.outcome[i]=tmp$beta.outcome[i]}
}

tmp<-transform(tmp,exposure=rep("g_303.2331",nrow(tmp)))
tmp<-transform(tmp,outcome=rep("Crohn's_disease",nrow(tmp)))

mr(tmp)


names(tmp)
names(tmp)<-c("SNP","EA","beta.outcome","pval.outcome","se.outcome","id.outcome","MinorAllele","beta.exposure","pval.exposure","se.exposure","id.exposure","mr_keep.exposure","mr_keep.outcome","mr_keep","outcome","exposure")
mr(tmp)



#ARA,EPA,DPA (sig >>CD)
#LA(sig <<CD)