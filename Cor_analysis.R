setwd("/Users/xiaojing/Documents/working/database/500FG/molecular_phenotypes/")

gm<-read.table("/Users/xiaojing/Desktop/metabolites/gQTL/preprocess/raw_matrix.txt",header=T,row.names = 1,sep="\t")

gm<-log10(gm)

bm<-read.table("500FG_log2_brainshake_metabolites.txt",header=T,row.names = 1)
bm<-na.omit(bm)
gm<-bm


phe<-read.table("../phenotypes_xj.txt",header=T,row.names = 1,sep="\t")
phe<-phe[,c(1,2,11)]
flag<-which(is.na(phe$contraceptive))
phe[flag,3]<-0

phe<-na.omit(phe)

ol<-intersect(rownames(gm),rownames(phe))
gm<-gm[ol,]
phe<-phe[ol,]

tmp<-gm

for (i in names(gm)){
  lm<-lm(gm[,i]~phe$gender+phe$age+phe$contraceptive)
  tmp[,i]<-lm$residuals
}

gm<-tmp


cc<-read.table("~/Documents/working/projects/300DM/NMDS_500FG/ccount_500FG.txt",header=T,row.names = 1,sep="\t",check.names = F)

cc<-cc[,c(1:15,20:73)]

cc<-log2(cc)

cp<-read.table("500FG_log2_circulating_protein_levels.txt",header=T,row.names = 1)

ig<-read.table("500FG_log2_immunoglobulin_levels.txt",header=T,row.names = 1)

plt<-read.table("500FG_log2_platelet_count.txt",header=T,row.names = 1)

mAbu<-read.table("~/Desktop/metabolites/gQTL/figC/500FG_arcsine_transformed_taxonomy_abundance.txt",header=T,row.names = 2)
mAbu<-mAbu[,-1]

mPat<-read.table("~/Desktop/metabolites/gQTL/figC/500FG_microbiome_pathways.txt",header=T,row.names = 1)

cyto<-read.csv("../pheno_91cytokines_4Raul_log2.csv",header=T,row.names = 1)
cyto<-as.data.frame(t(cyto))

cc<-cyto

###########
ol<-intersect(rownames(gm),rownames(cc))
tmp_gm<-gm[ol,]
tmp_cc<-cc[ol,]

res<-matrix(NA,nrow=ncol(gm),ncol=ncol(cc))
res<-as.data.frame(res)
rownames(res)<-names(gm)
colnames(res)<-names(cc)

p<-res

for (i in names(gm)){
  for (j in names(cc)){
    test<-cor.test(tmp_gm[,i],tmp_cc[,j],method = "spearman")
    
    res[i,j]<-test$estimate
    p[i,j]<-test$p.value
  }
  
}

fdr<-as.matrix(p)
fdr<-as.numeric(fdr)
fdr<-p.adjust(fdr,method = "fdr")
fdr<-matrix(fdr,nrow=nrow(p))
fdr<-as.data.frame(fdr)
rownames(fdr)<-rownames(p)
names(fdr)<-names(p)

min<-apply(fdr,1,min)
flag<-which(min<0.05)

r<-res[flag,]
  
library(pheatmap)

col=colorRampPalette(c("blue","white", "tomato"))(10)

anno_cc<-read.table("~/Documents/working/projects/300DM/correlation_cc_cyto/anno_cellcounts.txt",check.names = F,row.names = 1,header = T,sep="\t")
names(anno_cc)<-"SubPopulation"

mycol_CellSubpopulation=brewer.pal(9, c("Set1"))
names(mycol_CellSubpopulation)<-levels(anno_cc$SubPopulation)
annotation_colors<-list(SubPopulation=mycol_CellSubpopulation)

load("~/Desktop/metabolites/cyto_colors.Rdata")
mycol_StimType=brewer.pal(5, c("Set2"))
mycol_Cytokine=c("#FBB4AE","#B3CDE3","#CCEBC5","#FDDAEC","#F2F2F2","aquamarine")
mycol_Time=c("lightblue","pink","lightgreen")
mycol_Stim_Assay=brewer.pal(3, c("Set3"))

names(mycol_Cytokine)<-c("IL1b","IL6","TNFA","IFNy","IL17","IL22")
names(mycol_StimType)<-c("Bacteria","Fungi","Non_microbial","Virus","TLR_ligands")
names(mycol_Time)<-as.character(unique(anno$Time))
names(mycol_Stim_Assay)<-as.character(unique(anno$Stim_Assay))

annotation_colors<-list(Cytokine=mycol_Cytokine,StimType=mycol_StimType,Time=mycol_Time,Stim_Assay=mycol_Stim_Assay)

r<-t(r)

pheatmap(r,color = col,annotation_row = anno_cc,annotation_colors = annotation_colors,show_colnames = F,border_color = F)

cl<-pheatmap(r,col=col,breaks = seq(-0.5,0.5,0.1),show_colnames = F)

res<-melt(res)
p<-melt(p)
fdr<-melt(fdr)

table<-data.frame(metabolite=rep(names(tmp_gm),ncol(cc)),phenotype=res$variable,r=res$value,p=p$value,fdr=fdr$value)

write.table(table,file="~/Desktop/metabolites/Dec_2020/BM_cyto.txt",quote=F,sep="\t",row.names = F)



