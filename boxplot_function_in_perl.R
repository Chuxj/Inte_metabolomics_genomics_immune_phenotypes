library(getopt)

command=matrix(c('effect_allele',	'e',	1,	'character','alternative_allele',	'a',	1,	'character','working_directory',       'p',    1,      'character'),byrow=T,ncol=4)

opt=getopt(command)

EA<-opt$effect_allele
AA<-opt$alternative_allele
pwd<-opt$working_directory


dat<-read.table(paste(pwd,"/tmp",sep=""),row.names=1);

dat[1,]<-round(dat[1,])

dat<-t(dat)

dat<-as.data.frame(dat)

geno<-dat[,1]
geno<-gsub(geno,pattern="2",replacement=paste(EA,EA,sep=""))
geno<-gsub(geno,pattern="1",replacement=paste(EA,AA,sep=""))
geno<-gsub(geno,pattern="0",replacement=paste(AA,AA,sep=""))

dat[,1]<-factor(geno,levels=c(paste(EA,EA,sep=""),paste(EA,AA,sep=""),paste(AA,AA,sep="")))


library(ggplot2)
pdf(paste(pwd,"/",names(dat)[1],"_",names(dat)[2],"boxplot.pdf",sep=""))
ggplot(data=dat,aes(x=dat[,1],y=dat[,2],fill=dat[,1]))+geom_boxplot()+geom_point()+xlab(paste("genotype at",names(dat)[1]))+ylab(paste(names(dat)[2]))+guides(fill = guide_legend(title = NULL))+theme_bw()+scale_x_discrete(labels=c(paste(EA,EA,"(",table(dat[,1])[[1]],")",sep=""),paste(EA,AA,"(",table(dat[,1])[[2]],")",sep=""),paste(AA,AA,"(",table(dat[,1])[[3]],")",sep="")))
dev.off()
