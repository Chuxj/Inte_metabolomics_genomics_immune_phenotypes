library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd("~/Desktop/metabolites/Dec_2020/Manhattn/")

dat<-read.table("New_manha.txt",header=T,sep="\t")

d=transform(dat,logp=-log10(p))

d$platform<-factor(d$platform,levels = c("BM","GM","UM"))

df.tmp<- d %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(d, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


mypalette <- c("lightblue","royalblue")

col<-mypalette

#title<-"Manhattan Plot"

png("manhattan_combine_new.png",height = 540,width = 1080)


ggplot(df.tmp, aes(x=BPcum, y=logp)) +
  geom_point(data=subset(df.tmp, gwSig=="Yes"), color="red", size=0.5,alpha=0.5) +

  # Show all points
  geom_point(data=subset(df.tmp, gwSig!="Yes"),aes(color=as.factor(chr)), size=0.5) +
  scale_color_manual(values = rep(col, 22 )[1:22]) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  #  ggtitle(paste0(title)) +
  labs(x = "Chromosome") +
  labs(y = "-log10(p)")+
  
  # add genome-wide sig and sugg lines
#  geom_hline(yintercept = 7.30103, colour="grey") +
  #  geom_hline(yintercept = -7.30103, colour="red") +
  #  geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  
  # Add label using ggrepel to avoid overlapping
  #  geom_text_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=2, force=1.3) +
  
  
  # ADD facet
facet_grid(platform~.,scale="free")+
  
  
  # Custom the theme:
  theme_bw(base_size = 18) +
  theme( 
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18)
  )

dev.off()