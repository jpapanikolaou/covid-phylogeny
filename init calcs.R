library(ape)
library(data.table)
library(ggplot2)
library(magrittr)

library(rentrez)
library(DECIPHER)
library(ggpubr)
library(tictoc)
source('scripts/utils.R')
theme_set(theme_bw(base_size=15))

data("RESTRICTION_ENZYMES")
# example <- data.table(rep(c('x','y'),each=4),sample(1:100,size=8,replace=F),c('l','l','l','o'))
# colnames(example) <- c('a','b','c')
# example[,-sum(b/100*log(b/100)),by='a']
# help('data.table')

all_re_Fragments <- fread('data/coronavirus_all_re_fragments.csv')[,2:6]
Fragments <- fread('data/coronavirus_bsmbi_bsai_fragment_lengths.csv')[,2:6]

entropy_frags = Fragments[,list(entropy=-sum(fragment_lengths/unique(genome_length)*log(fragment_lengths/unique(genome_length))),
                                no_fragments=.N,genome_length=unique(genome_length)),by=c('tip.label')]
entropy_all_re_frags = all_re_Fragments[,list(entropy=-sum(fragment_lengths/unique(genome_length)*log(fragment_lengths/unique(genome_length))),
                                no_fragments=.N,genome_length=unique(genome_length)),by=c('tip.label','Restriction_Enzyme')]
Engineered_CoVs <- read.csv('data/CoV Infectious Clones entropy.csv')[c(2,4),] %>% as.data.table
#Engineered_CoVs[,species:=rVirus]
cls <- gg_color_hue(nrow(Engineered_CoVs))


g_recomb=ggplot(entropy_all_re_frags[no_fragments<=30 & !grepl('SARS2',tip.label)],aes(no_fragments,entropy))+
  geom_boxplot(aes(x=factor(no_fragments)),lwd=2,col='darkgrey')+
  # geom_jitter(alpha=0.07)+  ### makes plot too busy
  scale_x_discrete('Number of Fragments')+
  scale_y_continuous(name='Entropy')+
  geom_point(data=entropy_frags[tip.label=='SARS2-WHu1'],color='red',cex=7,pch=18)+
#   annotate(geom='segment',x=21.8,xend=6,y=.57,yend=.254,color='red',lwd=1.5)+
#   annotate(geom='text',x=24.5,y=.6,label='SARS-CoV-2',color='red',size=12)+
  geom_point(data=Engineered_CoVs,aes(color=`virus`),cex=5,pch=18)+
  scale_color_manual(values=cls)+
  ggtitle('SARS-CoV-2 + Known Reverse-Engineering Viruses')+
#   theme(legend.position=c(0.6,0.8))+
#   geom_line(data=ideal,col='darkred',lwd=1.5)+
#   annotate(geom='segment',x=4.4,xend=4.4,y=1/4.4,yend=9/30,lwd=1.5,col='darkred')+
#   annotate(geom='segment',x=8.6,xend=8.6,y=1/8.6,yend=9/30,lwd=1.5,col='darkred')+
#   annotate(geom='segment',x=4.4,xend=8.6,y=9/30,yend=9/30,lwd=1.5,col='darkred')+
  geom_point(data=entropy_frags[tip.label=='SARS2-WHu1'],color='red',cex=7,pch=18)

  
g_recomb
