library(ggplot2)
library(reshape2)
require(gtable)
library(gridExtra)

setwd('C:/Users/Jennifer/GoogleDrive/Nina_new/clone_size_paper/code_for_figures/')

# subplot #1 - Fraction of mixed clones as function of cutoff - merged trees of 2 naive IgM donors
df <- read.csv('merged_perc_by_cut.csv')
p1 <- ggplot(df, aes(x=cutoff, y=X))+
  geom_line()+
  coord_cartesian(xlim = c(1,20), ylim=c(0,0.25))+
  ylab('Fraction of mixed clones')+
  xlab('Cutoff')+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), # remove grid
        panel.grid.major=element_blank())

# subplot #2 - clone size distribution as function of cutoff - single IgG donor
df <- read.csv('clone_size_dist_by_cut_IGG_ALL_UNI.csv')

p2 <- ggplot(df, aes(x=CloneSize, y=Frequency, group=Cutoff, colour=Cutoff))+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~Seq, ncol=1)+
  coord_cartesian(xlim = c(1,1000))+
  ylab('Clone probability')+
  xlab('Clone size')+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), # remove grid
        panel.grid.major=element_blank(),
        strip.background = element_rect(fill="white"))

# subplot #3 - Hamming distance between ancestral sequences in clones (naive IgM - pfizer)
df <- read.csv('hamming_fig1_data.csv')
p3 <- ggplot(df, aes(x=avg.hamm))+
  geom_histogram()+
  coord_cartesian(ylim = c(0,900))+
  theme_bw()+
  xlab('Average hamming distance')+
  theme(panel.grid.minor=element_blank(), # remove grid
        panel.grid.major=element_blank())
 
# subplot #4 - Hamming distance between ancestral sequences in clones (naive IgM - pfizer)
df <- read.csv('frac_sep_synthetic.csv')
df <- df[df$sample==unique(df$sample)[4],]

p4 <- ggplot(df, aes(x=cutoff, y=value))+
  geom_line()+
  theme_bw()+
  ylab('Fraction of separated sequences')+
  xlab('Cutoff')+
  coord_cartesian(xlim = c(0,20), ylim=c(0,1))+
  theme(panel.grid.minor=element_blank(), # remove grid
        panel.grid.major=element_blank())

grid.arrange(p1, p2, p3, p4, ncol = 2)
