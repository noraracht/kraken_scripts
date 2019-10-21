require(ggplot2); require(scales); require(reshape2)

d=read.csv("summary_report_kraken_viral_dist_unmasked_0.0_corrected_k35.csv")
names(d)

qplot(value,C_sequences_prct/100,color=prct_species/100,data=m[m$variable=="mash_10K_distance" & m$value<.3,])+
  #facet_wrap(~variable,scales = "free_x",nrow=2,labeller = function(x) {x$"variable"=c("Marker-based phylogenetic distance","Genome-wide Skmer distance");x})+
  theme_classic()+theme(legend.position = c(.6,.9),legend.direction="horizontal",panel.border  = element_rect(fill=NA,size = 1))+
  geom_smooth(se =F,color="red", method="lm")+
  ylab(expression("Recall"~(1-f[n])))+
  scale_y_continuous(labels=percent)+
  scale_x_continuous(labels=percent,name="Distance to the closest match (M)")+
  scale_color_continuous(name="Species\n classification",labels=percent)
ggsave("classified-vs-match.pdf",width = 4*1, height = 4*1)

qplot(value,C_sequences_prct/100,color=prct_species,data=m)+
  facet_wrap(~variable,scales = "free_x",nrow=1,labeller = function(x) {x$"variable"=c("Marker-based phylogenetic distance","Genome-wide Skmer distance");x})+
  theme_bw()+theme(legend.position = "bottom")+
  geom_smooth(se =F,color="red", method="lm",data=m[m$value<0.3,])+
  ylab("Percent classified")+xlab("Distance to the closest reference species")+
  scale_y_continuous(labels=percent)+
  scale_color_continuous(name="Species Classification",labels=percent)
ggsave("classified-vs-match-full.pdf",width = 6, height = 4.5)

qplot(C_sequences_prct/100,prct_species/100,data=m[m$variable=="mash_10K_distance" & m$value<.3,])+
  geom_abline(linetype=2,color="red")+theme_classic()+
  scale_y_continuous(labels=percent, name="Classified at any level")+
  scale_x_continuous(labels=percent, name ="Classified at the species level")
ggsave("domain-species.pdf",width = 4, height = 4)

