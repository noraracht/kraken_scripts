require(ggplot2); require(scales); require(reshape2)

d=read.csv("summary_report_kraken_viral_dist_unmasked_0.0_corrected_k35.csv")
names(d)

qplot(value,C_sequences_prct/100,color=prct_species,data=m[m$variable=="mash_10K_distance" & m$value<.3,])+
  #facet_wrap(~variable,scales = "free_x",nrow=2,labeller = function(x) {x$"variable"=c("Marker-based phylogenetic distance","Genome-wide Skmer distance");x})+
  theme_bw()+theme(legend.position = "bottom")+
  geom_smooth(se =F,color="red", method="lm")+
  ylab(expression("Reads classified by Krarken"~(1-f[n])))+xlab("Distance to the closest reference species")+
  scale_y_continuous(labels=percent)+
  scale_color_continuous(name="Species Classification")
ggsave("classified-vs-match.pdf",width = 4*1.25, height = 3.5*1.25)

qplot(value,C_sequences_prct/100,color=prct_species,data=m)+
  facet_wrap(~variable,scales = "free_x",nrow=1,labeller = function(x) {x$"variable"=c("Marker-based phylogenetic distance","Genome-wide Skmer distance");x})+
  theme_bw()+theme(legend.position = "bottom")+
  geom_smooth(se =F,color="red", method="lm",data=m[m$value<0.3,])+
  ylab("Percent classified")+xlab("Distance to the closest reference species")+
  scale_y_continuous(labels=percent)+
  scale_color_continuous(name="Species Classification")
ggsave("classified-vs-match-full.pdf",width = 6, height = 4.5)

