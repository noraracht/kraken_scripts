require(ggplot2); require(scales); require(reshape2)

d=read.csv("summary_report_kraken_viral_dist_unmasked_0.0_corrected_k35.csv")
names(d)
#melt((d[,c(2,8,30,3,38)]),id.vars = c("Exclude","C_sequences_prct","prct_species"))
qplot(value,C_sequences_prct/100,color=prct_species,data=m[m$variable!="mash_10K_distance" | m$value<0.5,])+facet_wrap(~variable,scales = "free_x",labeller = function(x) {x$"variable"=c("Phylogenetic","Sequence (Mash computed)");x})+theme_bw()+geom_smooth(se =F,color="red")+ylab("Percent classified")+scale_y_continuous(labels=percent)+xlab("Distance to the closest reference species")+scale_color_continuous(name="Classified at species level")+theme(legend.position = "bottom")
ggsave("graph3.pdf",width = 5.5, height = 3.5)
