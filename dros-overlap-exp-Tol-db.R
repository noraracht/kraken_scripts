require(ggplot2); require(scales); require(reshape2); require(Hmisc);
getwd()
setwd("/Users/admin/Documents/Skimming/tree_of_life/dros_contam_test/test_output")


dh = read.csv('my-data-output_true_H.csv',sep=",")
qplot(Overlap/100,error,data=dh,color=Filtered,shape=Filtered,linetype=Filtered)+geom_line()+
  facet_wrap(~percent(round(D,3)),ncol=1)+
  scale_y_continuous(name="Relative error in Skmer distance", labels=percent)+ #, breaks= c(-sqrt(c(0,0.01,0.05,0.2,0.5,1))) ,labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=percent,name=expression("Overlap (H)"))+
  scale_color_manual(values=c("black","#e31a1c"),labels=c("Before Kraken","After Kraken"),name="")+  
  scale_shape_manual(name="",values=c(6,2),labels=c("Before Kraken","After Kraken") )+
  scale_linetype_manual(name="",values=c(2,1),labels=c("Before Kraken","After Kraken"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
ggsave("E2-h-sqrt.pdf",width = 3,height = 7)


