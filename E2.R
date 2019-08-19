require(ggplot2); require(scales); require(reshape2)
d=read.csv('E2.csv')
head(d)
dm=cbind(melt(d[,c(1:3,4,6)],id.vars = 1:3),melt(d[,c(1:3,5,7)],id.vars = 1:3)[,5])
names(dm)[6] <- c("error")
head(dm)
qplot(CL/100,(error),data=dm,linetype=variable,shape=bin,color=variable)+geom_line()+
  facet_grid(Dist~.,scales="free_y")+theme_classic()+
  scale_y_continuous(labels=percent,name="Relative error in Skmer distance")+scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]))+
  scale_color_brewer(name="Filtering", palette = "Set2")+  
  scale_linetype_manual(name="Filtering",values=c(2,1))+scale_shape(name="M")+
  theme(legend.position = c(.45,.93),legend.direction = "horizontal")
ggsave("E2.pdf",width = 4.5,height = 6.5)

