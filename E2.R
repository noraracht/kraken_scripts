require(ggplot2); require(scales); require(reshape2)
d=read.csv('E2.csv')
head(d)
dm=cbind(melt(d[,c(1:3,4,6)],id.vars = 1:3),melt(d[,c(1:3,5,7)],id.vars = 1:3)[,5])
names(dm)[6] <- c("error")
head(dm)
qplot(CL/100,(error),data=dm,linetype=variable,shape=bin,color=variable)+geom_line()+
  facet_grid(Dist~.,scales="free_y")+
  scale_y_continuous(labels=percent,name="Relative error in Skmer distance")+scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]))+
  scale_color_brewer(name="Filtering", palette = "Set2")+  
  scale_linetype_manual(name="Filtering",values=c(2,1))+scale_shape_manual(name="M",values=c(5,3,4,15))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = c(.449,.91),legend.direction = "horizontal",panel.grid.major.y = element_line(linetype = 1,size=0.24,color="gray"))
ggsave("E2.pdf",width = 4.5,height = 6.5)

qplot(CL/100,sign(error)*sqrt(abs(error)),data=dm,linetype=variable,shape=bin,color=variable)+geom_line()+
  facet_grid(Dist~.,scales="free_y")+
  scale_y_continuous(name="Relative error in Skmer distance",breaks=c(-sqrt(0.05),sqrt(c(0,0.05,0.2,0.5,1,2,5))),labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]))+
  scale_color_brewer(name="", palette = "Set2")+  
  scale_linetype_manual(name="",values=c(2,1))+scale_shape_manual(name="M",values=c(5,3,4,15))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.24,color="gray"))
ggsave("E2-sqrt.pdf",width = 4.5,height = 6.5)

