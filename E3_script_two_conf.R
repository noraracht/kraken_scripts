require(ggplot2); require(scales); require(reshape2)
getwd()
setwd("/Users/admin/Documents/Skimming/tree_of_life/dros_contam_test")

d= read.csv('Drosophila_contam_both_species3.csv')
print (d)
dm = (melt(d[,c(1,2,grep(pattern = "*Dros*", names(d)))],id.vars = 1:2))
dm$Pair=""
dm[grep("sim_WXD1",dm$variable),"Pair"] = "Small"
dm[grep("sech_plant",dm$variable),"Pair"] = "Medium"
dm[grep("yakuba",dm$variable),"Pair"] = "Large"
dm$Kraken = grepl ("*ucseq*", dm$variable)
dm$k = 35
dm[grep (".1$", dm$variable),"k"] = 32
dm[grep (".2$", dm$variable),"k"] = 29
dm[grep (".3$", dm$variable),"k"] = 28
dm[!dm$Kraken,"k"] = "None"
names(dm)
print (dm)
dm = (merge(dm, dm[!dm$Kraken & dm$cl == 0 & dm$bin == "00-00",c(5,4)], by.x = "Pair", by.y = "Pair"))
print (dm)
names(dm)[8] = "D"
names(dm)[5] = "est"
dm$error = with(dm, (est-D)/D)
print (dm)
write.csv(dm,"/Users/admin/Documents/Skimming/tree_of_life/dros_contam_test/Drosophila_contam_both_species_formatted.csv", row.names = FALSE)

#This code should be used on formatted input.

install.packages("readxl")
library(readxl)
dm= read_excel('Drosophila_contam_both_species_formatted_withAlpha.xls')
print(dm)


qplot(cl/100,error,data=dm[dm$k != 29,],linetype=as.factor(confidence),shape=bin,color=Kraken)+geom_line()+
  facet_grid(D~k,scales="free_y")+
  scale_y_continuous(labels=percent,name="Relative error in Skmer distance")+scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]))+
  scale_color_brewer(name="Filtering", palette = "Set2")+  
  scale_shape_manual(name="M",values=c(5,3,4,15))+
  scale_linetype_manual(name="Filtering",values=c(1,2,3))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = c(.449,.91),legend.direction = "horizontal",panel.grid.major.y = element_line(linetype = 1,size=0.24,color="gray"))


qplot(cl/100,(error),data=dm[dm$k != 29,],linetype=as.factor(confidence),color=k,shape=k)+geom_line()+
  facet_grid(percent(round(D,3))~bin,scales="free",space = "free_x")+
  scale_y_continuous(name="Relative error in Skmer distance", labels = percent)+ #,breaks=c(-sqrt(0.05),sqrt(c(0,0.05,0.2,0.5,1,2,5))),labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=scales::percent_format(accuracy = 1),name=expression("Contamination level"~c[l]),breaks=c(0,0.05,.1,0.2,0.4,0.6))+
  scale_color_manual(values=c("#fecc5c","#fd8d3c","#e31a1c","black"))+  
  scale_linetype_manual(name="confidence",values=c(1,2,3),labels=c("0.00","0.05","None")) + 
  scale_shape_manual(values=c(0,2,6,19))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.24,color="gray"))
ggsave("E2_with_conf.pdf",width = 8,height = 7)

qplot(cl/100,sign(error)*sqrt(abs(error)),data=dm[dm$k != 29,],linetype=as.factor(confidence),color=k,shape=k)+geom_line()+
  facet_grid(percent(round(D,3))~bin,scales="free",space = "free_x")+
  scale_y_continuous(name="Relative error in Skmer distance",breaks=c(-sqrt(0.05),sqrt(c(0,0.05,0.2,0.5,1,2,5))),labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=scales::percent_format(accuracy = 1),name=expression("Contamination level"~c[l]),breaks=c(0,0.05,.1,0.2,0.4,0.6))+
  scale_color_manual(values=c("#fecc5c","#fd8d3c","#e31a1c","black"))+  
  scale_linetype_manual(name="confidence",values=c(1,2,3),labels=c("0.00","0.05","None"))+scale_shape_manual(values=c(0,2,6,19))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
ggsave("E2-sqrt_with_conf.pdf",width = 8,height = 7)

ds = (read.csv("drosophilaskims.csv"))
names(ds)

qplot(abs(bk_no_clean_up-fna_dist)/fna_dist-abs(ak_no_clean_up-fna_dist)/fna_dist,data=ds,binwidth=0.02)

qplot(abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist,data=ds,binwidth=0.02)
qplot(fna_dist,abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist,data=ds)
ggplot(aes(x=fna_dist,y=abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist),data=ds)+
  #geom_violin(aes(group=cut(ds$fna_dist,breaks=c(0,0.03,0.08,0.11,0.126,0.145,0.2))),scale="width")+
  geom_point(color="blue")+geom_smooth(se=F,method="lm",color="red")+
  theme_light()+
  scale_y_continuous(labels=percent,name="Decrease in relative error after Kraken")+
  scale_x_continuous(name=("D"),labels=percent)+
  facet_wrap(~s1)

ggplot(aes(x=abs(bk_cleaned-fna_dist)/fna_dist,y=abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist),data=ds)+
  #geom_violin(aes(group=cut(ds$fna_dist,breaks=c(0,0.03,0.08,0.11,0.126,0.145,0.2))),scale="width")+
  geom_point(color="blue")+geom_smooth(se=F,method="lm",color="red")+
  theme_light()+
  scale_y_continuous(labels=percent,name="Decrease in relative error after Kraken")+
  scale_x_continuous(name=("D"),labels=percent)#+
  facet_wrap(~s1)

qplot(abs(bk_no_clean_up-fna_dist)/fna_dist-abs(bk_cleaned-fna_dist)/fna_dist,data=ds)


