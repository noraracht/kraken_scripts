require(ggplot2); require(scales); require(reshape2)

d= read.csv('Drosophila_contam_both_species.csv')

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
dm = (merge(dm, dm[!dm$Kraken & dm$cl == 0 & dm$bin == "00-00",c(5,4)], by.x = "Pair", by.y = "Pair"))
names(dm)[8] = "D"
names(dm)[5] = "est"
dm$error = with(dm, (est-D)/D)



qplot(cl/100,error,data=dm[dm$k != 29,],linetype=Kraken,shape=bin,color=Kraken)+geom_line()+
  facet_grid(D~k,scales="free_y")+
  scale_y_continuous(labels=percent,name="Relative error in Skmer distance")+scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]))+
  scale_color_brewer(name="Filtering", palette = "Set2")+  
  scale_linetype_manual(name="Filtering",values=c(2,1))+scale_shape_manual(name="M",values=c(5,3,4,15))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = c(.449,.91),legend.direction = "horizontal",panel.grid.major.y = element_line(linetype = 1,size=0.24,color="gray"))


qplot(cl/100,(error),data=dm[dm$k != 29,],linetype=Kraken,color=k,shape=k)+geom_line()+
  facet_grid(percent(round(D,3))~bin,scales="free",space = "free_x")+
  scale_y_continuous(name="Relative error in Skmer distance", labels = percent)+ #,breaks=c(-sqrt(0.05),sqrt(c(0,0.05,0.2,0.5,1,2,5))),labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]),breaks=c(0,0.05,.1,0.2,0.4,0.6))+
  scale_color_manual(values=c("#fecc5c","#fd8d3c","#e31a1c","black"))+  
  scale_linetype_manual(name="",values=c(2,1),labels=c("Before Kraken","After Kraken"))+scale_shape_manual(values=c(0,2,6,19))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.24,color="gray"))
ggsave("E2.pdf",width = 8,height = 7)

qplot(cl/100,sign(error)*sqrt(abs(error)),data=dm[!dm$k %in% c(29,32),],linetype=Kraken,color=k,shape=k)+geom_line()+
  facet_grid(percent(round(D,3))~bin,scales="free",space = "free_x")+
  scale_y_continuous(name="Relative error in Skmer distance",breaks=function(x) {if (x[2]<1.1) c(-sqrt(0.01),sqrt(c(0,0.01,0.05,0.2,0.5))) else c(-sqrt(0.05),sqrt(c(0,0.05,0.2,0.5,2,5)))},labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=percent,name=expression("Contamination level"~c[l]),breaks=c(0,0.05,.1,0.2,0.4,0.6))+
  scale_color_manual(values=c("#fecc5c",#"#fd8d3c",
                              "#e31a1c","black"))+  
  scale_linetype_manual(name="",values=c(2,1),labels=c("Before Kraken","After Kraken"))+scale_shape_manual(values=c(0,2,6,19))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
ggsave("E2-sqrt.pdf",width = 8,height = 7)



qplot(cl/100,(abs(error)),data=dm[!dm$k %in% c(29,28,32) & dm$Pair == "Medium"  & dm$cl <30,],linetype=Kraken,color=bin,shape=bin)+geom_line()+
  #facet_grid(percent(round(D,3))~bin,scales="free",space = "free_x",labeller = label_both)+
  scale_y_continuous(name="Relative error in Skmer distance",label=percent)+
  scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")"))+
  scale_color_brewer(palette = "RdBu",direction = -1)+ #Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
  scale_linetype_manual(name=element_blank(),values=c(2,1),labels=c("Before Filtering","After Filtering"))+scale_shape_manual(values=c(0,2,6,1))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = c(.25,.62),panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
ggsave("E2-Medium.pdf",width = 3.7,height = 3.7)









####################################

bk = (read.csv('ref_dist_mat_bk_kraken_std_conf0.0_k35.txt.xls',sep="\t")); bk$Filtered=FALSE
ak = (read.csv('ref_dist_mat_ak_kraken_std_conf0.0_k35.txt.xls',sep="\t")); ak$Filtered=TRUE
names(ak)=names(bk)
dh=rbind(bk,ak)

dh=(melt(dh))

db = dh[dh$variable %in%  levels(dh$variable)[1:5] & dh$sample %in%  levels(dh$variable)[1:5] ,]
dh = dh[!dh$variable %in%  levels(dh$variable)[1:5],]

dh$sample = sub(pattern = "ucseq_", "", dh$sample)
dh = dh[grep("s1_",dh$sample,),]
dh = dh[grep("s2_",dh$variable,),]

dh$Overlap =  as.numeric(sub(".*overlap_","",dh$variable))
dh$Overlap2 =  as.numeric(sub(".*overlap_","",dh$sample))


dh$s2 = (sub("_overlap.*","",sub("s2_","",dh$variable)))
dh$s1 = (sub("_overlap.*","",sub("s1_","",dh$sample)))

dh = dh[dh$Overlap == dh$Overlap2, ]

dh$D = 0
dh[dh$s2=="Dros_sech","D"] = 0.021215276
dh[dh$s2=="Drosophila_sim_WXD1","D"] = 0.001933919
dh[dh$s2=="Drosophila_yakuba","D"] = 0.062705933


dh$error = with(dh, (value-D)/D)

qplot(Overlap/100,error,data=dh,color=Filtered,shape=Filtered)+geom_line()+
  facet_wrap(~percent(round(D,3)))+
  scale_y_continuous(name="Relative error in Skmer distance")+ #, breaks= c(-sqrt(c(0,0.01,0.05,0.2,0.5,1))) ,labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=percent,name=expression("Overlap"))+
  scale_color_manual(values=c(#"#fecc5c",#"#fd8d3c",
                              "#e31a1c","black"),labels=c("Before Kraken","After Kraken"),name="")+  
  scale_shape_manual(name="",values=c(0,2,1),labels=c("Before Kraken","After Kraken") )+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
