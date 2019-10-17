require(ggplot2); require(scales); require(reshape2); require(Hmisc);
getwd()
setwd("/Users/admin/Documents/Skimming/tree_of_life/dros_contam_test/test_output")
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

write.csv(dh,'my-data-output.csv', row.names = FALSE)
dh = read.csv('my-data-output_true_H.csv',sep=",")
qplot(Overlap/100,error,data=dh,color=Filtered,shape=Filtered,linetype=Filtered)+geom_line()+
  facet_wrap(~percent(round(D,3)),ncol=1)+
  scale_y_continuous(name="Relative error in Skmer distance", labels=percent)+ #, breaks= c(-sqrt(c(0,0.01,0.05,0.2,0.5,1))) ,labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=percent,name=expression("Overlap (H)"))+
  scale_color_manual(values=c("black","#e31a1c"),labels=c("Before Kraken","After Kraken"),name="")+  
  scale_shape_manual(name="",values=c(6,2),labels=c("Before Kraken","After Kraken") )+
  scale_linetype_manual(name="",values=c(2,1),labels=c("Before Kraken","After Kraken"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
ggsave("E2-h.pdf",width = 3,height = 7)

qplot(Overlap/100,sqrt(abs(error))*sign(error),data=dh,color=Filtered,shape=Filtered,linetype=Filtered)+geom_line()+
  facet_wrap(~percent(round(D,3)),ncol=1, scales="free",strip.position = "right")+
  scale_y_continuous(name="Relative error in Skmer distance", breaks=c(-sqrt(c(1, 0.2)),sqrt(c(0,0.2,1,2,5))),labels=function(x) percent(sign(x)*x^2))+
  scale_x_continuous(labels=scales::percent_format(accuracy = 1),name=expression("Overlap (H)"))+
  scale_color_manual(values=c("black","#e31a1c"),labels=c("Before filtering","After filtering"),name="")+  
  scale_shape_manual(name="",values=c(6,2),labels=c("Before filtering","After filtering") )+
  scale_linetype_manual(name="",values=c(2,1),labels=c("Before filtering","After filtering"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom",panel.grid.major.y = element_line(linetype = 1,size=0.2,color="gray"))
ggsave("E2-h-sqrt-corrected.pdf",width = 3,height = 7)




