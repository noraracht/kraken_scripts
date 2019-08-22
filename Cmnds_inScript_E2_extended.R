require(ggplot2); require(scales); require(reshape2)
r=read.csv("export_df_all_data_extended.csv")
names(r)
options(digits=2)


qplot(FPR,recall,shape=as.factor(-k),color=as.factor(bin),group=bin,data=r[r$FPR<0.2&r$confidence_level<0.4,])+facet_wrap(~confidence_level,nrow=1)+geom_line()+scale_shape_discrete(name="k",labels=rev(unique(r$k)))+scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+theme_bw()+theme(legend.position = "bottom")
ggsave("ROC-peralpha2.pdf",width = 9,height = 3.5)
qplot(FPR,recall,shape=as.factor(k),color=factor(bin),group=bin,data=r[r$FPR<1.1&r$confidence_level<0.4&r$k!=20,])+facet_wrap(~confidence_level,nrow=1)+geom_line()+scale_shape_discrete(name="k")+scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+theme_bw()+theme(legend.position = "bottom")
ggsave("ROC-peralpha-supplementary2.pdf",width = 9,height = 3.5)



qplot(FPR,recall,shape=as.factor(-k),color=factor(bin),group=bin,data=r[r$FPR<0.2&r$confidence_level<0.4&r$k!=20,])+
  facet_wrap(~confidence_level,nrow=1)+geom_line()+scale_shape_discrete(name="k",labels=c("35","32","29","28","27","26"))+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-peralpha_extended.pdf",width = 9,height = 3.5)

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(k),group=k,data=r[r$FPR<0.2&r$confidence_level<0.4& !r$k %in% c(20,27,29),])+
  facet_wrap(~bin,nrow=2,labeller = function(x) {x$bin=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25");x})+
  scale_shape_manual(name="confidence",values=c(1,0,2,6,5,4))+scale_color_brewer(palette = "Dark2", name="k")+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")+
  scale_x_continuous(breaks=c(0.01,0.09,0.17))
ggsave("ROC-perbin_extended.pdf",width=8,height = 4.8)

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(bin),group=bin,data=r[r$FPR<0.2&r$confidence_level<0.4& !r$k %in% c(20,27,29),])+
  facet_wrap(~k,nrow=1,labeller = label_both)+geom_line()+
  scale_shape_manual(name="confidence",values=c(1,0,2,6,5,4))+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-perk_extended.pdf",width = 9,height = 3.5)



qplot(FPR,recall,shape=as.factor(-k),color=factor(bin),group=bin,data=r[r$FPR<1.1&r$confidence_level<0.4&r$k!=20&r$k!=23,])+
  facet_wrap(~confidence_level,nrow=1)+geom_line()+scale_shape_discrete(name="k",labels=c("35","32","29","28","27","26"))+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-peralpha-full_extended.pdf",width = 9,height = 3.5)

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(k),group=k,data=r[r$FPR<1.2&r$confidence_level<0.4&r$k!=20,])+
  facet_wrap(~bin,nrow=2,labeller = function(x) {x$bin=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25");x})+
  scale_shape_manual(name="confidence",values=c(1,0,2,6,5,4))+scale_color_brewer(palette = "Dark2", name="k")+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")+
  scale_x_continuous(breaks=c(0.2,0.5,0.8))
ggsave("ROC-perbin-full_extended.pdf",width=8,height = 4.8)


qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(bin),group=bin,data=r[r$FPR<1.1&r$confidence_level<0.4&r$k!=20,])+
  facet_wrap(~k,nrow=1)+geom_line()+
  scale_shape_discrete(name="confidence")+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-perk-full_extended.pdf",width = 9,height = 3.5)
