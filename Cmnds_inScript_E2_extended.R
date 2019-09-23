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

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(k),group=k,data=r[r$FPR<0.3&r$confidence_level<0.4& !r$k %in% c(20,27,29),])+
  facet_wrap(~bin,nrow=2,labeller = function(x) {x$bin=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25");x})+
  scale_shape_manual(name="confidence",values=c(1,0,2,6,5,4))+scale_color_brewer(palette = "Dark2", name="k")+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")+
  scale_x_continuous(breaks=c(0.01,0.09,0.17))
ggsave("ROC-perbin_extended.pdf",width=8,height = 4.8)

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(k),group=k,data=r[r$FPR<1.2&r$confidence_level<0.4&r$k!=20,])+
  facet_wrap(~bin,nrow=2,labeller = function(x) {x$bin=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25");x})+
  scale_shape_manual(name="confidence",values=c(1,0,2,6,5,4))+scale_color_brewer(palette = "Dark2", name="k")+
  theme_light() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")+
  scale_x_continuous(breaks=c(0.2,0.5,0.8))
ggsave("ROC-perbin-full_extended.pdf",width=8,height = 4.8)

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(k),group=k,data=r[r$FPR<1.3&r$confidence_level<0.4& !r$k %in% c(20,27,29),])+
  geom_line()+
  facet_wrap(~bin,ncol=6,labeller = function(x) {x$bin=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25");x})+
  scale_shape_manual(name=expression(alpha),values=c(1,0,2,6,5,4))+scale_color_brewer(palette = "Dark2", name="k", direction = -1)+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = c(.84,.23),legend.direction = "horizontal",legend.box.background =  element_rect(linetype = 2))+
  scale_x_continuous(breaks=c(0.2,0.5,0.8))+scale_y_continuous(labels=percent, name="Recall")
ggsave("ROC-perbin_extended-col.pdf",width=8*1.2,height =3*1.2)

qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(bin),group=bin,data=r[r$FPR<0.2&r$confidence_level<0.4& !r$k %in% c(20,27,29),])+
  facet_wrap(~k,nrow=1,labeller = label_both)+geom_line()+
  scale_shape_manual(name="confidence",values=c(1,0,2,6,5,4))+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-perk_extended.pdf",width = 9,height = 3.5)

ggplot(aes(x=as.factor(bin),group=interaction(confidence_level,k),shape=as.factor(confidence_level),color=interaction(k,sep=" / ")),data=r[r$confidence_level %in% c(0,0.05)& r$k %in% c(28,32,35),])+
  geom_line(aes(y=recall,linetype="Recall"))+ geom_point(aes(y=recall,group=1))+
  geom_line(aes(y=FPR,linetype="FPR"))+ geom_point(aes(y=FPR,group=1))+
  scale_shape_manual(name=expression(alpha),values=c(1,0))+scale_linetype_manual(name="",values=c(3,1))+
  scale_x_discrete(name="Distance to the closest match (M)",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1),axis.text.x = element_text(angle=45,hjust = 1,vjust=1), legend.box.background =  element_rect(linetype = 2), legend.position = c(.29,.3),legend.direction = "horizontal", legend.margin = margin(.0,10,.0,5), legend.spacing.y = unit(.5,"pt") )+
  scale_y_continuous(name="Recall",labels=percent,sec.axis = sec_axis(~.*1, name = "FPR", labels = percent, breaks = c(0.007,0.048,0.086,.539)))+
  scale_color_brewer(palette = "Dark2", name=expression(k),labels=function(x) {gsub(pattern = "(.*) / (.*)","\\2 / \\1",sub(pattern = "-","",x))}, direction = -1)
ggsave("Recall-FPR.pdf",width = 5*0.95, height = 4.2*0.95)

ggplot(aes(x=as.factor(bin),group=interaction(k),color=interaction(k,sep=" / ")),data=r[r$confidence_level <0.4 & !r$k %in% c(23,20),])+
  geom_line(aes(y=recall,linetype="Recall"))+ geom_point(aes(y=recall,group=1,shape="Recall"))+
  facet_wrap(~confidence_level,nrow=1)+
  geom_line(aes(y=FPR,linetype="FPR"))+ geom_point(aes(y=FPR,group=1,shape="FPR"))+
  scale_shape_discrete(name="")+scale_linetype_manual(name="",values=c(3,1))+
  scale_x_discrete(name="Distance to the closest match (M)",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_light() +theme(panel.border  = element_rect(fill=NA,size = 1),axis.text.x = element_text(angle=45,hjust = 1,vjust=1), legend.position = "bottom")+
  scale_y_continuous(name="Recall",labels=percent,sec.axis = sec_axis(~.*1, name = "FPR", labels = percent))+
  scale_color_brewer(palette = "Paired", name=expression(k~"/"~alpha),labels=function(x) {gsub(pattern = "(.*) / (.*)","\\2 / \\1",sub(pattern = "-","",x))})
ggsave("Recall-FPR-full.pdf",width = 9*1.2, height = 3.5*1.2)


qplot(FPR,recall,shape=as.factor(-k),color=factor(bin),group=bin,data=r[r$FPR<1.1&r$confidence_level<0.4&r$k!=20&r$k!=23,])+
  facet_wrap(~confidence_level,nrow=1)+geom_line()+scale_shape_discrete(name="k",labels=c("35","32","29","28","27","26"))+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_light() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-peralpha-full_extended.pdf",width = 9,height = 3.5)



qplot(FPR,recall,shape=as.factor(confidence_level),color=factor(bin),group=bin,data=r[r$FPR<1.1&r$confidence_level<0.4&r$k!=20,])+
  facet_wrap(~k,nrow=1)+geom_line()+
  scale_shape_discrete(name="confidence")+
  scale_color_brewer(palette = "Paired", name="dist to closest match",labels=c("0","(0-1]","(1-2]","(2-3]","(3-5]","(5-10]","(10-15]","(15-20]","(20-25]",">25"))+
  theme_light() +theme(panel.border  = element_rect(fill=NA,size = 1), legend.position = "bottom")
ggsave("ROC-perk-full_extended.pdf",width = 9,height = 3.5)
