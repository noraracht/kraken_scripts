require(ggplot2)
require(scales)

# assume:
#  1- all clamination is unshared

j = function(x,fp=0,cl=0) { (1-x)*(1-cl)/(2-(1-x)*(1-cl)) }; # Jaccard with no correction

#  2- FP is randomly distributed across the genome but the same exact positions will be FP in the two genomes. 
j2 = function(x,fp=0,cl=0,fn=0) { ((1-cl)*(1-x)*(1-fp))/((1-cl)*(1+x)*(1-fp)+2*fn*cl) };                               # Jaccard with correction

#  3- FP is randomly distributed across both genomes. 
j3 = function(x,fp=0,cl=0,fn=0) { ((1-cl)*(1-x)*(1-2*fp+fp^2))/((1-cl)*(1+x)*(1-fp^2)+2*fn*cl) };                               # Jaccard with correction
# 2((1-c_l)(1-f_p)+c_lf_n)-(1-c_l)(1-f_p)^2(1-d)
j3 = function(x,fp=0,cl=0,fn=0) { ((1-cl)*(1-x)*(1-2*fp+fp^2)) / (2*(1-cl)*(1-fp) +2*fn*cl - (1-cl)*(1-fp)^2*(1-x)) };                               # Jaccard with correction

dj = function (j) 1-(2*j/(1+j))^(1/31)
d = function (D) (1-(1-D)^31)

cp=c(0, .05,.2); 
r=c(12:20)/20; # percent claminant
a=expand.grid(tp=r, fp=(0:8)*5/100, d=c(0, 50)/100, cl=cp); 
a = data.frame(a, rbind(
			data.frame(Filtered="RR Assumption", est=apply(a,1, FUN=function(x) {j2(x=x[[3]],cl=x[[4]], fp=x[[2]],fn=(1-x[[1]]))})),
			data.frame(Filtered="AR Assumption", est=apply(a,1, FUN=function(x) {j3(x=x[[3]],cl=x[[4]], fp=x[[2]],fn=(1-x[[1]]))}))
			#data.frame(M="N",                    est=apply(a,1, FUN=function(x)            { j(x=x[[3]],cl=x[[4]], fp=x[[2]],fn=(1-x[[1]]))}))
			),
	       cor=apply(a,1,function(x) j(x[[3]])))
a$D=round(1-(1-a$d)^(1/31),2)

qplot((1-tp),fp,fill=est/cor,label=round(est/cor,digits=2),geom="tile",data=a)+
	geom_text(color="white",size=2.0)+
	facet_grid(cl+d+D~Filtered,scales="free_x",shrink=T,space="free_x",labeller=label_both)+theme_classic()+
	scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"))+scale_y_continuous(labels=percent,name=expression("FP portion of genome ("~f[p]~")"))+scale_fill_continuous(name="estimated/correct J",guide="none")+
	geom_text(aes(x=0.2,y=-0.045,label=round(j(x=d,cl=cl)/j(x=d),2)),color="red",size=3.2)


ggsave("backenvlope-jacard.pdf",width=6,height=9)

qplot((1-tp),fp,fill=est/cor,label=round(est/cor,digits=2),geom="tile",data=a)+
	geom_text(color="white",size=2.2)+
	facet_grid(Filtered~cl+d+D,scales="free_y",shrink=T,space="free_y",labeller=label_both)+theme_classic()+
	scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"))+scale_y_continuous(labels=percent,name=expression("FP portion of genome ("~f[p]~")"))+scale_fill_continuous(name="estimated/correct J",guide="none")+
	geom_text(aes(x=0.2,y=-0.045,label=round(j(x=d,cl=cl)/j(x=d),2)),color="red",size=3.2)

ggsave("backenvlope-jacard-h.pdf",width=15,height=6.7)

cp=c(0.02, 0.05, 0.1, .2, .40);  
r=c(c(5:10)*8/80,0.95,0.99); # percent claminant
a=expand.grid(tp=r, fp=(0:7)*5/100, d=d( c(0.002,0.007,0.02,0.06,0.18)), cl=cp); 
a = data.frame(a, #rbind(
			data.frame(M="Filtered (RR Assumption)", est=apply(a,1, FUN=function(x) {j2(x=x[[3]],cl=x[[4]], fp=x[[2]],fn=(1-x[[1]]))})),
			#data.frame(M="N",                    est=apply(a,1, FUN=function(x)            { j(x=x[[3]],cl=x[[4]], fp=x[[2]],fn=(1-x[[1]]))}))
		 #	),
	       cor=apply(a,1,function(x) j(x[[3]])))
a$D=round(1-(1-a$d)^(1/31),3)
qplot((1-tp),fp,fill=-log(est/cor),label=round(est/cor,digits=2),geom="tile",data=a[xor(a$tp==0 , a$M!="N"),])+ geom_text(color="white",size=2.0)+ facet_grid(D~cl,scales="free_x",shrink=T,space="free_x",labeller=label_both)+theme_classic()+ scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"))+scale_y_continuous(labels=percent,name=expression("FP portion of genome ("~f[p]~")"))+scale_fill_continuous(name="estimated/correct J",guide="none")+geom_text(aes(x=0.3,y=-0.045,label=round(j(x=d,cl=cl)/j(x=d),2),color=j(x=d,cl=cl)/j(x=d)),size=3.1)+scale_color_continuous(name="estimated/correct J",guide="none")

qplot((1-tp),fp,fill=-(abs((cor-est)/cor))^1/10,label=percent((cor-est)/cor),geom="tile",data=a[xor(a$tp==0 , a$M!="N"),])+ geom_text(color="white",size=2.0)+
  facet_grid(D~cl,scales="free_x",shrink=T,space="free_x",labeller=label_both)+theme_classic()+
  scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"))+scale_y_continuous(labels=percent,name=expression("FP portion of genome ("~f[p]~")"))+
  scale_fill_continuous(name="estimated/correct J",guide="none")+
  geom_text(aes(x=0.3,y=-0.045,label=percent((j(x=d)-j(x=d,cl=cl))/j(x=d)),color=-(j(x=d)-j(x=d,cl=cl))/j(x=d)^1/10),size=3.1)+scale_color_continuous(name="estimated/correct J",guide="none")
ggsave("backenvlope-RR.pdf",width=10,height=9)

qplot((1-tp),y=abs(((est)-(cor))/(cor)),color=as.factor(-fp),geom="line",data=a[xor(a$tp==0 , a$M!="N"),])+ 
  facet_grid(D~cl,scales="free_y",labeller=label_both)+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1))+
  scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"),breaks = c(0.01,0.2,0.4))+
  scale_y_continuous(labels=percent,name="Relative Jaccard error")+
  scale_fill_continuous(name="estimated/correct J",guide="none")+
  geom_hline(aes(yintercept=(j(x=d)-j(x=d,cl=cl))/j(x=d)),color="red",linetype=2)+
  scale_color_brewer(name=expression("FP ("~f[p]~")"),palette = 1,label=function(x) percent(-as.numeric(x)))
ggsave("backenvlope-RR-lines.pdf",width=6.5,height=5.5)

qplot((1-tp),fp,fill=-log10(((dj(est)-dj(cor))/dj(cor))+0.2),label=percent((dj(est)-dj(cor))/dj(cor)),geom="tile",data=a[xor(a$tp==0 , a$M!="N"),])+ geom_text(color="white",size=2.0)+ facet_grid(D~cl,scales="free_x",shrink=T,space="free_x",labeller=label_both)+theme_classic()+ scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"))+scale_y_continuous(labels=percent,name=expression("FP portion of genome ("~f[p]~")"))+scale_fill_continuous(name="estimated/correct J",guide="none")+geom_text(aes(x=0.3,y=-0.045,label=percent((dj(j(x=d,cl=cl))-dj(j(x=d)))/dj(j(x=d))),color=-log10((dj(j(x=d,cl=cl))-dj(j(x=d)))/dj(j(x=d))+0.2)),size=3.1)+scale_color_continuous(name="estimated/correct J",guide="none")
ggsave("backenvlope-RR-d.pdf",width=10,height=9)

qplot((1-tp),y=abs((dj(est)-dj(cor))/dj(cor)),color=as.factor(-fp),geom="line",data=a[xor(a$tp==0 , a$M!="N"),])+ 
  facet_grid(D~cl,scales="free_y",labeller=label_both)+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1))+
  scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"),breaks = c(0.01,0.2,0.4))+
  scale_y_continuous(labels=percent,name="Relative distance error")+
  scale_fill_continuous(name="estimated/correct J",guide="none")+
  geom_hline(aes(yintercept=((dj(j(x=d,cl=cl))-dj(j(x=d)))/dj(j(x=d)))),color=2,linetype=2)+
  scale_color_brewer(name=expression("FP ("~f[p]~")"),  palette = 1,label=function(x) percent(-as.numeric(x)))
ggsave("backenvlope-RR-d-lines.pdf",width=6.5,height=5.5)

qplot((1-tp),y=abs((dj(est)-dj(cor))/dj(cor)),color=as.factor(-fp),linetype=as.factor(cl),group=interaction(fp,cl),geom="line",
      data=a[xor(a$tp==0 , a$M!="N") &a$tp<1 & a$cl %in%c(0.02, 0.05,0.2,0.4) ,])+ 
  facet_grid(.~D,labeller=label_both)+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1),legend.position="bottom")+
  scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"),breaks = c(0.01,0.2,0.4))+
  scale_y_log10(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_fill_continuous(name="estimated/correct J",guide="none")+
  geom_hline(aes(yintercept=((dj(j(x=d,cl=cl))-dj(j(x=d)))/dj(j(x=d))),linetype=as.factor(cl)),color=2)+
  scale_linetype_manual(name=expression(c[l]),values=c(2,3,4,5))+
  scale_color_brewer(name=expression("FP ("~f[p]~")"),  palette = 1,label=function(x) percent(-as.numeric(x)))
ggsave("backenvlope-RR-d-lines-compact.pdf",width=11,height = 3.4)

qplot(h, (estd-cord)/cord,color=factor((cl),levels=sort(cp,decreasing = T),labels = percent(sort(cp,decreasing = T))),
      group=cl, data=a,geom="line") +facet_wrap(~D,nrow=1,labeller = label_both)+
  theme_classic()+ scale_color_brewer(palette = "Spectral",name=expression(c[l])) + coord_cartesian(ylim=c(1.5,-1))+
  scale_y_continuous(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_x_continuous(labels=percent,name="Contaminant Jaccard (H)")+
  theme(legend.position=c(.86,.83),legend.direction = "horizontal",panel.border  = element_rect(fill=NA,size = 1),
        panel.grid.major.y = element_line(size = 0.08,color = "gray70"))



qplot((1-tp),y=abs((dj(est)-dj(cor))/dj(cor)),color=as.factor(-fp),geom="line",data=a[xor(a$tp==0 , a$M!="N"),])+ 
  facet_grid(D~cl,scales="free_y",labeller=label_both)+
  theme_classic() +theme(panel.border  = element_rect(fill=NA,size = 1))+
  scale_x_continuous(labels=percent,name=expression("FN rate ("~f[n]~")"),breaks = c(0.01,0.2,0.4))+
  scale_y_log10(labels=percent,name="Relative distance error")+
  scale_fill_continuous(name="estimated/correct J",guide="none")+
  geom_hline(aes(yintercept=((dj(j(x=d,cl=cl))-dj(j(x=d)))/dj(j(x=d)))),color=2,linetype=2)+
  scale_color_brewer(name=expression("FP ("~f[p]~")"),  palette = 1,label=function(x) percent(-as.numeric(x)))
ggsave("backenvlope-RR-d-lines-log.pdf",width=6.5,height=5)


cp=c((1:80)/100);
Ds = c(0.001,c(1:20)/100)
a = expand.grid(d=Ds, cl=cp)
a = data.frame(a, rbind(
			data.frame(M="Unfiltered", est=apply(a,1, FUN=function(x) {j(x=d(x[[1]]),cl=x[[2]])}))
		  ),
	       cor=apply(a,1,function(x) j(d(x[[1]]))))
a$estd = dj (a$est)
a$cord = dj (a$cor)

qplot(cl, est/cor,color=sqrt(d),group=d, data=a,geom="line") +theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")"))+scale_y_continuous(name=expression(frac("estimated","true")~Jaccard))+scale_color_continuous(name=element_blank(),breaks=sqrt(c(0.2,0.12,0.001,0.02,0.06)),labels=function(x) {paste("D=",x^2,"; d= ",round(d(x^2),2),sep="")})+theme(legend.position=c(.7,.8))
ggsave("j-vs-cont.pdf",width=5.5,height=4)

qplot(cl, (cor-est)/cor,color=sqrt(d),group=d, data=a,geom="line") +theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")"))+scale_y_log10(name=expression(frac("true"-"estimated","true")~Jaccard))+scale_color_continuous(name=element_blank(),breaks=sqrt(c(0.2,0.12,0.001,0.02,0.06)),labels=function(x) {paste("D=",round(x^2,3),"; d= ",round(d(x^2),2),sep="")})+theme(legend.position=c(.7,.28))+geom_hline(yintercept=0.1,linetype=3,color="gray50")
ggsave("j-vs-cont-rel.pdf",width=5.5,height=4)

qplot(cl, (estd-cord)/cord,color=sqrt(d),group=d, data=a[a$cord!=0,],geom="line")+theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")"))+scale_y_log10(name=expression(frac("estimated"-"true","true")~" genomic distance (D)"),breaks=c(0.001,0.01, 0.1,1,10),labels=scales::comma)+scale_color_continuous(name=element_blank(),labels=function(x) {paste("D=",round(x^2,3),"; d= ",round(d(x^2),2),sep="")},breaks=sqrt(c(0.2,0.12,0.001,0.02,0.06)))+theme(legend.position=c(.87,.22))+geom_text(aes(label=d,color=sqrt(d)),data=a[a$cl==max(a$cl) & a$d %in%c(0.001,0.01,0.02,0.03,0.05,0.08,0.12,0.2),],nudge_x=0.028,size=3.2)+geom_hline(yintercept=0.1,linetype=3,color="gray50")
ggsave("d-vs-cont.pdf",width=5.5,height=4)




jh = function(d,h=0,cl=0) {c=2*cl/(1-cl); hx=h*c/(1+h); ((1-d)+hx)/(1+d+c-hx) }; # Jaccard with no correction


cp=c(0.01,0.02,0.05,(1:5)/10);
Ds = c(0.002,0.007,0.02,0.06,0.18)
Hs = (0:100)/100
a = expand.grid(D=Ds, cl=cp, h=Hs)
a = data.frame(a, rbind(
  data.frame(M="Unfiltered", est=apply(a,1, FUN=function(x) {jh(d=d(x[[1]]),h=x[[3]],cl=x[[2]])}))
),
cor=apply(a,1,function(x) j(d(x[[1]]))))
a$estd = dj (a$est)
a$cord = dj (a$cor)
qplot(cl, abs(cor-est)/cor,color=sqrt(h),group=h, data=a,geom="line") +facet_wrap(~D,scales="free_y")+
 theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")")) +
  scale_y_log10(name=expression(frac("true"-"estimated","true")~Jaccard))#+
  scale_color_continuous(name=element_blank(),breaks=sqrt(c(0.2,0.12,0.001,0.02,0.06)),labels=function(x) {paste("D=",round(x^2,3),"; d= ",round(d(x^2),2),sep="")})+theme(legend.position=c(.7,.28))+geom_hline(yintercept=0.1,linetype=3,color="gray50")

p= qplot(cl, (estd-cord)/cord,color=(h),group=h, data=a,geom="line") +
    theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")")) + coord_cartesian(ylim=c(1.5,-1)) +
  scale_y_continuous(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_color_continuous(name="H",labels=percent)+
  geom_hline(yintercept=0,linetype=1,color="red")

p +  facet_wrap(~D,nrow=1,labeller = label_both)+ 
  theme(legend.position=c(.89,.78),legend.direction = "horizontal",panel.border  = element_rect(fill=NA,size = 1),
          panel.grid.major.y = element_line(size = 0.05,color = "gray80"))
ggsave("j-h-cl.pdf",width=11,height = 2.8)  

p +  facet_wrap(~D,nrow=2,labeller = label_both)+ coord_cartesian(xlim=c(0,0.32),ylim=c(1.5,-1)) +
  theme(legend.position=c(.89,.18),legend.direction = "vertical",panel.border  = element_rect(fill=NA,size = 1),
        panel.grid.major.y = element_line(size = 0.05,color = "gray80"))
ggsave("j-h-cl2.pdf",width=5.6,height = 5.6)  

qplot(h, (estd-cord)/cord,color=factor((cl),levels=sort(cp,decreasing = T),labels = percent(sort(cp,decreasing = T))),
      group=cl, data=a,geom="line") +facet_wrap(~D,nrow=1,labeller = label_both)+
    theme_classic()+ scale_color_brewer(palette = "Spectral",name=expression(c[l])) + coord_cartesian(ylim=c(1.5,-1))+
    scale_y_continuous(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_x_continuous(labels=percent,name="Contaminant Jaccard (H)")+
    theme(legend.position=c(.86,.83),legend.direction = "horizontal",panel.border  = element_rect(fill=NA,size = 1),
          panel.grid.major.y = element_line(size = 0.08,color = "gray70"))
ggsave("j-h.pdf",width=11,height = 2.8)  
  

correctd = function(c1=0.15,c2=0.1,h=0.05,dold=0.08886,k=31) {
  a=c1/(1-c1)+c2/(1-c2); 
  #print(c( (1-dold)^k,(1+a/2),(1-dold)^k*(1+a/2),(a*h)/(1+h)))
  1- ( (1-dold)^k*(1+a/2) - (a*h)/(1+h) )^(1/k)
}

qplot(h, (correctd(c1=cl/y,c2=cl/y,h=h,dold=estd)-cord)/cord, linetype=as.factor(y),
      color=factor((cl),levels=sort(cp,decreasing = T),labels = percent(sort(cp,decreasing = T))),
      group=interaction(cl,y), data=merge(a,c(0.95,0.99,1.01,1.05)), geom="line") +facet_wrap(~D,nrow=1,labeller = label_both)+
  theme_classic()+ 
  scale_color_brewer(palette = "Spectral",name=expression(c[l])) + coord_cartesian(ylim=c(1,-1))+
  scale_y_continuous(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_x_continuous(labels=percent,name="Contaminant Jaccard")+
  scale_linetype_manual(name="True/Estimated\ncontamination",values=c(3,4,2,5))+
  theme(legend.position="bottom",legend.direction = "horizontal",panel.border  = element_rect(fill=NA,size = 1),
        panel.grid.major.y = element_line(size = 0.08,color = "gray70"))


qplot(cl, (correctd(c1=cl/y,c2=cl/y,h=h,dold=estd)-cord)/cord,color=(h),group=interaction(h,y), linetype=as.factor(y),
      data=merge(a[(a$h*100) %%10==0,],c(0.95,0.99,1.01,1.05)), geom="line") +
  theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")")) + 
  coord_cartesian(ylim=c(1,-1)) +
  scale_y_continuous(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_color_continuous(name="H",labels=percent,guide="legend")+
  geom_hline(yintercept=0,linetype=1,color="red")+ 
  facet_wrap(~D,nrow=1,labeller = label_both)+ 
  scale_linetype_manual(name="True/Estimated\ncontamination",values=c(3,4,2,5))+
  theme(legend.position="bottom",legend.direction = "horizontal",panel.border  = element_rect(fill=NA,size = 1),
        panel.grid.major.y = element_line(size = 0.05,color = "gray80"))
ggsave("filter-free-sensitivity.pdf",width=11,height = 3.4)  

qplot(cl, (correctd(c1=cl,c2=cl,h=h/y,dold=estd)-cord)/cord,color=(h),group=interaction(h,y), linetype=as.factor(y),
      data=merge(a[(a$h*100) %%10==0,],c(0.95,0.99,1.01,1.05)), geom="line") +
  theme_classic()+ scale_x_continuous(labels=percent,name=expression("Contamination level ("~c[l]~")")) + 
  coord_cartesian(ylim=c(1,-1)) +
  scale_y_continuous(name=expression(frac("estimated"-"true","true")~D),label=percent)+
  scale_color_continuous(name="H",labels=percent,guide="legend")+
  geom_hline(yintercept=0,linetype=1,color="red")+ 
  facet_wrap(~D,nrow=1,labeller = label_both)+ 
  scale_linetype_manual(name="True/Estimated\ncontaminant Jaccard (H)",values=c(3,4,2,5))+
  theme(legend.position="bottom",legend.direction = "horizontal",panel.border  = element_rect(fill=NA,size = 1),
        panel.grid.major.y = element_line(size = 0.05,color = "gray80"))
ggsave("filter-free-sensitivity-2.pdf",width=11,height = 3.4)  


# correctd2 = function(c1=0.15,c2=0.1,cc=1,h=0.05,dold=0.08886,k=31) {
#   a=c1/(1-c1)+c2/(1-c2); 
#   1- (((1-dold)^k)/cc*(1+a/2)-(a*h)/(1+h))^(1/k)
# }


