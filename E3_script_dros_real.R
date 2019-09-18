require(ggplot2); require(scales); require(reshape2)

ds = (read.csv("drosophilaskims.csv"))
names(ds)
kf= (read.csv('/Users/admin/Documents/Skimming/dros_real_data/14_dros_genomes/report_skmer_14dros_real_cleaned_kraken_std_unmasked_0.0_k35/sum_report_skmer_14dros_real_cleaned_kraken_std_unmasked_0.0_k35.csv'))
names(kf)

ds=(merge(merge(ds,kf[,c(2,6)],by.x = "s1", by.y = "species"),kf[,c(2,6)],by.x = "s2", by.y = "species"))


qplot(abs(bk_no_clean_up-fna_dist)/fna_dist-abs(ak_no_clean_up-fna_dist)/fna_dist,data=ds,binwidth=0.02)

qplot(abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist,data=ds,binwidth=0.02)
qplot(fna_dist,abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist,data=ds)
 ggplot(aes(x=fna_dist,y=abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist,
           color=(C_sequences_prct.x+C_sequences_prct.y)/200),data=ds)+
  #geom_violin(aes(group=cut(ds$fna_dist,breaks=c(0,0.03,0.08,0.11,0.126,0.145,0.2))),scale="width")+
  geom_point()+geom_smooth(se=F,method="lm",color="grey50")+
  theme_light()+theme(legend.position = c(.3,.7))+
  scale_color_continuous(name="Filtered",label=percent)+
  scale_y_continuous(labels=percent,name="Decrease in relative error after Kraken")+
  scale_x_continuous(name=("D"),labels=percent)+geom_hline(yintercept = 0, color="red",linetype=2)
ggsave("Drosophila-errreduction.pdf",width=5,height = 4.5)

ggplot(aes(x=fna_dist,y=abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist,
           color=(C_sequences_prct.y)/100),data=ds)+
  #geom_violin(aes(group=cut(ds$fna_dist,breaks=c(0,0.03,0.08,0.11,0.126,0.145,0.2))),scale="width")+
  geom_point()+geom_smooth(se=F,method="lm",color="grey50")+
  theme_light()+theme(legend.position = c(.8,.1))+
  scale_color_continuous(name="Filtered",label=percent)+
  scale_y_continuous(labels=percent,name="Decrease in relative error after Kraken")+
  scale_x_continuous(name=("D"),labels=percent)+geom_hline(yintercept = 0, color="red",linetype=2)+
  facet_wrap(~s1+C_sequences_prct.x)
ggsave("Drosophila-errreduction-all.pdf",width=9,height = 8)

ggplot(aes(x=abs(bk_cleaned-fna_dist)/fna_dist,y=abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist, color=(C_sequences_prct.x+C_sequences_prct.y)/200),
       data=ds)+
  geom_point()+geom_smooth(se=F,method="lm",color="red")+
  theme_light()+
  scale_color_continuous(name="Filtered",label=percent)+
  scale_y_continuous(labels=percent,name="Reduction in error after Kraken")+
  scale_x_continuous(name=("Error before Kraken"),labels=percent)+geom_hline(color="red",linetype=2,yintercept = 0)
ggsave("Drosophila-filtered--error_before-red_in_error.pdf",width=9,height = 8)

ggplot(aes(color=abs(bk_cleaned-fna_dist)/fna_dist,y=abs(bk_cleaned-fna_dist)/fna_dist-abs(ak_cleaned-fna_dist)/fna_dist, 
           x=(C_sequences_prct.x+C_sequences_prct.y)/200),
       data=ds)+
  geom_point()+geom_smooth(se=F,method="lm",color="red")+
  theme_light()+
  theme(legend.position = c(.20,.77))+
  scale_color_continuous(name="Error before Kraken",label=percent)+
  scale_y_continuous(labels=percent,name="Delta relative error after Kraken")+
  scale_x_continuous(name=("Proportion filtered by Kraken"),labels=percent)+geom_hline(color="red",linetype=2,yintercept = 0)
ggsave("Drosophila--proportion_filtered-delta_in_error_updated.pdf",width=5,height = 4)

ds2= ds
ds2$d2 = 0
ds2[as.character(ds2$s1)<as.character(ds2$s2),"d2"] =ds2[as.character(ds2$s1)<as.character(ds2$s2),"bk_cleaned"]
ds2[as.character(ds2$s1)>as.character(ds2$s2),"d2"] =ds2[as.character(ds2$s1)>as.character(ds2$s2),"ak_cleaned"]

ggplot(aes(fill=abs(d2-fna_dist)/fna_dist, 
           x=s1,y=s2),
       data=ds2)+
  geom_tile()+#geom_smooth(se=F,method="lm",color="red")+
  theme_light()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
  theme(axis.title.x=element_blank(), axis.title.y = element_blank())+
  scale_fill_gradient2(name="Relative error",label=percent)#+
  scale_y_continuous(labels=percent,name="Delta relative error after Kraken")+
  scale_x_continuous(name=("Proportion filtered by Kraken"),labels=percent)+geom_hline(color="red",linetype=2,yintercept = 0)
  ggsave("Drosophila_tile_plot.pdf",width=9,height = 8)

qplot(C_sequences_prct.y/100, abs(bk_cleaned-fna_dist)/fna_dist, data=ds)+
  theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
  geom_smooth(method="lm")+
  stat_summary(color="red")+
  scale_y_continuous(name="Error before Kraken",label=percent)+
  scale_x_continuous(name="Proportion filtered by Kraken",label=percent)
ggsave("Drosophila--proportion_filtered-error_before_Kraken.pdf",width=9,height = 8)

qplot(abs(bk_no_clean_up-fna_dist)/fna_dist-abs(bk_cleaned-fna_dist)/fna_dist,data=ds)

