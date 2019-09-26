require(ggplot2); require(scales); require(reshape2); require(Hmisc); require(RColorBrewer)
getwd()
setwd("/Users/admin/Documents/Skimming/tree_of_life/dros_contam_test")
df_dros_domains = (read.csv("Drosophila_C_to_different_domains.csv"))
colnames(df_dros_domains)

ggplot(df_dros_domains, aes(x=species, y = values/100, fill = legend_val))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
  theme(plot.margin = margin(0.4, 0.2, -0.2, 0.7, "cm"))+
  theme (legend.position = c(.5,0.94),legend.direction = "horizontal")+
  scale_x_discrete(name="")+
  scale_y_continuous(label = percent, name="Classified reads (%)", limits = c(0/100, 17/100))+
  scale_fill_brewer(name = "Classification group",labels = c("archaea", "bacteria", "human", "other", "virus"), palette = "Set2")
ggsave("Classified_at_D.pdf",width=9,height = 8)  
  

