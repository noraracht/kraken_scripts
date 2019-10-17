require(ggplot2); require(scales); require(reshape2); require(Hmisc); require(RColorBrewer)
getwd()
setwd("/Users/admin/Documents/Skimming/tree_of_life/dros_contam_test")
df_dros_domains = (read.csv("Drosophila_C_to_different_domains.csv"))
colnames(df_dros_domains)

ggplot(df_dros_domains, aes(x=species, y = values/100, fill = legend_val))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+
  theme(plot.margin = margin(0.4, 0.2, -1.0, 0.5, "cm"))+
  theme (legend.position = c(.41,1.1), legend.direction = "horizontal",
         legend.title = element_text(size = 9), 
         legend.text = element_text(size = 9))+
  guides(shape = guide_legend(override.aes = list(size = 0.0001)))+
  guides(color = guide_legend(override.aes = list(size = 0.0001)))+
  theme(aspect.ratio = 0.55)+
  scale_x_discrete(name="")+
  scale_y_continuous(label = percent, name="Classified reads (%)", limits = c(0/100, 16/100))+
  scale_fill_brewer(name = "Classification group",labels = c("archaea", "bacteria", "human", "other", "virus"), palette = "Set2")
ggsave("Classified_at_D2.pdf",width=5,height = 4)  
  

#guides(fill=guide_legend(ncol =2,byrow=FALSE))+
