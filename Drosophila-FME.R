require(ggplot2)
require(scales)
require(ggrepel)

x = read.table("../data/FME_per_species.csv", header=F)
am = read.table("../data/amount.csv",header=F)
ds = merge(x=x,y=am,by="V1")
colnames(ds) = c("species","abs_error", "error_perc", "method", "cleaned_perc")
ds$species <- gsub("_", " ", ds$species)
df = ds[ds$method =="wo_kraken",]
df$err_bk = abs((ds[ds$method=="wo_kraken",]$abs_error - ds[ds$method=="assembly",]$abs_error)/ds[ds$method=="assembly",]$abs_error)
df$aft_bk = ((ds[ds$method=="with_kraken",]$abs_error - ds[ds$method=="assembly",]$abs_error)/ds[ds$method=="assembly",]$abs_error)

percent=function(x) paste0(as.double(x), "%")

ggplot(aes(color=err_bk, y=err_bk-aft_bk, 
           x=cleaned_perc, label=species),
       data=df)+
  geom_hline(color="red",linetype=2,yintercept = 0)+
  geom_smooth(se=F,method="lm",color="red", data = df[!df$species %in% c("Drosophila virilis") ,])+
  geom_point()+ 
  geom_text_repel(size=3, force = 2)+
  theme_light()+
  scale_color_continuous(name="FME error before\n filtering", label=percent)+
  scale_y_continuous(name="Change in relative FME error after filtering", labels=percent)+
  scale_x_continuous(name="Proportion filtered", labels=percent)+
  theme(legend.position = c(0.2, 0.75))

ggsave("Drosophila-FME-error.pdf",width=5,height = 4)

