
tiff("Multi_Species.tiff", res=300, height=300*4, width=300*8)
windowsFonts(TNR=windowsFont("Arial"))

mymodel <- c("LMM","BSLMM","BayesR","KAML")
#mycolor <- RColorBrewer::brewer.pal(5,"Set1")[-1]
mycolor <- c("red3",
"darkgoldenrod1",
"purple",
"dodgerblue1"
)
#mytrait <- c("MFP","MY","SCS","CC","Coat Color","DD","YWK","SSK","GDD")
mytrait <- c("mfp","my","scs","cc","coat colour","dd","ywk","ssk","gdd")
myspecies <- c("Cattle","Horse","Maize")

library(plyr)
data$Accuracy <- round(data$Accuracy, 3)
data$Time <- round(data$Time)
datay <- ddply(data,
c("Species","Trait","Model"),summarise,
N=sum(!is.na(Accuracy)),
min=min(Accuracy),
max=max(Accuracy),
mean=mean(Accuracy),
sd=sd(Accuracy),
se=sd/sqrt(N)
)
datay[is.na(datay)] <- 0

datay$Model <- factor(datay$Model,levels=mymodel)
datay$Trait <- factor(datay$Trait,levels=mytrait)
datay$Species <- factor(datay$Species,levels=myspecies)

p <- ggplot(
data=datay,
aes(x=Trait,y=mean)
)

p <- p+xlab("")+ylab("Prediction Performance (Correlation)")
#p <- p+facet_grid(.~Species,scales="free",space="free")
p <- p+facet_wrap(~Species,scales="free")
p <- p+theme(
axis.title=element_text(family="TNR",face=2,size=12),
strip.text.x=element_text(size=14,face=2,family="TNR"),
legend.title=element_blank(),
#strip.text.y=element_text(size=12,face=2,family="TNR"),
#axis.ticks.x=element_blank(),
#axis.text.x=element_blank(),
axis.text.x=element_text(face=2,size=10.5,family="TNR"),
axis.text.y=element_text(face=2,size=9.5,family="TNR"),
axis.title=element_text(family="TNR",face=1,size=12),
strip.background=element_rect(fill = "gray80", colour = "gray80"),
panel.grid.major.x=element_blank()
)

p <- p+scale_fill_manual(
values=mycolor
)

p <- p+geom_boxplot(aes(fill = Model,
lower = mean - sd, 
upper = mean + sd, 
middle = mean, 
ymin = min, 
ymax = max
),
lwd=0.5,
stat="identity"
)

xx <- data.frame(Species=c("Cattle","Cattle","Maize","Maize"), value=c(1.5, 2.5, 1.5, 2.5))
p <- p + geom_vline(data=xx, aes(xintercept=value), colour="white")

# xx <- cbind(tail(datay, n=12), tail(ggplot_build(p)$data[[1]], n=12))
# p <- p + geom_segment(data=xx, aes(x=xmin, xend=xmax, y=mean, yend=mean), colour=rep(mycolor[c(1,2,4,3)], 3),lwd=1.2)
p
dev.off()


tiff("Multi_Species_Time.tiff", res=300, height=300*4, width=300*8)
windowsFonts(TNR=windowsFont("Arial"))

mymodel <- c("LMM","BSLMM","BayesR","KAML")
#mycolor <- RColorBrewer::brewer.pal(5,"Set1")[-1]
mycolor <- c("red3",
"darkgoldenrod1",
"purple",
"dodgerblue1"
)
#mytrait <- c("MFP","MY","SCS","Coat Color","YWK","SSK","GDD")
mytrait <- c("mfp","my","scs","coat colour","ywk","ssk","gdd")
myspecies <- c("Cattle","Horse","Maize")

library(plyr)
data$Accuracy <- round(data$Accuracy, 3)
data$Time <- log10(round(data$Time))
datay <- ddply(data,
c("Species","Trait","Model"),summarise,
N=sum(!is.na(Time)),
min=min(Time),
max=max(Time),
mean=mean(Time),
sd=sd(Time),
se=sd/sqrt(N)
)
datay[is.na(datay)] <- 0
datay$Model <- factor(datay$Model,levels=mymodel)
datay$Trait <- factor(datay$Trait,levels=mytrait)
datay$Species <- factor(datay$Species,levels=myspecies)
datay <- datay[!is.na(datay[,2]), ]

p <- ggplot(
data=datay,
aes(x=Trait, y=mean)
)

p <- p+xlab("")+ylab("Calculating Time(Log10, s)")
#p <- p+facet_grid(.~Species,scales="free",space="free")
p <- p+theme(
axis.title=element_text(family="TNR",face=2,size=13),
strip.text.x=element_text(size=14,face=2,family="TNR"),
legend.title=element_blank(),
#strip.text.y=element_text(size=12,face=2,family="TNR"),
#axis.ticks.x=element_blank(),
#axis.text.x=element_blank(),
axis.text.x=element_text(face=2,size=11,family="TNR"),
axis.text.y=element_text(face=2,size=9.5,family="TNR"),
axis.title=element_text(family="TNR",face=1,size=12),
strip.background=element_rect(fill = "gray80", colour = "gray80")
)

p <- p+scale_fill_manual(
values=mycolor
)

p <- p+geom_bar(aes(fill=Model), stat = "identity", position=position_dodge())+scale_y_continuous(expand = c(0, 0))+geom_blank(aes(y = 6))
#+geom_errorbar(aes(ymin=mean-se, ymax=mean+se),  width=.5, size=0.3)
p

dev.off()



data=read.delim("KAML_posterior.txt", head=T)
tiff("Human_post.tiff", res=300, height=300*4, width=300*12)
#windowsFonts(TNR=windowsFont("Arial"))

myX <- c("LMM", "Half", "Adaptive")
mycolor <- c("red3",
"palegreen3",
"dodgerblue1"
)
mytrait <- c("cad","ht","t2d","bd","cd","ra","t1d","mfp","my","scs","coat colour","ywk","ssk","gdd")
myspecies <- c("Human","Cattle","Horse","Maize")
library(plyr)
datay <- ddply(data,
c("SP","Ref","trait"),summarise,
N=sum(!is.na(Accuracy)),
min=min(Accuracy),
max=max(Accuracy),
mean=mean(Accuracy),
sd=sd(Accuracy),
se=sd/sqrt(N)
)

datay$SP <- factor(datay$SP,levels=myspecies)
datay$trait <- factor(datay$trait,levels=mytrait)
datay$Ref <- factor(datay$Ref,levels=myX)

p <- ggplot(
data=datay,
aes(x=trait,y=mean)
)

p <- p+xlab("")+ylab("Prediction Performance (AUC/Cor)")
#p <- p+facet_grid(.~SP,scales="free")
p <- p+facet_wrap(~SP,scales="free",nrow=1)
p <- p+theme(
text=element_text(family="Arial"),
axis.title=element_text(face=2,size=13),
strip.text.x=element_text(size=13,face=2),
legend.title=element_blank(),
#strip.text.y=element_text(size=12,face=2),
#axis.ticks.x=element_blank(),
#axis.text.x=element_blank(),
axis.text.x=element_text(face=2,size=10.5),
axis.text.y=element_text(face=2,size=10.5),
strip.background=element_rect(fill = "gray80", colour = "gray80"),
panel.grid.major.x=element_blank()
)

p <- p+scale_fill_manual(
values=mycolor
)

a<-0.7
p <- p+geom_boxplot(aes(fill = Ref,
lower = mean - sd, 
upper = mean + sd, 
middle = mean, 
ymin = min, 
ymax = max,
width=0.7
),
lwd=0.5,
stat="identity"
#position = position_dodge(0.5)
)
xx <- data.frame(SP=c(rep("Human",6),"Cattle","Cattle","Maize","Maize"), value=c(1.5,2.5,3.5,4.5,5.5,6.5,1.5,2.5,1.5,2.5))
p <- p + geom_vline(data=xx, aes(xintercept=value), colour="white")

# convert ggplot object to grob object
gp <- ggplotGrob(p)

# optional: take a look at the grob object's layout
gtable::gtable_show_layout(gp)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x, function(l) length(l$range$range))
x.var[x.var == 1] <- 1.2
# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

# plot result
grid::grid.draw(gp)

dev.off()




data$trait <- factor(data$trait, levels=c("cad","ht","t2d","bd","cd","ra","t1d"))
data$trait <- factor(data$trait, levels=c("mfp","my","scs","coat colour","ywk","ssk","gdd"))
data$model <- factor(data$model, levels=c("pQ","pQ+Kw","pQ+Ks","Kw","Ks")
a=tapply(data$count, data$trait, function(x){x=x[x!=0];y=NULL;for(i in 1:length(x)){y=c(y, ifelse(i==1,x[i]/2,sum(x[1:(i-1)])+x[i]/2))};return(y)})
txt = data.frame(trait=rep(names(a),sapply(a,length)), 
	y=unlist(a), 
	lab=data$count[data$count!=0],
	model=data$model[data$count!=0]
	)

p<-ggplot(data=data, aes(x=trait, y=count, fill=model)) +
  geom_bar(stat="identity")

p <- p + xlab("") + ylab("Count of selected model")
p <- p + theme(
	legend.title=element_blank(),
	legend.key = element_rect(fill = NA, colour = NA),
	legend.text=element_text(face=3),
	axis.text.x=element_text(face=2,size=8),
	axis.text.y=element_text(face=2,size=8),
	axis.title.y=element_text(face=2,size=12),
	axis.title.x=element_blank(),
	strip.text.x = element_text(size = 8, face=2),
	strip.background=element_rect(fill = "gray80", colour = "gray80")
) + scale_y_continuous(expand = c(0, 0))

p <- p + scale_fill_viridis_d(direction=-1, limits=c("pQ","pQ+Kw","pQ+Ks","Kw","Ks"))
p <- p + geom_text(data=txt, aes( x=trait, y=y, label=lab), vjust=0.5, 
            color="white", size=3.5, fontface=2)

# ggsave("human_model.jpg", p, width=5, height=5, dpi = 600)
ggsave("species_model.jpg", p, width=5, height=5, dpi = 600)




