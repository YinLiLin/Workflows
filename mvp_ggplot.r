library(ggplot2)
data=read.table("data.txt", sep="\t", head=T)
data$Model <- factor(data$Model,levels=c("GLM","MLM","FarmCPU"))
data$Software <- factor(data$Software,levels=c("PLINK v1.9","PLINK v2.0(64threads)","GEMMA","FarmCPU_pkg","rMVP(1thread)","rMVP(4threads)","rMVP(16threads)","rMVP(64threads)"))
#data$Pop_SNPs <- factor(data$Pop_SNPs,levels=c("N1_M1", "N2_M2", "N4_M4", "N8_M8", "N16_M16"))
#p <- ggplot(data=data, aes(as.factor(Pop_SNPs), as.numeric(Time)/3600,colour=Software, group=Software))
p <- ggplot(data=data, aes(Num, as.numeric(Time)/3600, group=Software))+facet_wrap(.~Model,scale="free")
p <- p+geom_line(aes(colour=Software, linetype=Software),size=0.6)+geom_point(aes(colour=Software), size=1, show.legend=FALSE)
p <- p + theme(
legend.position=c(0.765,0.7),
legend.spacing.x = unit(0, 'cm'),
legend.key.height = unit(0.43, "cm"),
legend.key.width = unit(0.7, "cm"),
axis.text=element_text(face=2,size=8),
axis.title=element_text(face=2,size=12),
strip.text.x=element_text(face=2,size=12),
legend.title=element_blank(),
legend.key = element_rect(fill = NA, colour = NA),
legend.background = element_rect(colour = NA, fill=NA),
legend.text=element_text(face=2, size=8),
# panel.grid.major.x=element_blank(),
panel.grid.minor.x=element_blank()
)
p <- p + xlab("Data units")+ylab("Computing Time(h)")+scale_x_continuous(breaks=c(1,16,64,256),labels=c(1,4,8,16))
p <- p + scale_colour_manual(
	values=c(RColorBrewer::brewer.pal(5,"Set1")[1], RColorBrewer::brewer.pal(5,"Set1")[5],
	RColorBrewer::brewer.pal(5,"Set1")[4],
	RColorBrewer::brewer.pal(5,"Set1")[2],
	"#B8DE29FF","#73D055FF","#3CBB75FF","#1F968BFF"
)) + scale_linetype_manual(values = c(rep(1,4), rep(6,4)))
ggsave("Parallel.tiff", p, width=11, height=3, dpi = 300)


library(ggplot2)
data=read.table("data.txt", sep="\t", head=T)
data$model <- factor(data$model,levels=c("GLM","MLM","FarmCPU"))

p <- ggplot() + geom_polygon(data=data, aes(x=time/3600, y=log2(mem2), alpha=I(rep(0.3, nrow(data))), fill=thread), show.legend=TRUE)+facet_wrap(.~model,scale="free")
p <- p + geom_polygon(data=data, aes(x=time/3600, y=log2(mem1), fill=thread), show.legend=TRUE)
p <- p + theme(
axis.text=element_text(face=2,size=8),
axis.title=element_text(face=2,size=12),
strip.text.x=element_text(face=2,size=12),
legend.title=element_blank(),
legend.position=c(0.78,0.85),
legend.key = element_rect(fill = NA, colour = NA),
legend.background = element_rect(colour = NA, fill=NA),
legend.text=element_text(face=2, size=8),
)+guides(fill=guide_legend(ncol=2))
p <- p + xlab("Computing Time(h)")+ylab("Memory Usage(Mb, log2)")
p <- p+scale_fill_manual(values=RColorBrewer::brewer.pal(5,"Set1")[-1], limits=c("1thread",  "4threads" ,"16threads"  ,"64threads"))
p
ggsave("Memory.tiff", p, width=11, height=3, dpi = 300)

