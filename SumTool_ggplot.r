library(ggplot2)
dt <- read.table("uk10k+wtccc.ldsc.txt", head=TRUE)
trait <- c("CAD", "HT", "T2D", "BD", "CD", "RA", "T1D")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
dt <- dt[sample(1:nrow(dt), nrow(dt)), ]
dt$trait <- factor(dt$trait, levels=trait)
p <- ggplot(data = dt, aes(x = uk10k, y = wtccc1), group = trait) + 
  geom_abline(intercept=0, slope=1, col="red", linetype=2) + 
  geom_point(aes(colour=trait)) + 
  theme_classic() +
  theme(
    text=element_text(size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = .5, face="bold"),
    legend.position = c(0.1, 0.8),
    legend.text = element_text(size = 10),
    legend.key.height=unit(1.2,"line"),
    legend.key.width=unit(1.2,"line"),
    legend.text.align = 0,
    legend.background = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    #axis.line.x  =element_blank(),
    legend.key = element_blank(),
    # axis.title.x = element_blank(),
    # axis.ticks = element_line(size = 0.5),
    # axis.line.y = element_line(size = 1),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_text(angle = 45, vjust = -0.01),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  #coord_cartesian(ylim = c(0, 0.21), expand = T) +
  # scale_y_continuous(breaks=seq(0,0.21,0.03)) +
  scale_color_manual(values = c(cbPalette)) +
  scale_fill_manual(values = c(cbPalette)) +
  xlab("UK10K") + ylab("WTCCC1") +
  ggtitle("LD score")
  ggsave("LDscore.jpg", p, width=5, height=5, dpi = 600)


SNP <- read.table("exclude.snplist")[, 1]
WTCCC_SNP <- read.table("wtccc.snplist")[, 1]
index <- match(SNP, WTCCC_SNP)
# index <- c(1:length(WTCCC_SNP))[-match(SNP, WTCCC_SNP)]
data <- NULL

for(i in c("BD","CAD","RA","CD","T1D","T2D","HT")){
    
  print(i)

  true <- read.table(paste0("/data2/lxl/llyin/wtccc/sumdata/", i, "_allele_updated_corclean_HMP3.beta"),head=TRUE, stringsAsFactors=FALSE)[index, ]
  
  beagle <- read.table(paste0(i, "_beagle.qassoc"),head=TRUE, stringsAsFactors=FALSE)[index, ]
  freq <- read.table(paste0(i, "_beagle.frq"),head=TRUE, stringsAsFactors=FALSE)$MAF
  freq <- freq[index]

  # SImpute <- read.table(paste0(i, "_New_D.txt"), head=TRUE, stringsAsFactors=FALSE)[index, ]
  # SImpute_40 <- read.table(paste0(i, "_New_D_LD.txt"), head=TRUE, stringsAsFactors=FALSE)[index, ]
  SImpute_p <- read.table(paste0(i, "_phased_New_D_0.01.txt"), head=TRUE, stringsAsFactors=FALSE)[index, ]
  # SImpute_p_40 <- read.table(paste0(i, "_phased_New_D_LD.txt"), head=TRUE, stringsAsFactors=FALSE)[index, ]
  # SImpute_p_40_new <- read.table(paste0(i, "_phased_New_D_40Mb.txt"), head=TRUE, stringsAsFactors=FALSE)[index, ]

  data <- rbind(data, data.frame(trait = i, model="Beagle", maf=freq, x=true$BETA, y=beagle$BETA, stringsAsFactors=FALSE))
  data <- rbind(data, data.frame(trait = i, model="SImpute", maf=freq, x=true$BETA, y=SImpute_p$BETA, stringsAsFactors=FALSE))
  # data <- data[order(data$maf, decreasing=FALSE), ]
}
data$maf <- ifelse(data$maf > 0.5, 1 - data$maf, data$maf)
trait <- c("CAD", "HT", "T2D", "BD", "CD", "RA", "T1D")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
data <- data[sample(1:nrow(data), nrow(data)), ]
data <- data[order(data$maf, decreasing=FALSE), ]
data$trait <- factor(data$trait, levels=trait)

library(ggplot2)
p <- ggplot(data = data, aes(x = x, y = y, colour = trait)) + 
  facet_wrap(~model, scales="free", ncol = 2) + 
  geom_abline(intercept=0, slope=1, col="red", linetype=2) + 
  geom_point(aes(alpha=I(maf*9/5+0.1))) + 
  theme_classic() +
  theme(
    text=element_text(size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = .5, face="bold"),
    # legend.position = c(0.1, 0.8),
    legend.text = element_text(size = 10),
    legend.key.height=unit(1.2,"line"),
    legend.key.width=unit(1.2,"line"),
    legend.text.align = 0,
    legend.background = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face="bold"),
    panel.border = element_blank(),
    #axis.line.x  =element_blank(),
    legend.key = element_blank(),
    # axis.title.x = element_blank(),
    # axis.ticks = element_line(size = 0.5),
    # axis.line.y = element_line(size = 1),
    #axis.ticks.x = element_blank(),
    axis.title.x = element_text(vjust = -0.1),
    # axis.text.x = element_text(angle = 45, vjust = -1.5)
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_cartesian(xlim = c(-0.3, 0.2), ylim = c(-0.3, 0.2), expand = TRUE) +
  # scale_y_continuous(breaks=seq(0,0.21,0.03)) +
  # scale_color_gradient2(midpoint=2.5, low="darkgreen", mid="yellow", high="red", space ="Lab" ) + 
  # scale_color_manual(values = c(cbPalette)) +
  # scale_fill_manual(values = c(cbPalette)) +
  scale_color_manual(values = c(cbPalette)) +
  scale_fill_manual(values = c(cbPalette)) +
  xlab(expression(paste("TRUE ", italic("BETA"), sep= ""))) + ylab(expression(paste("IMPUTED ", italic("BETA"), sep= "")))
  
  ggsave("Beagle_vs_SImpute.jpg", p, width=8, height=4, dpi = 600)


library(plyr)
library(ggplot2)
data1=read.table("data.txt", head=T)
data <- ddply(data1,
c("type","model","trait"),summarise,
N=sum(!is.na(auc)),
min=min(auc),
max=max(auc),
mean=mean(auc),
sd=sd(auc),
se=sd/sqrt(N)
)

type <- c("Individual-level", "Summary-level")
trait <- c("CAD", "HT", "T2D", "BD", "CD", "RA", "T1D")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
model <- c("LMM", "BayesR", "C+T", "SBLUP", "Ldpred2_Grid", "SBayesR_Chr_Dense", "SBayesR_C+Top50K")
data$trait <- factor(data$trait, levels=trait)
data$model <- factor(data$model, levels=model)
data$type <- factor(data$type, levels=type)

p <- ggplot(data = data, aes(x = trait, y = mean)) +
  geom_errorbar(aes(group=model, ymin = mean, ymax= mean + sd), position=position_dodge(width = 0.8), width=0.3) +
  geom_bar(position=position_dodge(0.8),stat = "identity",aes(color = model, fill = model), width = 0.7) +
  theme_classic() +
  theme(
    text=element_text(size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_blank(),
    # legend.position = c(0.1, 0.84),
    legend.text = element_text(size = 7),
    legend.key.height=unit(0.8,"line"),
    legend.key.width=unit(0.8,"line"),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.text.align = 0,
    legend.background = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face="bold"),
    panel.border = element_blank(),
    #axis.line.x  =element_blank(),
    legend.key = element_blank(),
    axis.title.x = element_blank(),
    # axis.ticks = element_line(size = 0.5),
    # axis.line.y = element_line(size = 1),
    #axis.ticks.x = element_blank(),
    # axis.title.x = element_text(vjust = -0.1),
    # axis.text.x = element_text(angle = 45, vjust = -1.5)
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  ) +
  # coord_cartesian(xlim = c(-0.4, 0.6), ylim = c(-0.4, 0.6), expand = TRUE) +
  scale_y_continuous(expand = c(0, 0.03)) + 
  scale_color_manual(values = c(cbPalette)) +
  scale_fill_manual(values = c(cbPalette)) +
  ggtitle("") + xlab("") + ylab(expression(paste("Prediction Accuracy (", italic(AUROC), ")", sep= "")))
  ggsave("wtccc1_model_compare.jpg", p, width=8, height=4, dpi = 600)
