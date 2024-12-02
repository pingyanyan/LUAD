## Step 0: Set up the working environment ----
setwd("E:/LUAD_CTS/Result4/data/")
library(CellChat); library(patchwork); options(stringsAsFactors = FALSE)
library(mindr); library(Seurat); library(ggalluvial)
cellchat <- readRDS("output/cellchat.rds")
cells.level <- levels(cellchat@idents)


## Step 1: Parameters and data preparation ----
# 1) Set ggplot2 theme
tttheme <- theme(axis.text.x = element_text(angle = 90, 
                                            hjust = 1, 
                                            vjust = 0.5), 
                 axis.title.y = element_text(size = 12),
                 axis.title.x = element_text(size = 12)) + 
                 theme(plot.title = element_text(size=14, 
                                  hjust = 0.5,
                                  vjust = 0.5))
# 2) Set cell ordering
TME <- c("naive B", "memory_B", "memory_B(stress)", "IgA+ PC", "IgG+ PC", "Igm+ PC(stress)",
         "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex", "CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex",
         "mo_mac", "THBS2+ Myofibroblasts", "tumor ECs")
Tumor <- c("NEF", "Immune", "Cellcycle", "EMT")
# 3) Data preparation for plotting the violin plot
w10x <- CreateSeuratObject(cellchat@data.signaling, 
                           meta.data = cellchat@meta)
Idents(w10x) <- w10x$celltype
w10x_tumor <- subset(w10x, idents = Tumor)
Idents(w10x_tumor) <- factor(w10x_tumor$celltype,
                             levels = Tumor)
w10x_TME <- subset(w10x, idents = TME)
Idents(w10x_TME) <- factor(w10x_TME$celltype,
                           levels = TME)
tumor_color <- c("#90bee0", "#4b74b2", "#db3124", "#ffdf92") # 肿瘤细胞的颜色

## Step 2: Common receptor-ligand pairs: MDK-NC ----
pair <- data.frame(interaction_name = "MDK_NCL")
MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                              sources.use = cells.level[14:17], targets.use = cells.level[1:13], 
                              thresh = 0.05)
MN$prob <- -1/log(MN$prob)
MN$source <- factor(MN$source, levels = rev(Tumor))
MN$target <- factor(MN$target, levels = TME)
g1 <- ggplot(MN, 
             aes(x = source, y = target,
                 color = prob, size = prob)) + 
  geom_point(pch = 16) + 
  theme_linedraw() + # 使用"linedraw"主题，该主题可以使图形线条看起来更加手绘
  theme(panel.grid.major = element_blank()) + # 隐藏主要网格线
  scale_x_discrete(position = "bottom") +# 设置x轴在底部 参数"left", "right", "top", "bottom"
  geom_vline(xintercept = seq(1.5, length(unique(MN$source)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(MN$target)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90")+
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(MN$prob, 0, na.rm = TRUE), quantile(MN$prob, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(MN$prob, 0, na.rm = TRUE), quantile(MN$prob, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  # 设置颜色和对应的colour bar
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + 
  guides(size = guide_legend("Commun. Prob.")) + xlab("") + ylab("")+
  labs(title = "MDK - NCL") + tttheme + theme(legend.position = "none") +
  coord_flip()

x <- MN
x$source <- factor(x$source, levels = Tumor)
g2 <- ggplot(x)+
  geom_boxplot(aes(source, prob, color = source)) +
  geom_jitter(aes(source, prob, color = source), width = 0.01) +
  scale_color_manual(values = tumor_color) +
  # theme_classic()+ 
  labs(title = "MDK - NCL") + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  tttheme +  
  theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

g3 <- VlnPlot(w10x_tumor, features = "MDK", pt.size = 0, cols = tumor_color) + 
      theme_classic() + tttheme +  theme(legend.position = "none") + xlab("") +
      ylab("Expression")

g4 <- VlnPlot(w10x_TME, features = "NCL", pt.size = 0) + 
      theme_classic() + tttheme +  theme(legend.position = "none")+ xlab("") +
      ylab("Expression")

layout <- "AAABB" 
g1 + g2 + plot_layout(design = layout)


## Step 3: Plot the early receptor-ligand pairs: SCGB3A2-MARCO ----
pair <- data.frame(interaction_name = "SCGB3A2_MARCO")
MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                          sources.use = cells.level[14:17], targets.use = cells.level[1:13], 
                          thresh = 0.05)
MN$prob <- -1/log(MN$prob)
MN$source <- factor(MN$source, levels = rev(Tumor))
MN$target <- factor(MN$target, levels = TME)
x <- MN
# Since cellcycle is not expressed, add a row with cellcycle prob = NA
x <- rbind(x,x)[1:4, ]; x[4,1] <- "Cellcycle"; x[4,5] <- NA
cell_rece <- c("naive B", "memory_B", "memory_B(stress)", "IgA+ PC", "IgG+ PC", "Igm+ PC(stress)",
               "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex", "CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex",
               "THBS2+ Myofibroblasts", "tumor ECs")
Tumor <- c("NEF", "Immune", "Cellcycle", "EMT")
xx <- x[, 1:5]
x0 <- data.frame(
  source = rep(Tumor, time=12),
  target = rep(cell_rece, each=4),
  ligand = "SCGB3A2",
  receptor = "MARCO",
  prob = NA
)
xall <- rbind(xx, x0)
g1 <- ggplot(xall, 
             aes(x = source, y = target,
                 color = prob, size = prob)) + 
  geom_point(pch = 16) + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank()) + 
  scale_x_discrete(position = "bottom") +
  geom_vline(xintercept = seq(1.5, length(unique(xall$source)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(xall$target)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90")+
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(xall$prob, 0, na.rm = TRUE), quantile(xall$prob, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(xall$prob, 0, na.rm = TRUE), quantile(xall$prob, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + 
  guides(size = guide_legend("Commun. Prob.")) + xlab("") + ylab("")+
  labs(title = "SCGB3A2_MARCO") + tttheme + theme(legend.position = "none") + coord_flip()
g1
gdotplot1 <- g1 
x$source <- factor(x$source, levels = Tumor)
g2 <- ggplot(x, aes(source, prob, fill = source))+
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c("#90bee0", "#4b74b2", "#db3124", "#ffdf92")) +
  labs(title = "SCGB3A2_MARCO") + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  tttheme +  
  theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")
g2
g3 <- VlnPlot(w10x_tumor, features = "SCGB3A2", pt.size = 0, cols = tumor_color) + 
  theme_classic() + tttheme +  theme(legend.position = "none") + xlab("") +
  ylab("Expression")
g4 <- VlnPlot(w10x_TME, features = "MARCO", pt.size = 0) + 
  theme_classic() + tttheme +  theme(legend.position = "none")+ xlab("") +
  ylab("Expression")
layout <- "ABBCCDDDDD"
g1 + g2 + g3 + g4 + plot_layout(design = layout)


## Step 4:Plot the celcycle receptor-ligand pairs: SPP1 − (ITGA5+ITGB1) ----
pair <- data.frame(interaction_name = "SPP1_ITGA5_ITGB1")
MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                          sources.use = cells.level[14:17], targets.use = cells.level[1:13], 
                          thresh = 0.05)
MN$prob <- -1/log(MN$prob)
MN$source <- factor(MN$source, levels = rev(Tumor))
MN$target <- factor(MN$target, levels = TME)
x <- MN 
cell_rece <- c("naive B", "memory_B", "memory_B(stress)", "IgA+ PC", "IgG+ PC", "Igm+ PC(stress)",
               "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex", "CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex",
               "mo_mac", "THBS2+ Myofibroblasts","tumor ECs")
xx <- x[, 1:5]
x0 <- data.frame(
  source = rep(Tumor, time=13),
  target = rep(cell_rece, each=4),
  ligand = "SPP1",
  receptor = "ITGA5_ITGB1",
  prob = NA
)
xall <- rbind(xx, x0)
xall <- xall[!duplicated(xall[,1:2]), ] # Remove the duplicate rows of cellcycle
g5 <- ggplot(xall, 
             aes(x = source, y = target,
                 color = prob, size = prob)) + 
  geom_point(pch = 16) + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank()) + 
  scale_x_discrete(position = "bottom") +
  geom_vline(xintercept = seq(1.5, length(unique(xall$source)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(xall$target)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90")+
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(xall$prob, 0, na.rm = TRUE), quantile(xall$prob, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(xall$prob, 0, na.rm = TRUE), quantile(xall$prob, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  # 设置颜色和对应的colour bar
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + 
  guides(size = guide_legend("Commun. Prob.")) + xlab("") + ylab("")+
  labs(title = "SPP1_ITGA5_ITGB1") + tttheme + theme(legend.position = "none") + coord_flip()
xall$source <- factor(xall$source, levels = Tumor)
g6 <- ggplot(xall, aes(source, prob, color = source))+
  geom_boxplot() +
  geom_jitter(aes(source, prob, color = source), width = 0.01) +
  scale_color_manual(values = c("#db3124", "#ffdf92","#90bee0", "#4b74b2")) +
  # theme_classic()+ 
  labs(title = "SPP1_ITGA5_ITGB1") + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  tttheme +  
  theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")

# vlnplot
g7 <- VlnPlot(w10x_tumor, features = "SPP1", pt.size = 0, cols = tumor_color) + 
  theme_classic() + tttheme +  theme(legend.position = "none") + xlab("") +
  ylab("Expression")
g8 <- VlnPlot(w10x_TME, features = c("ITGA5"), pt.size = 0) + 
  theme_classic() + tttheme +  theme(legend.position = "none")+ xlab("") +
  ylab("Expression") + theme(axis.text.x = element_blank())
g9 <- VlnPlot(w10x_TME, features = c("ITGB1"), pt.size = 0) + 
  theme_classic() + tttheme +  theme(legend.position = "none")+ xlab("") +
  ylab("Expression")


## Step 5:Plot the EMT receptor-ligand pairs: SEMA3C − (NRP2+PLXNA1) ----
pair <- data.frame(interaction_name = "SEMA3C_NRP2_PLXNA1")
MN <- subsetCommunication(cellchat, slot.name = "net", pairLR.use = pair,
                          sources.use = cells.level[14:17], targets.use = cells.level[1:13], 
                          thresh = 0.05)
MN$prob <- -1/log(MN$prob)
MN$source <- factor(MN$source, levels = rev(Tumor))
MN$target <- factor(MN$target, levels = TME)
x <- MN 
x <- rbind(x,x,x,x); x[c(2,3,4),5] <- NA
x[2,1] <- "NEF"; x[3,1] <- "Immune"; x[4,1] <- "Cellcycle"
cell_rece <- c("naive B", "memory_B", "memory_B(stress)", "IgA+ PC", "IgG+ PC", "Igm+ PC(stress)",
               "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex", "CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex",
               "mo_mac", "tumor ECs")
xx <- x[, 1:5]
x0 <- data.frame(
  source = rep(Tumor, time=12),
  target = rep(cell_rece, each=4),
  ligand = "SEMA3C",
  receptor = "NRP2_PLXNA1",
  prob = NA
)
xall <- rbind(xx, x0)
g5 <- ggplot(xall, 
             aes(x = source, y = target,
                 color = prob, size = prob)) + 
  geom_point(pch = 16) + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank()) + 
  scale_x_discrete(position = "bottom") +
  geom_vline(xintercept = seq(1.5, length(unique(xall$source)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(xall$target)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(xall$prob, 0, na.rm = TRUE), quantile(xall$prob, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(xall$prob, 0, na.rm = TRUE), quantile(xall$prob, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + 
  guides(size = guide_legend("Commun. Prob.")) + xlab("") + ylab("")+
  labs(title = "SEMA3C_NRP2_PLXNA1") + tttheme + theme(legend.position = "none") + coord_flip()
x$source <- factor(x$source, levels = Tumor)
g6 <- ggplot(x, aes(source, prob, fill = source))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#90bee0", "#4b74b2", "#db3124", "#ffdf92")) +
  # theme_classic()+ 
  labs(title = "SEMA3C_NRP2_PLXNA1") + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  tttheme +  
  theme(legend.position = "none") + ylab("Commun. Prob") + xlab("")
g7 <- VlnPlot(w10x_tumor, features = "SEMA3C", pt.size = 0, cols = tumor_color) + 
  theme_classic() + tttheme +  theme(legend.position = "none") + xlab("") +
  ylab("Expression")
g8 <- VlnPlot(w10x_TME, features = c("NRP2"), pt.size = 0) + 
  theme_classic() + tttheme +  theme(legend.position = "none")+ xlab("") +
  ylab("Expression") + theme(axis.text.x = element_blank())
g9 <- VlnPlot(w10x_TME, features = c("PLXNA1"), pt.size = 0) + 
  theme_classic() + tttheme +  theme(legend.position = "none")+ xlab("") +
  ylab("Expression")
layout <- "
AAABBCCDDDD
AAABBCCDDDD
EFFFGGGHHHH
EFFFGGGIIII" 
g1 + g2 + g3 +g4 + g5 + g6 +g7+ g8 + g9+plot_layout(design = layout) 




