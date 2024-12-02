## Step 0: Set up the working environment
setwd("E:/LUAD_CTS/Result4/data/")
library(Seurat); library(CellChat); library(patchwork); options(stringsAsFactors = FALSE)
library(mindr); library(Seurat); library(ggalluvial); library(reshape2) 


## Step 1: Data preprocessing: Obtain the emission and reception matrices (excluding interactions within tumor cells)
cellchat <- readRDS("cellchat.rds")
tumor_send_net_weight <- cellchat@net$weight[14:17, 1:13] 
tumor_recevier_net_weight <- cellchat@net$weight[1:13, 14:17]
tsnw <- melt(tumor_send_net_weight)
trnw <- melt(tumor_recevier_net_weight)
tall <- rbind(tsnw, trnw) # Both plots share a common legend, so use a single matrix when generating the legend
TME <- c("naive B", "memory_B", "memory_B(stress)", "IgA+ PC", "IgG+ PC", "Igm+ PC(stress)",
         "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex", "CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex",
         "mo_mac", "THBS2+ Myofibroblasts", "tumor ECs")
Tumor <- c("NEF", "Immune", "Cellcycle", "EMT")
tsnw$Var1 <- factor(tsnw$Var1,
                    levels = Tumor)
tsnw$Var2 <- factor(tsnw$Var2,
                    levels = rev(TME))
trnw$Var2 <- factor(trnw$Var2,
                    levels = Tumor)
trnw$Var1 <- factor(trnw$Var1,
                    levels = rev(TME))
save(list=c("tsnw", "trnw"), file="FigureA.RData")


## Step 2: Plot the dotplot
load("FigureA.RData")
tall <- rbind(tsnw, trnw)
g1 <- ggplot(tsnw, 
            aes(x = Var1, y = Var2,
                    color = value, size = value)) + 
  geom_point(pch = 16) + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0, 
                                   vjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_x_discrete(position = "bottom") +
  geom_vline(xintercept = seq(1.5, length(unique(tsnw$Var1)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(tsnw$Var2)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90")+
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  # 设置颜色和对应的colour bar
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))+
  scale_x_discrete(position = "top")
g2 <- ggplot(trnw, 
            aes(x = Var2, y = Var1,
                color = value, size = value)) + 
  geom_point(pch = 16) + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 0, 
                                   vjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_x_discrete(position = "bottom") +
  geom_vline(xintercept = seq(1.5, length(unique(trnw$Var2)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(trnw$Var1)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90")+
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  # 设置颜色和对应的colour bar
  guides(color = guide_colourbar(barwidth = 0.7, title = "Weight"),
         size = guide_legend("Weight")) + 
  scale_x_discrete(position = "top") 
G1 <- g1 + scale_radius(range =  c(2,6), 
                        limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE))) +
                        theme(legend.position = "none")
G2 <- g2 + scale_radius(range =  c(2,6), 
                        limits = c(quantile(tall$value, 0, na.rm = TRUE), quantile(tall$value, 1, na.rm = TRUE)),
                        labels = c("weight > 0.0", "weight > 0.1", "weight > 0.2", "weight > 0.3", "weight > 0.4")  )
G1+G2 
                  
               
