## Step 0: Set up the working environment ----
setwd("E:/LUAD_CTS/Result4/data/")
library(Seurat); library(CellChat); library(patchwork); options(stringsAsFactors = FALSE)
library(mindr); library(Seurat); library(ggalluvial); library(reshape2) #长、宽数据转换


## Step 1: All Ligand-Receptor Pairs Emitted by Tumor ----
cellchat <- readRDS("output/cellchat.rds")
cells.level <- levels(cellchat@idents)
df.net <- subsetCommunication(cellchat, slot.name = "net", 
                              sources.use = cells.level[14:17], targets.use = cells.level[1:17], 
                              thresh = 0.05)
length(unique(df.net$interaction_name_2)) 
df.net$prob.original <- df.net$prob
df.net$prob <- -1/log(df.net$prob)
df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")

# sort
TME <- c("naive B", "memory_B", "memory_B(stress)", "IgA+ PC", "IgG+ PC", "Igm+ PC(stress)",
         "CD8+ GZMA+ cyto", "CD8+ LAG3+ ex", "CD4+ FOXP3+ treg", "CD4+ TIGIT+ ex",
         "mo_mac", "THBS2+ Myofibroblasts", "tumor ECs","NEF", "Immune", "Cellcycle", "EMT")
Tumor <- c("NEF", "Immune", "Cellcycle", "EMT")
source_target <- paste(rep(Tumor, each = 17), rep(TME, times = 4), sep = " -> ")
df.net$source.target <- factor(df.net$source.target,
                                levels = source_target)
LR_sort <- c("MIF - (CD74+CD44)", "MIF - (CD74+CXCR4)", "MDK - NCL", "MDK - (ITGA4+ITGB1)", "MDK - (ITGA6+ITGB1)", 
             "MDK - SDC1", "MDK - SDC2", "MDK - SDC4", "MDK - LRP1", "C3 - C3AR1", "C3 - (ITGAM+ITGB2)", 
             "C3 - (ITGAX+ITGB2)", "GDF15 - TGFBR2", "GRN - SORT1", "AREG - EGFR", "AREG - (EGFR+ERBB2)", 
             "HBEGF - EGFR", "HBEGF - (EGFR+ERBB2)", "EREG - EGFR", "EREG - (EGFR+ERBB2)", "NAMPT - INSR", 
             "NAMPT - (ITGA5+ITGB1)", "VEGFA - VEGFR1", "VEGFA - VEGFR2", "VEGFB - VEGFR1", "VEGFA - VEGFR1R2", 
             "ANXA1 - FPR1", "SCGB3A2 - MARCO", "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2)", 
             "TGFB1 - (ACVR1B+TGFBR2)", "TGFB2 - (ACVR1B+TGFBR2)", "TGFB1 - (ACVR1+TGFBR1)", "TGFB2 - (ACVR1+TGFBR1)", 
             "ANGPTL4 - (ITGA5+ITGB1)", "ANGPTL4 - CDH5", "ANGPTL4 - CDH11", "ANGPTL4 - SDC1", "ANGPTL4 - SDC2", 
             "ANGPTL4 - SDC3", "ANGPTL4 - SDC4", "LGALS9 - CD44", "LGALS9 - CD45", "LGALS9 - HAVCR2", "WNT7B - (FZD4+LRP5)", 
             "WNT7B - (FZD5+LRP5)", "WNT7B - (FZD6+LRP5)", "WNT7B - (FZD1+LRP6)", "WNT7B - (FZD7+LRP6)", 
             "PDGFA - PDGFRA", "PDGFA - PDGFRB", "EDN1 - EDNRA", "CCL5 - CCR1", "CCL3 - CCR1", "CCL28 - CCR10", 
             "CCL5 - ACKR1", "CXCL1 - ACKR1", "CXCL2 - ACKR1", "CXCL3 - ACKR1", "CXCL8 - ACKR1", "CXCL16 - CXCR6", 
             "HC - C5AR1", "SPP1 - CD44", "SPP1 - (ITGAV+ITGB1)", "SPP1 - (ITGAV+ITGB5)", "SPP1 - (ITGAV+ITGB6)", 
             "SPP1 - (ITGA4+ITGB1)", "SPP1 - (ITGA8+ITGB1)", "SPP1 - (ITGA5+ITGB1)", "PTN - NCL", "PTN - SDC1", 
             "PTN - SDC2", "PTN - SDC3", "PTN - SDC4", "PRSS3 - F2R", "PRSS3 - F2RL1", "PRSS3 - PARD3", "CSF1 - CSF1R",
             "TGFA - EGFR", "TGFA - (EGFR+ERBB2)", "SEMA3B - (NRP1+PLXNA1)", "SEMA3B - (NRP1+PLXNA2)", 
             "SEMA3B - (NRP1+PLXNA3)", "SEMA3C - (NRP1+PLXNA1)", "SEMA3C - (NRP1+PLXNA2)", "SEMA3C - (NRP1+PLXNA3)", 
             "SEMA3B - (NRP2+PLXNA1)", "SEMA3B - (NRP2+PLXNA2)", "SEMA3B - (NRP2+PLXNA3)", "SEMA3C - (NRP2+PLXNA1)", 
             "SEMA3C - (NRP2+PLXNA2)", "SEMA3C - (NRP2+PLXNA3)", "SEMA3C - (NRP1+NRP2)", "SEMA3C - PLXND1")
df <- df.net 
df$interaction_name_2 <- factor(df$interaction_name_2, 
                                levels = rev(LR_sort))
# dotplot
g <- ggplot(df, aes(x = source.target, y = interaction_name_2,
                    color = prob, size = prob)) + 
     geom_point(pch = 16) + 
     theme_linedraw() + 
     theme(panel.grid.major = element_blank()) + 
     theme(axis.text.x = element_text(angle = 90, 
                                      hjust = 1, 
                                      vjust = 0.5), 
           axis.title.x = element_blank(),
           axis.title.y = element_blank()) + 
     scale_x_discrete(position = "bottom") + 
    geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) - 
                                       0.5, 1), lwd = 0.1, colour = "grey90") + 
    geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) - 
                                       0.5, 1), lwd = 0.1, colour = "grey90") + 
    scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                           na.value = "white", 
                           limits = c(quantile(df$prob, 0, na.rm = TRUE), quantile(df$prob, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                           breaks = c(quantile(df$prob, 0, na.rm = TRUE), quantile(df$prob, 1, na.rm = TRUE)), # labels在bar上的位置
                           labels = c("min", "max")) + 
    guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + 
  theme(text = element_text(size = 10), 
        plot.title = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 8), # 设置图例标题的字体大小为8
        legend.text = element_text(size = 6)) + # 设置图例文本的字体大小为6 
  ggtitle( "All Commun") + 
  theme(plot.title = element_text(hjust = 0.5)) # 水平居中标题
g
ggsave("../all_send.pdf", width = 20, height = 20)


## Step 2: All Ligand-Receptor Pairs Received by Tumor----
df.net <- subsetCommunication(cellchat, slot.name = "net", 
                              sources.use = cells.level[1:17], targets.use = cells.level[14:17], 
                              thresh = 0.05)
df.net$prob.original <- df.net$prob
df.net$prob <- -1/log(df.net$prob)
df.net$source.target <- paste(df.net$source, df.net$target, sep = " -> ")

# sort
source_target <- paste(rep(TME,times = 4), rep(Tumor, each = 17), sep = " -> ")
df.net$source.target <- factor(df.net$source.target,
                               levels = source_target)
lr_sort <- c("AREG - EGFR", "AREG - (EGFR+ERBB2)", "MIF - (CD74+CD44)", "IFNG - (IFNGR1+IFNGR2)", 
             "HBEGF - EGFR", "HBEGF - (EGFR+ERBB2)", "LGALS9 - CD44", "SPP1 - CD44", 
             "TNFSF12 - TNFRSF12A", "EREG - EGFR", "EREG - (EGFR+ERBB2)", "MDK - NCL", 
             "MDK - SDC1", "MDK - SDC4", "GDF15 - TGFBR2", "HGF - MET", "PTN - NCL", 
             "PTN - SDC1", "PTN - SDC4", "ANGPTL4 - SDC1", "ANGPTL4 - SDC4", "TGFA - EGFR", 
             "TGFA - (EGFR+ERBB2)", "NAMPT - INSR", "TNF - TNFRSF1A", "SPP1 - (ITGAV+ITGB1)", 
             "SPP1 - (ITGAV+ITGB6)", "WNT5A - FZD5", "WNT2 - (FZD5+LRP5)", 
             "WNT7B - (FZD5+LRP5)", "SEMA3C - (NRP1+PLXNA2)", "SEMA3C - (NRP2+PLXNA2)", 
             "SEMA3F - (NRP2+PLXNA2)", "SEMA3B - (NRP2+PLXNA2)", "SEMA3B - (NRP1+PLXNA2)",
             "MDK - SDC2", "ANGPTL4 - SDC2", "PTN - SDC2", "MIF - (CD74+CXCR4)", 
             "LGALS9 - CD45", "CXCL12 - CXCR4", "GRN - SORT1", "GZMA - PARD3", 
             "SPP1 - (ITGAV+ITGB5)", "MDK - (ITGA6+ITGB1)", "POSTN - (ITGAV+ITGB5)", 
             "SEMA3C - (NRP1+NRP2)", "SEMA3C - PLXND1", "SEMA3C - (NRP1+PLXNA1)", 
             "SEMA3C - (NRP1+PLXNA3)", "SEMA3C - (NRP2+PLXNA1)", "SEMA3C - (NRP2+PLXNA3)", 
             "SEMA3B - (NRP1+PLXNA1)", "SEMA3B - (NRP1+PLXNA3)", "SEMA3B - (NRP2+PLXNA1)", 
             "SEMA3B - (NRP2+PLXNA3)", "SEMA3F - (NRP2+PLXNA1)", "SEMA3F - (NRP2+PLXNA3)", 
             "IGF1 - (ITGA6+ITGB4)", "IGF1 - IGF1R", "WNT5A - FZD6", "WNT2 - (FZD6+LRP5)", 
             "WNT7B - (FZD6+LRP5)", "PRSS3 - PARD3", "TGFB1 - (ACVR1B+TGFBR2)", 
             "TGFB2 - (ACVR1B+TGFBR2)", "TGFB3 - (ACVR1B+TGFBR2)", "PTN - SDC3", 
             "WNT5A - MCAM", "ANGPTL4 - SDC3", "IL1B - (IL1R1+IL1RAP)", "PRSS3 - F2RL1",
             "GZMA - F2RL1")
df.net$interaction_name_2 <- factor(df.net$interaction_name_2,
                                    levels = rev(lr_sort))
# plot
g <- ggplot(df.net, aes(x = source.target, y = interaction_name_2,
                        color = prob, size = prob)) + 
  geom_point(pch = 16) + 
  theme_linedraw() + 
  theme(panel.grid.major = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_x_discrete(position = "bottom") +
  geom_vline(xintercept = seq(1.5, length(unique(df.net$source.target)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90") +
  geom_hline(yintercept = seq(1.5, length(unique(df.net$interaction_name_2)) - 
                                0.5, 1), lwd = 0.1, colour = "grey90")+
  scale_colour_gradientn(colors = colorRampPalette(c("#E6F598", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"))(99), 
                         na.value = "white", 
                         limits = c(quantile(df.net$prob, 0, na.rm = TRUE), quantile(df.net$prob, 1, na.rm = TRUE)), # limits: 定义了颜色映射的范围，即最小值和最大值。 
                         breaks = c(quantile(df.net$prob, 0, na.rm = TRUE), quantile(df.net$prob, 1, na.rm = TRUE)), # labels在bar上的位置
                         labels = c("min", "max")) +  # 设置颜色和对应的colour bar
  guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
g
ggsave("../all_rece.pdf", width = 20, height = 20)
