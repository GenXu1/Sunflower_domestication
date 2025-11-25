library(data.table)
library(ape)
library(ggplot2)
library(tidyverse)
library(ggtree)
library(treeio)
d <- fread("all_sample_info.txt", header = TRUE, data.table = FALSE)
d0=d[d[,2]=="Landrace"|d[,2]=="wild",]
d=d[-which(d[,1] %in% d0[,1]),]
d1=d[grep("^HA|^RHA",d[,2]),]
d2=d[-grep("^HA|^RHA",d[,2]),]
d2[,2]="Other"
d=rbind(d0,d1,d2)
 group_info <- d[, c(1, 2)]
colnames(group_info) <- c("sample", "group")

id  <- fread("Sunflower_SAM_HanXRQr2_dp1_miss05_maf001_332samples.mibs.id.txt", header = F,data.table = FALSE)
ibs <- fread("Sunflower_SAM_HanXRQr2_dp1_miss05_maf001_332samples.mibs.txt", header = F,data.table = FALSE)
ibs= as.matrix(ibs)
rownames(ibs) <- colnames(ibs) <- id$V2

dist_matrix <- as.dist(1 - ibs)
nj_tree <- nj(dist_matrix)

nj_tree <- ladderize(nj_tree, right = TRUE)

group_levels <- c("wild","Landrace","HA-Oil","HA-NonOil","HA","RHA-Oil","RHA-NonOil","RHA","Other")
group_info$group <- factor(group_info$group, levels = group_levels)

order_by_group <- group_info %>%
  arrange(group, sample) %>%           # 组间先按 group 排，组内按 sample 名称；也可换成其它指标
  pull(sample)

nj_tree_rot <- rotateConstr(nj_tree, order_by_group)

is_mono <- tapply(group_info$sample, group_info$group, function(tips) {
  tips <- intersect(tips, nj_tree_rot$tip.label)
  if (length(tips) >= 2) is.monophyletic(nj_tree_rot, tips = tips) else NA
})
print(is_mono)

group_colors <- c(
  "wild" = "#A020F080",
  "Landrace" = "chartreuse4",
  "HA-Oil"    = "darkorange1",
  "HA-NonOil"     = "gold",
  "HA" = "#e6ab02",
  "RHA-Oil"     = "darkcyan",
  "RHA-NonOil"     = "cyan1",
  "RHA"    = "cornflowerblue",
  "Other"  = "gray"
)
colours()[c(142,147,91)]

tree_unr  <- ape::unroot(nj_tree_rot)
tree_grp  <- groupOTU(tree_unr, split(group_info$sample, group_info$group))

p <- ggtree(tree_grp, layout="unrooted", branch.length="none",
       aes(color=group), linewidth=0.4) +
  geom_tippoint(aes(color=group), size=1.0, alpha=0.9) +
  scale_color_manual(values=group_colors, na.value="grey70", name="Heterotic Group") +
  theme_tree2() +
  ggtitle("Unrooted NJ tree (branches grouped by heterotic group)")
p
ggsave("IBS_tree_circular_grouped_all_SNPs_unrooted_all_samples.pdf", p, width = 5, height = 5)
