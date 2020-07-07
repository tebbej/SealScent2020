require(ggplot2)

load("RData/objects/pup_colonies_nmds_scaling.RData")

colony_all_pupsonly_16_17 <- ggplot(data = scent_nmds, aes(MDS1, MDS2, color = colony, shape = colony)) + 
  geom_point(size = 4.5) + 
  scale_shape_manual(values = c(15,20,17,15,17,18),
                     labels = c("FWB", "Johnson cove", "Landing beach", "Main bay", "Natural arch", "SSB")) +
  scale_color_manual(values = c("#D55E00", "#000000", "#E69F00", "#009E73", "#CC79A7", "#0072B2"),
                     labels = c("FWB", "Johnson cove", "Landing beach", "Main bay", "Natural arch", "SSB")) +
  theme_void() + 
  annotate("text", x = 0.6, y = -0.94, label = "2D Stress: 0.24", size = 4) +
  theme(panel.background = element_rect(colour = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "right", #c(0.1,0.87),
        legend.title = element_blank(),
        # legend.key.size = unit(0.5, "cm"),
        # legend.key.width = unit(0.5, "cm"),
        legend.background = element_rect(size = 0.3, linetype = "solid", color = "black"))

colony_all_pupsonly_16_17
