require(ggplot2)
load("RData/objects/mom_pup_nmds_scaling.RData")

## colony membership plot
mp_colony_gg <- ggplot(data = scent_nmds) + 
  geom_point(size = 4.5, aes(MDS1, MDS2, color = BeachAge, shape = BeachAge)) + 
  scale_shape_manual(values = c(19, 1, 19, 1), 
                     labels = c("FWB mothers ", "FWB pups ", 
                                "SSB mothers ", "SSB pups ")) +
  scale_color_manual(values = c("#D55E00", "#D55E00", "#56B4E9", "#56B4E9"), 
                     labels = c("FWB mothers ", "FWB pups ", 
                                "SSB mothers ", "SSB pups ")) +
  theme_void() + 
  ylim(-0.75,1.1) +
  annotate("text", x = 0.64, y = 1.1, label = "A", size = 5) +
  annotate("text", x = 0.47, y = -0.74, label = "2D Stress: 0.23", size = 4) +
  theme(panel.background = element_rect(colour = "black", size = 1, fill = NA),
        aspect.ratio = 1,
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(size = 0.3, linetype = "solid", color = "black")) 
# call colony membership plot
mp_colony_gg

## mother-offspring similarity plot
# create color palette for the plot
clr <- c("#D55E00", "red", "#56B4E9", "#009E73","#000000", "#CC79A7") 

# assign pch values for plotting
shp <- c(0,1,2,7,10,5,6,18,16,17,15) 

# create unique color-pch pairs
color_shape_pairs <- crossing(clr,shp) 

# randomly sample 50 unique pairs (sample without replacement)
set.seed(123) # always get same pairs in a run
color_shape_pairs <- color_shape_pairs[sample(nrow(color_shape_pairs), 50),] 

# assign new dataframes to transform scent_nmds$clr & shp with the unique values we created
color_shape_pairs_plot <- rbind(color_shape_pairs[1:25,],color_shape_pairs[1:7,]
                                ,color_shape_pairs[7,],  color_shape_pairs[8:25,], 
                                color_shape_pairs[26:50,], color_shape_pairs[26:50,])
scent_nmds$clr <- as.factor(color_shape_pairs_plot$clr)
scent_nmds$shp <- as.factor(color_shape_pairs_plot$shp)

# call family plot
mp_family_gg <- ggplot(data = scent_nmds,aes(MDS1,MDS2, color = clr, shape = shp)) + 
  geom_point(size = 4.5) +
  scale_shape_manual(values = as.numeric(levels(scent_nmds$shp))) +
  theme_void() + 
  ylim(-0.75,1.1) +
  scale_color_manual(values = levels(scent_nmds$clr)) +
  annotate("text", x = 0.64, y = 1.1, label = "B", size = 5) +
  annotate("text", x = 0.48, y = -0.74, label = "2D Stress: 0.23", size = 4) +
  theme(panel.background = element_rect(colour = "black", size = 1,
                                        fill = NA), aspect.ratio = 1, 
        legend.position = "none") 
mp_family_gg