## Effect sizes for old and new season data for PERMANOVA on colony membershipt and family ID effect for Antarctic fur seal
## chemical data

## required packages
library(ggplot2)
library(ggbeeswarm)

## Load data and assign data to data.frames
load("RData/objects/R2_initial_season_btrap.RData")

old_season_colony <- paov_r2_results[[1]][[2]]
old_season_family <- paov_r2_results[[1]][[3]]

load("RData/objects/R2_replication_season_btrap.RData")
new_season_colony <- paov_r2_results[[1]][[2]]
new_season_family <- paov_r2_results[[1]][[3]]

MP_effectsize <- c(old_season_colony, new_season_colony, old_season_family, new_season_family)

MP_effectsize.groups <- c(rep("Colony S1", 5000),
                          rep("Colony S2", 5000),
                          rep("Family S1", 5000), 
                          rep("Family S2", 5000))

MP_effectsize.df <- data.frame(btrap_combined_results = MP_effectsize,
                               btrap_subset_groups = MP_effectsize.groups)
save(MP_effectsize.df, file = "RData/objects/effect_size_df.RData")

# point estimates for PERMANOVA on non-bootstrapped (original) data
point_estimate <- c(mean(old_season_colony), mean(new_season_colony), mean(old_season_family), mean(new_season_family))
# point estimate groups for reasons of comprehensibility
point_estimate_groups <- c("Colony S1", "Colony S2", "Family S1", "Family S2")

# data is a dataframe with two columns: one with the values (adj_vals) and one which specifies to which species each value belongs to
MP_effectsize_gg <- ggplot(MP_effectsize.df, aes(y = btrap_combined_results, 
                                                 x = btrap_subset_groups, 
                                                 color = btrap_subset_groups)) + 
  # this arranges the points according to their density
  geom_quasirandom(alpha = 0.06, size = 3, width = 0.3, bandwidth = 1) + #width = 0.47, bandwidth = 2.5
  scale_color_manual(values = c("#E69F00" ,"#E69F00" ,"#CC79A7", "#CC79A7")) +
  # makes the boxplots 
  geom_boxplot(width = 0.35, outlier.shape = NA, color = "white", alpha = 0.1, lwd=0.8) +
  annotate("point", x = 1, y = point_estimate[4], colour = "#000000", fill = "#CCCCCC", size = 2, shape = 21) + 
  annotate("point", x = 2, y = point_estimate[3], colour = "#000000", fill = "#CCCCCC", size = 2, shape = 21) +
  annotate("point", x = 3, y = point_estimate[2], colour = "#000000", fill = "#CCCCCC", size = 2, shape = 21) +
  annotate("point", x = 4, y = point_estimate[1], colour = "#000000", fill = "#CCCCCC", size = 2, shape = 21) +
  # plots the point for the mode
  # this is a possible theme of the plot, there are many others
  theme_classic() +
  # changes the labels on the x axis
  scale_y_continuous(limits = c(-0.01 ,0.25),
                     breaks = seq(0, 0.25, 0.05)) +
  scale_x_discrete(labels = c("Family S2" = "Mother-offspring similarity\nreplication study",
                              "Family S1" = "Mother-offspring similarity\noriginal study",
                              "Colony S2" = "Colony membership\nreplication study",
                              "Colony S1" = "Colony membership\noriginal study"),
                   limits = c("Family S2",
                              "Family S1",
                              "Colony S2", 
                              "Colony S1")) +
  # geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("") +
  # label for y axis
  ylab("Explained variation [R²]") +
  # flips plot so everything is horizontal
  coord_flip() +
  theme(panel.background = element_rect(colour = "black", size = 1.25, fill = NA),
        axis.text = element_text(colour = "black"),
        text = element_text(size = 15),
        legend.position = "none")


MP_effectsize_gg
