library(tidyverse) #for general data processing
library(brms) #for handling model output
library(tidybayes) #for handling brms output for plotting
library(ggrepel) #for adding labels to PCA plot
library(patchwork) #for stitching together multiple panels
library(ggeffects) #for model predictions
source("ggplot_paper_theme.R") #custom theme script

#load in data------------------------------------------------------------------
#initial trends
init <- read_csv("../Data/main_sampling_clean.csv") %>% 
  #select only data from the initial time point
  filter(sampling_period == 0) %>% 
  #make sure treatment is a factor
  mutate(treatment = factor(treatment, 
                            levels = c("natural", "half", "zero"))) %>% 
  filter(treatment != 'dd' & treatment != 'ft') %>% 
  drop_na(shoot_count, mean_blade_height, above_ground, below_ground, 
          percent_total_n)

#full sampling data
props <- read_csv("../Data/sampling_data_averaged_proportions.csv") %>% 
  #remove treatments from other experiment
  filter(treatment != "dd" & treatment != "ft") %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("natural", "half", "zero"))) %>% 
  #remove all rows with NAs for pca to work
  filter_at(vars(mean_shoot_prop, mean_height_prop, mean_above_prop, 
                 mean_below_prop, mean_percent_total_n_prop),
            all_vars(!is.na(.))) %>% 
  #remove the three outliers that are clearly mistakes based on the
  #shoot count/height comparison to above and belowground (see manuscript for
  #full justification of these omissions)
  filter(patch != "Vancouver" | (patch == "Vancouver" & quadrat != 3)) %>% 
  filter(patch != 'Montreal' | patch == "Montreal" & quadrat != 1) %>% 
  filter(patch != '78S' | patch == '78S' & quadrat != 1)

#responses from sampling data to combine with the growth data
below <- read_csv("../Data/sampling_data_averaged_proportions.csv") %>% 
  #remove treatments from other experiment
  filter(treatment != "dd" & treatment != "ft")  %>%   
  dplyr::select(patch, quadrat, mean_shoot_count_offset:mean_total_n_offset) %>% 
  distinct() %>% 
  mutate(scale_mean_height_offset = scale(mean_height_offset),
         scale_mean_below_offset = scale(mean_below_offset))

#growth experiment data
growths <- read_csv("../Data/seagrass_growth_clean.csv") %>% 
  #remove treatments from other experiment
  filter(treatment != "dd" & treatment != "ft") %>% 
  mutate(scale_num_cut = scale(shoots_cut)) %>% 
  left_join(below, by = c("patch", "quadrat")) %>% 
  mutate(below_per_shoot = mean_below_offset/shoots_found,
         scale_below_per_shoot = scale(below_per_shoot))

#load in all the model objects
mod_growth <- readRDS("../Data/final_growth_rate_model.rds")
mod_init <- readRDS("../Data/initial_PC1_mod.rds")
mod_init_n <- readRDS("../Data/initial_nitrogen_mod.rds")
mod_main_density <- readRDS("../Data/main_model_density.rds")
mod_main_treatment <- readRDS("../Data/main_model_treatment.rds")

#FIGURE 1A - INITIAL PCA-------------------------------------------------
#extract the quantitative metrics
init_resp <- init %>% 
  select(shoot_count, mean_blade_height, above_ground, below_ground) %>% 
  drop_na()
#let's condense these variables into principal components
pca_values <- prcomp(init_resp, center = TRUE, scale = TRUE)
summary(pca_values) 

#turn the points into a df
pca_points <- 
  # first convert the pca results to a tibble
  as_tibble(pca_values$x) %>% 
  # now we'll add the full data
  bind_cols(init)

#and now extract the Eigenvectors
pca_load <- 
  as_tibble(pca_values$rotation, rownames = 'variable') %>% 
  # we can rename the variables so they look nicer on the figure
  mutate(variable = dplyr::recode(variable,
                                  "shoot_count" = "Shoot count", 
                                  "mean_blade_height" = "Mean leaf length", 
                                  "above_ground" = "Above-ground biomass", 
                                  "below_ground" = "Below-ground biomass"),
         #extend arrows so they're more visible but keep relative differences
         #*4 makes it identical to the default autoplot loadings
         segment_end1 = PC1 * 4,
         segment_end2 = PC2 * 4,
         #create points at which to place arrow labels
         label_edge1 = c(2.9,2.5,3.2,3.2),
         label_edge2 = c(-1,2.8,1.1,-2.4))

#and plot!
fig1a <- ggplot(pca_points, aes(x = PC1, y = PC2)) +
  geom_point(size = 0.2) +
  geom_segment(data = pca_load, 
               aes(x = 0, y = 0, 
                   xend = segment_end1,
                   yend = segment_end2),
               arrow = arrow(length = unit(1/4, 'picas')),
               size = 0.2) +
  annotate('text', 
           x = (pca_load$label_edge1), 
           y = (pca_load$label_edge2),
           label = pca_load$variable,
           size = 3) +
  theme_paper_small() +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(margin = margin(t = 0, r = -70, b = 0, l = 0),
                                    size = 10),
        plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12)) +
  labs(x = "PC1 (58.03%)", y = "PC2 (23.12%)",
       title = "A")
fig1a
#ggsave("../Figures/Figure1a.png", fig1a, device = "png",
#       height = 60, width = 80, units = "mm", dpi = 600)

#FIGURE 1B - INITIAL TRENDS COEF------------------------------------------------
#create nice labels for the figure
init_mod_labs <- c("b_scale_sand"="Sand depth",
                   "b_scale_patch_size"="Total patch size",
                   "b_scale_log_grunt_pee_m2"="log(Fish nutrient\ninput)", 
                   "b_scale_log_og_cukes:scale_log_grunt_pee_m2"="log(Sea cucumber\ndensity) * log(Fish\nnutrient input)",                   
                   "b_scale_log_og_cukes"="log(Sea cucumber\ndensity)",
                   "b_scale_log_distance"="log(Distance from reef)")
#custom colours
coef_cols <- c("#B7CDBB", "#537970")


fig1b <- mod_init %>% 
  #gather the posterior draws for all of the coefficient estimates for the fixed
  #effects
  gather_draws(`b_.*`, regex = TRUE) %>% 
  #we don't need the intercept though
  filter(`.variable` != "b_Intercept") %>%
  #reorder the levels so we're plotting from lowest to highest
  mutate(`.variable` = factor(`.variable`, 
                              levels = c("b_scale_sand",
                                         "b_scale_patch_size",
                                         "b_scale_log_grunt_pee_m2",
                                         "b_scale_log_og_cukes:scale_log_grunt_pee_m2",
                                         "b_scale_log_og_cukes",
                                         "b_scale_log_distance"))) %>%
  #now plot!
  ggplot(aes(y = .variable, x = .value, fill = stat(x < 0))) +
  stat_halfeye(size = 0.1) +
  geom_vline(xintercept = 0, size = 0.3) +
  scale_y_discrete(labels = init_mod_labs) +
  scale_fill_manual(values = coef_cols) +
  theme_paper_small() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 1,
                                   size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = 
          element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = -0.5, vjust = 2.5, size = 12)) +
  labs(x = "", y = "",
       title = "B")
fig1b
#ggsave("../Figures/Figure1b.png", fig1b, device = "png",
#       height = 60, width = 80, units = "mm", dpi = 600)

#FIGURE 1C - NITROGEN COEF PLOT-----------------------------------------
#same thing as above with the other model
init_mod_labs_n <- c("b_scale_sand"="Sand depth",
                     "b_scale_patch_size"="Total patch size",
                     "b_scale_log_grunt_pee_m2"="log(Fish nutrient\ninput)", 
                     "b_scale_log_og_cukes:scale_log_grunt_pee_m2"="log(Sea cucumber\ndensity) * log(Fish\nnutrient input)",                   
                     "b_scale_log_og_cukes"="log(Sea cucumber\ndensity)",
                     "b_scale_log_distance"="log(Distance from reef)")
coef_cols <- c("#B7CDBB", "#537970")

fig1c <- mod_init_n %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  filter(`.variable` != "b_Intercept") %>% 
  mutate(`.variable` = factor(`.variable`, 
                              levels = c("b_scale_sand",
                                         "b_scale_patch_size",
                                         "b_scale_log_grunt_pee_m2",
                                         "b_scale_log_og_cukes:scale_log_grunt_pee_m2",
                                         "b_scale_log_og_cukes",
                                         "b_scale_log_distance"))) %>% 
  ggplot(aes(y = .variable, x = .value, fill = stat(x < 0))) +
  stat_halfeye(size = 0.1) +
  geom_vline(xintercept = 0, size = 0.3) +
  scale_y_discrete(labels = init_mod_labs_n) +
  scale_fill_manual(values = coef_cols) +
  theme_paper_small() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 1,
                                   size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = 
          element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = -0.5, vjust = 2.5, size = 12)) +
  labs(x = "", y = "",
       title = "C")
fig1c

#FIGURE 1 COMBINED-------------------------------------------
#combine the three panels with patchwork
fig1 <- fig1a + fig1b + fig1c + plot_layout(design = "1111
                                    1111
                                    2233
                                    2233
                                    ")
fig1
#ggsave("../Figures/Figure1.png", fig1, device = "png",
#       height = 180, width = 180, units = "mm", dpi = 600)

#FIGURE S2------------------------------------------------
#take the same PCA plot from Fig1 and modify the axes to make it easier to 
#merge with the second PCA plot
figS2a <- fig1a +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 10),
        plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12)) +
  xlim(c(-3.5,5))
  

pca_values2 <- prcomp(props[,42:45], center = TRUE, scale = TRUE)
summary(pca_values2) 

#turn the points into a df
pca_points2 <- 
  # first convert the pca results to a tibble
  as_tibble(pca_values2$x) %>% 
  # now we'll add the full data
  bind_cols(props) 

#and now extract the Eigenvectors
pca_load2 <- 
  as_tibble(pca_values2$rotation, rownames = 'variable') %>% 
  # we can rename the variables so they look nicer on the figure
  mutate(variable = dplyr::recode(variable,
                                  "mean_shoot_prop" = "Proportional change\nin shoot count", 
                                  "mean_height_prop" = "Proportional change\nin mean leaf length", 
                                  "mean_above_prop" = "Proportional change in\nabove-ground biomass", 
                                  "mean_below_prop" = "Proportional change in\nbelow-ground biomass", 
                                  "mean_percent_total_n_prop" = "Proportional change in\npercent total nitrogen"))


figS2b <- ggplot(pca_points2, aes(x = PC1, y = PC2)) +
  geom_point(size = 0.2) +
  geom_segment(data = pca_load2, 
               aes(x = 0, y = 0, 
                   xend = PC1*5,
                   yend = PC2*5),
               arrow = arrow(length = unit(1/4, 'picas')),
               size = 0.2) +
  annotate('text', x = (pca_load2$PC1*6), y = (pca_load2$PC2*5.2),
           label = pca_load2$variable,
           size = 2) +
  theme_paper_small() +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 10),
        #plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 10),
        plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12)) +
  labs(x = "PC1 (47.37%%)", y = "PC2 (22.06%)",
       title = "B") +
  xlim(c(-3.5,5))
figS2b

#combine with patchwork
figS2 <- figS2a + figS2b
figS2

#ggsave("../Figures/FigureS2.png", figS2, device = "png",
#       height = 100, width = 180, units = "mm", dpi = 600)

#FIGURE S3-------------------------------------------------------
#create nice labels
main_mod_labs <- c("b_scale_patch_size" = "Total patch size",
                   "b_scale_log_cukes_final" = "log(Sea cucumber\ndensity)",
                   "b_scale_log_cukes_final:scale_log_grunt_pee_m2" = "log(Sea cucumber\ndensity) * log(Fish\nnutrient input)",  
                   "b_scale_log_grunt_pee_m2" = "log(Fish nutrient\ninput)", 
                   "b_scale_mean_PC1"="Initial quantity\nindex (PC1)",
                   "b_scale_distance"="Distance from reef")
main_mod_labs_treatment <- c("b_scale_patch_size" = "Total patch size",
                   "b_scale_cont_treatment"="Change in sea\ncucumber density",
                   "b_scale_cont_treatment:scale_log_grunt_pee_m2" = "Change in sea\ncucumber density *\nlog(Fish nutrient\ninput)",                   
                   "b_scale_log_grunt_pee_m2" = "log(Fish nutrient\ninput)", 
                   "b_scale_mean_PC1"="Initial quantity\nindex (PC1)",
                   "b_scale_distance"="Distance from reef")
#custom colours
coef_cols <- c("#B7CDBB", "#537970")

#extract the coefficients just like before
figS3a <- mod_main_density %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  filter(`.variable` != "b_Intercept") %>% 
  mutate(`.variable` = factor(`.variable`, 
                              levels = c("b_scale_patch_size",
                                         "b_scale_log_cukes_final",
                                         "b_scale_log_cukes_final:scale_log_grunt_pee_m2",
                                         "b_scale_log_grunt_pee_m2",
                                         "b_scale_mean_PC1",
                                         "b_scale_distance"))) %>% 
  ggplot(aes(y = .variable, x = .value, fill = stat(x < 0))) +
  stat_halfeye(size = 0.1) +
  geom_vline(xintercept = 0, size = 0.3) +
  scale_y_discrete(labels = main_mod_labs) +
  scale_fill_manual(values = coef_cols) +
  theme_paper_small() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 1,
                                   size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.x = 
          element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = -0.5, vjust = 2.5, size = 12)) +
  labs(x = "", y = "",
       title = "A - Quantity, Density")
figS3a
figS3b <- mod_main_treatment %>% 
  gather_draws(`b_.*`, regex = TRUE) %>% 
  filter(`.variable` != "b_Intercept") %>% 
  mutate(`.variable` = factor(`.variable`, 
                              levels = c("b_scale_patch_size",
                                         "b_scale_cont_treatment",
                                         "b_scale_cont_treatment:scale_log_grunt_pee_m2",
                                         "b_scale_log_grunt_pee_m2",
                                         "b_scale_mean_PC1",
                                         "b_scale_distance"))) %>% 
  ggplot(aes(y = .variable, x = .value, fill = stat(x < 0))) +
  stat_halfeye(size = 0.1) +
  geom_vline(xintercept = 0, size = 0.3) +
  scale_y_discrete(labels = main_mod_labs_treatment) +
  scale_fill_manual(values = coef_cols) +
  theme_paper_small() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 1,
                                   size = 7),
        axis.text.x = element_text(size = 7),
        #create same margin as panel A
        axis.title.x = 
          element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = -0.5, vjust = 2.5, size = 12)) +
  labs(x = "", y = "",
       title = "B - Quantity, Treatment")
figS3b

#combine
figS3 <- figS3a + figS3b
figS3
#ggsave("../Figures/FigureS3.png", figS3, device = "png",
#       height = 90, width = 180, units = "mm", dpi = 600)
#FIGURE 2A - CUKE EFFECT ON GROWTH---------------------------------------------
#create a dataset to predict on
#we want to set every predictor to its mean value, except for cukes (we want
#a sequence across its entire range) and fish (we want to predict across three 
#different levels to see the interaction)
newdata <- tibble(scale_log_cuke_density_final = rep(seq(from = -1.69, to = 1.71, 
                                                         by = 0.01), times = 3),
                  #we'll use the mean and the mean +/- 1SD to get a good range
                  scale_log_grunt_pee_m2 = c(rep(-1, times = 341), 
                                             rep(0, times = 341),
                                             rep(1, times = 341)),
                  scale_distance = rep(0),
                  scale_below_per_shoot = rep(0),
                  scale_patch_size = rep(0))

#now we'll add fitted draws, which means we calculate the predicted value for
#each observation in our newdata df, based on the parameter estimates for every
#iteration of the model
predict_growth <- newdata %>% 
  add_fitted_draws(mod_growth, re_formula = NA, seed = 123) %>% 
  #convert the cuke densities back into their unlogged and unscaled form
  mutate(final_cuke_density = exp((scale_log_cuke_density_final*1.52)-2.34) - 
           0.5*0.0150220,
         scale_log_grunt_pee_m2_factor = factor(scale_log_grunt_pee_m2),
         #increase contrast a little bit because the blues are a little
         #too similar - this won't impact the points
         scale_log_grunt_pee_m2 = 1.5*scale_log_grunt_pee_m2)

#now plot!
fig2a <- ggplot(predict_growth, aes(x = final_cuke_density, 
                           y = blade_height,
                           fill = rev(scale_log_grunt_pee_m2),
                           colour = rev(scale_log_grunt_pee_m2))) +
  #raw points, jittered slightly to help with the overlap
  geom_jitter(data = growths, size = 0.5, width = 0.02) +
  #stat_lineribbon plots the 95% credible interval of the predictions and the
  #mean estimate
  stat_lineribbon(aes(y = .value, 
                      group = scale_log_grunt_pee_m2_factor,
                      fill = after_scale(alpha(colour, 0.4))), 
                  .width = 0.95,
                  size = 0.8) +
  labs(x = expression(paste("Sea cucumber density (m"^-2,")")),
       y = "New growth per shoot (cm)",
       colour = "Standardized logged fish pee",
       title = "A") +
  theme_paper_small() + 
  theme(axis.title.x = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0),
                                    size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12),
        legend.position = "none"
  )
fig2a

#FIGURE 2B - GROWTH COEFFICIENTS-----------------------------------------------
#create nice labels again for the coefficient plots
growth_mod_labs <- c("b_scale_log_grunt_pee_m2"="log(Fish nutrient\ninput)", 
                     "b_scale_log_cuke_density_final"="log(Sea cucumber\ndensity)",
                     "b_scale_log_grunt_pee_m2:scale_log_cuke_density_final"= "log(Sea cucumber\ndensity) * log(Fish\nnutrient input)",
                     "b_scale_patch_size"="Total patch size",
                     "b_scale_below_per_shoot"="Belowground biomass\nper shoot",
                     "b_scale_distance"="Distance from patch")
coef_cols <- c("#B7CDBB", "#537970")

#and the same old posterior plot again!
fig2b <- mod_growth %>% 
  gather_draws(b_scale_distance,
               b_scale_log_cuke_density_final,
               b_scale_log_grunt_pee_m2,
               `b_scale_log_grunt_pee_m2:scale_log_cuke_density_final`,
               b_scale_patch_size,
               b_scale_below_per_shoot) %>% 
  mutate(`.variable` = factor(`.variable`, 
                              levels = c("b_scale_log_grunt_pee_m2",
                                         "b_scale_log_cuke_density_final",
                                         "b_scale_log_grunt_pee_m2:scale_log_cuke_density_final",
                                         "b_scale_patch_size",
                                         "b_scale_below_per_shoot",
                                         "b_scale_distance"))) %>% 
  ggplot(aes(y = .variable, x = .value, fill = stat(x < 0))) +
  stat_halfeye(size = 0.1) +
  geom_vline(xintercept = 0, size = 0.3) +
  scale_y_discrete(labels = growth_mod_labs) +
  scale_fill_manual(values = coef_cols) +
  theme_paper_small() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 1,
                                   size = 7),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.x = element_text(size = 7),
        #create same margin as panel A
        axis.title.x = 
          element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        plot.title = element_text(hjust = -0.5, vjust = 2.5, size = 12)) +
  labs(x = "", y = "",
       title = "B")
fig2b
#FIGURE 2 COMBINED---------------------------------------------------------
#combine into 1
fig2 <- fig2a + fig2b + plot_layout(design = "11122
                                    11122")
fig2
#ggsave("../Figures/Figure2.png", fig2, device = "png",
#       height = 100, width = 180, units = "mm", dpi = 600)
