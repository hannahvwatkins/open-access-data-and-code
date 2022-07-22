#if you open the RProject found in the Code folder before opening this, it will
#automatically set the working directory of the entire project and all of the
#filepaths will work on your own computer. If you open this just as a script, 
#you will need to set the working directory to wherever on your computer the 
#patch_movement/Code folder is located

#load necessary packages
library(tidyverse) #data wrangling/plotting
library(janitor) #data cleaning
library(performance) #model checks
library(see) #plotting model checks
library(MASS) #negative binomial modelling
library(brms) #Bayesian modelling to check for issues with MASS optimizer
library(tidybayes) #for plotting brms output
library(pscl) #zero inflation modelling
library(visreg) #for quick checks of the model outputs in response space
library(lme4) #for checking random effect of family
library(MuMIn) #for AICc comparison
library(ggrepel) #for figures 2
library(fishualize) #for colours on figures
library(ciTools) #for generating prediction intervals
library(boot) #for calculating iverse logit
library(patchwork) #for comparing plots side by side
source("ggplot_paper_theme.R") #ggplot theme for figures

#load data
df <- read_csv("../Data/Rate.csv") %>% 
  clean_names() %>% 
  mutate(x_rate = as.integer(rate_2),
         log_dist = log(distance_2)) %>% 
  filter(!is.na(x_rate))
df_l <- read_csv("../Data/rlength.csv")%>% 
  clean_names() %>% 
  filter(!is.na(species)) %>% 
  mutate(log_distance = log(distance),
         log_length = log(length))
spp <- read_csv("../Data/sppspecific.csv") %>% 
  dplyr::select(species:distnear)
sizes <- read_csv("../Data/species_sizes.csv")

#Fishbase and CEI record comparison---------------------------------------------
mod_size <- lm(Fishbase_length ~ CEI_length, data = sizes)
summary(mod_size)

#MODEL 1: determine effect of distance on frequency of reef-associated fish----
#with only one predictor of interest in this question, the only decisions
#we need to make when building this model are about the distribution we use
#since the response is a count, we'll try poisson first and then test for 
#overdispersion and zero inflation
dist_mod <- glm(x_rate ~ log_dist, data = df, family = "poisson")
summary(dist_mod)
#the ratio of null to residual deviance looks a little high, so we likely
#are dealing with overdispersion
check_overdispersion(dist_mod)
#are we EVER

#let's try a negative binomial instead
dist_mod_nb <- glm.nb(x_rate ~ log_dist, data = df)
summary(dist_mod_nb)
#MASS is struggling with this model, but is able to converge if we set a low 
#initial theta
dist_mod_nb2 <- glm.nb(x_rate ~ log_dist, data = df, init.theta = 1)
summary(dist_mod_nb2)
#theta is the shape parameter of the nb distribution, and glm.nb from the MASS
#packages appears to be struggling to estimate it because the overdispersion is 
#SO high. I believe MASS's estimator struggles when the shape of the nb 
#distribution is extremely different from Poisson, and theta gets pushed towards
#infinity. We can do a quick check to see if it's a MASS specific problem or a
#model specific problem by comparing it to the same model run in brms:
#note: if you don't have a C++ compiler, this model won't run in R for you, so
#skip the next line and uncomment the readRDS() line to load the model output
#from the Data folder
dist_mod_nb_brms <- brm(x_rate ~ log_dist, data = df, 
                        family = negbinomial)
#saveRDS(dist_mod_nb_brms, "../Data/brms_model.rds")
#dist_mod_nb_brms <- readRDS("../Data/brms_model.rds")
summary(dist_mod_nb_brms)
#the brms model converges nicely and does not have issues estimating theta 
#(i.e., the shape parameter), indicating it is more likely just an issue with
#the way the MASS package is able to handle extreme cases of overdispersion with 
#small sample sizes. The brms model is quite similar to the MASS model when
#init.theta is set <30, indicating this is an appropriate starting value

#let's see if there's also an issue of zero inflation
check_zeroinflation(dist_mod_nb2)
#there seems to be slight zero inflation (though its hard to tell because the
#sample size is really small)
#lets check if a zero inflated nb works better using the pscl package
dist_mod_zi <- zeroinfl(x_rate ~ log_dist, data = df, dist = "negbin")
summary(dist_mod_zi)
#ok that ran smoothly, so lets test whether zi is a better fit with the pscl 
#package
vuong(dist_mod_zi, dist_mod_nb2)
#it is not a better fit, meaning the glm.nb wasn't having issues because of 
#zero inflation

#in an ideal world, everything would be Bayesian all the time because the MCMC
#sampler can handle all the weird cases we bump into in ecology, but since the
#glm.nb() model (with an init.theta set <30) is both qualitatively and 
#quantitatively similar to the more accurate brms model, it appears to be a 
#reasonable approximation of the true trend. It is also the best frequentist 
#approximation available in R at this time!

#finally, we'll run a few model checks to make sure the model is appropriate
#for our data
check_model(dist_mod_nb2)
#none of the observations seem too worrisome other than observation 18, which 
#has an unusually high rate of fish for its distance, but its Cook's distance is 
#~0.5 so it's not TOO influential. It does seem to be driving the slight 
#abnormality in the qqplot, but overall there seems to be very little structure 
#in the residuals. Given the small sample size of the dataset, the model 
#diagnostics look good.


#MODEL 1 adjusted for varying angles of observation--------------------
#we'll first estimate how much of the circle (in degrees) a gopro should be
#able to see for each distance from the patch, assuming all fish within 10m
#of the camera can be recorded. Visibility during the surveys was 25-30m but
#it is very unlikely that we could see a fish that far from the camera. I picked
#10 m based on other studies using cameras underwater, and since we're only 
#interested in the relative differences between cameras at different distances,
#I think this approximation is ok
#this calculation is done based on trigonometry used to calculate the angle away 
#from the straight line between nearest neighbours that the camera would be
#able to capture. See appendix for diagram
df <- df %>% 
  mutate(angle_viz = ((atan(10/distance_2) * 180) / pi) * 2,
         x_rate_per_deg = x_rate/angle_viz)

#and now use the angle as an offset in the model or as the denominator for the
#response. These should give identical results but with differences in the 
#units of the response
mod_offset <- glm.nb(x_rate ~ log_dist + offset(log(angle_viz)), 
                  data = df)
mod_offset_brm <- brm(x_rate ~ log_dist + offset(log(angle_viz)), 
                     data = df, family = "negbinomial")
#saveRDS(mod_offset_brm, "../Data/offset_brms_model.rds")
summary(mod_offset_brm)
#MASS does NOT like the offset approach here (probably because of the massive
#range of offset values) and won't run, and brms also hates it and won't 
#converge (probably also because of the small sample size, and the correlation
#between the offset and the fixed effect)
#let's go with using it in the denominator of our response instead
#since turning our response into a rate means it is no longer an integer,
#we'll change the distribution to a hurdle Gamma (since there are zeros)
#we want log_dist to influence both the prob of getting a zero and the value
#of non-zero responses
mod_response <- brm(bf(x_rate/angle_viz ~ scale(log_dist), 
                    hu ~ scale(log_dist)), 
                    data = df, 
                    family = hurdle_gamma(),
                    #set seed for reproducibility
                    seed = 123)

mod_response <- brm(x_rate/angle_viz + 0.01 ~ log_dist, 
                    data = df, 
                    family = Gamma(),
                    #set seed for reproducibility
                    seed = 123)

#saveRDS(mod_response, "../Data/gamma_hurdle_model.rds")
mod_response <- readRDS("../Data/gamma_hurdle_model.rds")
summary(mod_response)
#that is able to run smoothly, and it looks like regardless of the assumptions 
#we make about the directionality of movement between patches, the trend stays 
#the same

#finally, we'll do some basic model checks
pp_check(mod_response) + xlim(0, 10)
#the density plot is a little iffy, but the major fluctuations seem most likely
#to be the result of our really small sample size
fitted_gamma <- fitted(mod_response)[,1]
resid_gamma <- resid(mod_response)[,1]

#we can look at the residuals too
gamma_check <- tibble(fitted = fitted_gamma,
                      resid = resid_gamma) %>% 
  cbind(df)
ggplot(gamma_check, aes(fitted, resid)) +
  geom_point()
ggplot(gamma_check, aes(x_rate_per_deg, resid)) +
  geom_point()
ggplot(gamma_check, aes(x_rate_per_deg, fitted)) +
  geom_point()
#again, not as perfect as we'd hope, but there doesn't seem to be much structure
#(remember that since this is a Gamma model, values below zero are not possible
#so it's normal to see a diagonal edge to the residuals representing their upper
#bound)
#given the small sample size, it is unlikely that we can find a better 
#alternative (I already tried reparameterizing the hurdle component of the model
#with no improvement)

#PLOT 1---------------------------------------
#we'll use predict to extract our model predictions and CI to plot the exact
#values from our model rather than ggplot's default glm.nb, since we set a
#different init than the default
#create a new dataframe to predict on
#note that this has to be done with the log_dist term - doing it with log(dist)
#in the model means the model will work fine but it'll mess up the predictions
newdata <- tibble(log_dist = seq(from = 0.1, to = 7.7, by = 0.02))
predict_nb <- as_tibble(predict(dist_mod_nb2, newdata = newdata, 
        #keep in link space to calculate CIs first, then backtransform                        
                                 type = "link", se.fit = TRUE)) %>% 
  cbind(newdata) %>% 
  dplyr::select(-residual.scale) %>% 
  #use the inverse of natural log to backtransform
  mutate(lowerCI = exp(fit - (1.96 * se.fit)),
         upperCI = exp(fit + (1.96 * se.fit)),
         x_rate_adj = exp(fit),
         #use inverse of log to backtransform for viz
         distance_2 = exp(log_dist),
         exp_fit = exp(fit))

#to estimate our effect size, we can see what the proportional change in the
#estimated rate per degree is at each order of magnitude increase\
#we'll check between a few levels just to make sure this rate is constant
e1 <- predict_nb %>% 
  filter(log_dist == 1) %>% 
  dplyr::select(exp_fit) %>% 
  pull()
e2 <- predict_nb %>% 
  filter(log_dist == 2) %>% 
  dplyr::select(exp_fit) %>% 
  pull()
e3 <- predict_nb %>% 
  filter(log_dist == 3) %>% 
  dplyr::select(exp_fit) %>% 
  pull()
e4 <- predict_nb %>% 
  filter(log_dist == 4) %>% 
  dplyr::select(exp_fit) %>% 
pull()

rate1 <- (e1 - e0)/e0
rate2 <- (e2 - e1)/e1
rate3 <- (e3 - e2)/e2
rate4 <- (e4 - e3)/e3
#excellent

#can also just use this formula lol
effect_size <- exp(coef(dist_mod_nb2)[2]) - 1

#for plotting purposes, we're going to need to add a small value to the zeros
#so they still show up on the figure. If we add it to the df rather than
#in ggplot, the model fit and CI will be completely unaffected
df <- df %>% 
  mutate(x_rate_ch = as.character(x_rate),
         x_rate_adj = as.numeric(case_when(x_rate_ch == '0' ~ '0.01',
                                TRUE ~ x_rate_ch)))

#here is the model in link (log) space which I think is the best
#way to present it. Since the line and CI are coming from our model
#output and not ggplot's defaults, adding the 0.01 to keep the zeros 
#in the figure doesn't change that line
negbinplot <- ggplot(df, aes(distance_2, (x_rate_adj))) +
  geom_point() +
  labs(x = "Distance to nearest patch reef (m)", 
       "Total number of individuals observed/hour") +
  theme_paper_large() + 
  theme(legend.text = element_text(size=15)) +
  geom_ribbon(data = predict_nb, aes(ymin = lowerCI, ymax = upperCI), 
              alpha = 0.1) +
  geom_smooth(data = predict_nb, se = FALSE, col = "black")  +
  scale_y_continuous(tran = "log10", 
                     labels = c("0.1", "1", "10", "100", "1000"),
                     breaks = c(0.1, 1, 10, 100, 1000)) +
  #note that since we inversed the log in our prediction dataframe, we can
  #transform this scale to whatever we want, since all of the data we are using
  #for this figure is neither log nor log10 transformed 
  scale_x_continuous(tran = "log10") +
  #note that we're going to call this the rate and not rate + 0.01
  #because the line and CI are truly the rate predictions, and it's only the few
  #points at 0 that needed the +0.01 to show up on the figure
  labs(y = "Number of individuals observed per hour",
       x = "Distance to nearest patch reef (m)")
negbinplot
#ggsave("../Figures/distance_vs_rate_log.png",
#       width = 180, height = 110, units = "mm", dpi = 600)

#here it is in real space with CIs
ciplot <- ggplot(df, aes(log_dist, x_rate_adj)) +
  geom_point() +
  labs(x = "Distance to nearest patch reef (m)", 
       "Total number of individuals observed/hour") +
  theme_paper() + 
  theme(legend.text = element_text(size=15)) +
  geom_ribbon(data = predict_nb, aes(ymin = lowerCI, ymax = upperCI), 
              alpha = 0.1) +
  geom_smooth(data = predict_nb, se = FALSE, col = "black") 
ciplot
#ggsave("../Figures/distance_vs_rate_real.png",
#       width = 180, height = 110, units = "mm", dpi = 600)
#PLOT 1 adjusted--------------
newdata2 <- tibble(log_dist = seq(from = 0, to = 7.7, by = 0.02), angle_viz=50)
predict_gamma <- newdata2 %>% 
  add_fitted_draws(mod_response) %>% 
  mutate(log_fit = log(`.value`)) %>% 
  group_by(`.row`) %>% 
  summarize(log_mean_fit = mean(log_fit),
            log_lowerCI = quantile(log_fit, probs = 0.025),
            log_upperCI = quantile(log_fit, probs = 0.975),
            log_dist = mean(log_dist)) %>% 
  ungroup() %>% 
  mutate(distance_2 = exp(log_dist),
         mean_fit = exp(log_mean_fit),
         lowerCI = exp(log_lowerCI),
         upperCI = exp(log_upperCI))
#to estimate our effect size, we can see what the proportional change in the
#estimated rate per degree is at each order of magnitude increase
#to make it directly comparable to our original negative binomial model, we have
#to use base e. We can find these values in our predict_gamma dataframe
e0 <- predict_gamma %>% 
  filter(log_dist == 0) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e1 <- predict_gamma %>% 
  filter(log_dist == 1) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e2 <- predict_gamma %>% 
  filter(log_dist == 2) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e3 <- predict_gamma %>% 
  filter(log_dist == 3) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e4 <- predict_gamma %>% 
  filter(log_dist == 4) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e5 <- predict_gamma %>% 
  filter(log_dist == 5) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e6 <- predict_gamma %>% 
  filter(log_dist == 6) %>% 
  dplyr::select(mean_fit) %>% 
  pull()
e7 <- predict_gamma %>% 
  filter(log_dist == 7) %>% 
  dplyr::select(mean_fit) %>% 
  pull()

rate1 <- (e1 - e0)/e0
rate2 <- (e2 - e1)/e1
rate3 <- (e3 - e2)/e2
rate4 <- (e4 - e3)/e3
rate5 <- (e5 - e4)/e4
rate6 <- (e6 - e5)/e5
rate7 <- (e7 - e6)/e6

mean(c(rate1, rate2, rate3, rate4, rate5, rate6, rate7))

#the relationship isn't perfectly exponential anymore because the model has both
#a log link and logit link, but on average, the rate declines by about 35% for 
#ever order of magnitude increase in the distance

#for plotting purposes, we're going to need to add a small value to the zeros
#so they still show up on the figure. If we add it to the df rather than
#in ggplot, the model fit and CI will be completely unaffected
df <- df %>% 
  mutate(x_rate_ch = as.character(x_rate),
         x_rate_adj = as.numeric(case_when(x_rate_ch == '0' ~ '0.01',
                                           TRUE ~ x_rate_ch)),
         x_rate_deg_true = as.character(x_rate_adj/angle_viz),
         x_rate_deg = as.numeric(case_when(x_rate_adj == "0.01" ~ "0.0001",
                                TRUE ~ x_rate_deg_true)))

#here is the model in link (log) space which I think is the best
#way to present it. Since the line and CI are coming from our model
#output and not ggplot's defaults, adding the 0.01 to keep the zeros 
#in the figure doesn't change that line
gammaplot <- ggplot(df, aes(distance_2, (x_rate_deg))) +
  geom_point() +
  labs(x = "Distance to nearest patch reef (m)", 
       y = "Total number of individuals observed/degree/hour") +
  theme_paper_large() + 
  theme(legend.text = element_text(size=15)) +
  geom_ribbon(data = predict_gamma, aes(y = mean_fit, 
                                      ymin = lowerCI, ymax = upperCI), 
              alpha = 0.1) +
  geom_smooth(data = predict_gamma, aes(y = mean_fit), 
              se = FALSE, col = "black") +
  #stat_lineribbon(data = predict_gamma, aes(x = exp(log_dist),
  #                                          y = `.value`)) +
  scale_y_continuous(tran = "log10") +
  scale_x_continuous(tran = "log10")
  #note that we're going to call this the rate and not rate + 0.01
  #because the line and CI are truly the rate predictions
gammaplot
#ggsave("../Figures/distance_vs_rate_adjusted.png",
#       width = 180, height = 110, units = "mm", dpi = 600)

#since the CI for log_dist overlaps zero, but the relationship is still
#quite strong, we can also just examine the probability of observing zero fish
#on a trap
prob_zero <- newdata2 %>% 
  mutate(logistic = (1/(1 + exp(-(-5.51 + 1.2*log_dist)))),
         distance = exp(log_dist))

ggplot(prob_zero, aes(distance, logistic)) +
  geom_line()
#this is a super strong relationship! So it makes sense that most of the
#variation isn't explained by the rate part (i.e., the CI on the parameter 
#estimate for distance overlaps 0), because it's already being explained by
#the hurdle

#MODEL 2: determine effect of body size on max distance----------
#length_mod_brms <- brm(distance ~ length, data = df_l, family = Gamma(link = "log"))
#summary(length_mod)

#the continuous, 0 bounded data we have can likely be represented by either a 
#Gamma distribution or a lognormal distribution. We can fit both models and 
#exmaine the residuals to determine which is a better fit
length_mod <- lm(log_distance ~ length, data = df_l)
summary(length_mod)
visreg(length_mod)

length_mod_gamma <- glm(distance ~ length, data = df_l, 
                        family = Gamma(link = "log"))
summary(length_mod_gamma)
visreg(length_mod_gamma, scale = "linear")

#both models yield the same qualitative output, although the lognormal model
#is only marginally significant and has a slightly lower effect size. However, 
#by examining the visreg plots, the lognormal model seems as though it may
#be a better fit to the data, as in the Gamma model the 3 extremely 
#high-distance points appear to have a lot of leverage, leading most of the 
#other observations to fall below the line
#we can do a more formal observation of the residuals as well
plot(length_mod)
plot(length_mod_gamma)

#in both cases, the "outlier" points don't look great, but appear to be causing
#more problems in the Gamma model. In the resid vs fitted plot, the points are
#much more evenly spread around 0 than in the Gamma model. Since both models
#yield the same qualitative results and the lognormal model appears to be a
#better fit for the data, we will carry on with this model

#let's compare that simple lm to a mixed effects model accounting
#for the potential effect of family, just in case that is driving the trends
#we see
length_mod_fam <- glmer(distance ~ length + (1|family), data = df_l,
                        family = Gamma(link = "log"))
summary(length_mod_fam)

#can compare the goodness-of-fit using anova
anova(length_mod_fam, length_mod)
#huge difference here, but since our sample size is small we can also compare 
#the models using AIC corrected for small sample size
AICc(length_mod, length_mod_fam)
#adding in a random effect of family really does not improve the model, as it is
#~119 AIC greater than the basic lm. It also doesn't change the effect of length

#we also have those three points that look like they may be outliers (although
#they do not appear to have too much leverage based on their Cook's distance)
#we can quickly make sure they aren't changing our results by dropping them and
#seeing if the results change much
outliers <- filter(df_l, distance < 80)
length_mod_out <- glm(distance ~ length, data = outliers, 
                      family = Gamma(link = "log"))
summary(length_mod_out)
#the relationship is a little less steep, but there is still a clear, 
#significant effect of body size on distance traveled, with or without those 
#extreme values. Since we have no biological reason to drop these values, we 
#should keep them in the analysis, but I think it's also valuable to show the
#model without these points since it changes the effect size

#let's just do some model checks on our selected model and the outlier one
plot(length_mod_gamma)
plot(length_mod_out)
#the residuals look good apart from the three very high points affecting the
#tail of the qqplot and the density plot. Without dropping those points 
#(which we've already confirmed does not substantially impact the analysis), 
#there is no way to deal with these, and they seem like reasonable violations

#to explain the effect size in real space rather than log space, we just need
#to inverse log the coefficient
effect_size_length <- exp(coef(length_mod)[2]) - 1
effect_size_out <- exp(coef(length_mod_out)[2]) - 1
#this means for every unit increase in length, we expect a 5% increase in 
#distance travelled if we keep the outliers in, and a 3% increase if we drop 
#them

#PLOT 2-------------------------------
#use predictions from Gamma model
#newdata3 <- tibble(length = seq(from = 5, to = 30, by = 0.02))
#predict_dist <- as_tibble(predict(length_mod_gamma, newdata = newdata3, 
#                                #keep in link space to calculate CIs first, 
#                                #then backtransform                        
#                                type = "link", se.fit = TRUE)) %>% 
#  cbind(newdata3) %>% 
#  dplyr::select(-residual.scale) %>% 
#  #use the inverse of natural log to backtransform
#  mutate(lowerCI = exp(fit - (1.96 * se.fit)),
#         upperCI = exp(fit + (1.96 * se.fit)),
#         distance = exp(fit))
#
#dist_plot <- ggplot(df_l, aes(length, distance)) +
#  geom_point() +
#  theme_paper() + 
#  geom_ribbon(data = predict_dist, aes(ymin = lowerCI, ymax = upperCI), 
#              alpha = 0.1) +
#  geom_smooth(data = predict_dist, se = FALSE, col = "black")  +
#  scale_y_continuous(tran = "log10") + 
#  #note that since we inversed the log in our prediction dataframe, we can
#  #transform this scale to whatever we want, since all of the data we are using
#  #for this figure is neither log nor log10 transformed 
#  labs(y = "Maximum observed distance travelled (m)",
#       x = "Mean species body length (cm)")
#dist_plot
#
##now do the same thing for the outlier model
#predict_out <- as_tibble(predict(length_mod_out, newdata = newdata3, 
#                                  #keep in link space to calculate CIs first, 
#                                  #then backtransform                        
#                                  type = "link", se.fit = TRUE)) %>% 
#  cbind(newdata3) %>% 
#  dplyr::select(-residual.scale) %>% 
#  #use the inverse of natural log to backtransform
#  mutate(lowerCI = exp(fit - (1.96 * se.fit)),
#         upperCI = exp(fit + (1.96 * se.fit)),
#         distance = exp(fit))
##and plot
#out_plot <- ggplot(outliers, aes(length, distance)) +
#  geom_point() +
#  theme_paper() + 
#  geom_ribbon(data = predict_out, aes(ymin = lowerCI, ymax = upperCI), 
#              alpha = 0.1) +
#  geom_smooth(data = predict_out, se = FALSE, col = "black")  +
#  scale_y_continuous(tran = "log10") + 
#  #note that since we inversed the log in our prediction dataframe, we can
#  #transform this scale to whatever we want, since all of the data we are using
#  #for this figure is neither log nor log10 transformed 
#  labs(y = "Maximum observed distance travelled (m)",
#       x = "Mean species body length (cm)")
#out_plot
#
#dist_plot + out_plot
#
##we could also try putting both predictions on the same plot
#combined_plot <- ggplot(df_l, aes(length, distance)) +
#  geom_point() +
#  theme_paper() + 
#  geom_ribbon(data = predict_dist, aes(ymin = lowerCI, ymax = upperCI), 
#              alpha = 0.1) +
#  geom_smooth(data = predict_dist, se = FALSE, col = "black")  +
#  geom_ribbon(data = predict_out, aes(ymin = lowerCI, ymax = upperCI), 
#              alpha = 0.1) +
#  geom_smooth(data = predict_out, se = FALSE, col = "black")  +
#  scale_y_continuous(tran = "log10") + 
#  #note that since we inversed the log in our prediction dataframe, we can
#  #transform this scale to whatever we want, since all of the data we are using
#  #for this figure is neither log nor log10 transformed 
#  labs(y = "Maximum observed distance travelled (m)",
#       x = "Mean species body length (cm)")
#
#combined_plot
#
#use predictions from linear model
#option 1: colour by family, with legend
ggplot(df_l, aes(length, distance)) +
  stat_smooth(method = "lm", col = "black") +
  geom_point(aes(colour = family)) +
  labs(x = "Mean fish length (cm)", 
       y = "Maximum observed distance travelled (m)") +
  theme_paper_large() +
  scale_y_continuous(tran = "log10") +
  scale_colour_fish_d(option = "Hypsypops_rubicundus")

#ggsave("../Figures/length_vs_distance1.png",
#       width = 180, height = 110, units = "mm", dpi = 600)

#option 2: colour by species, no legend, with labels
ggplot(df_l, aes(length, distance)) +
  stat_smooth(method = "lm", col = "black") +
  geom_point(aes(colour = species)) +
  labs(x = "Mean fish length (cm)", 
       y = "Maximum observed distance travelled (m)") +
  theme_paper_large() +
  scale_y_continuous(tran = "log10") +
  scale_colour_fish_d(option = "Hypsypops_rubicundus") +
  theme(legend.position = "none") +
  geom_text_repel(aes(label=species), size = 2, 
                  box.padding = unit(0.75, "lines"),
                  max.overlaps = 100, max.time = 10,
                  force_pull = 10)

#ggsave("../Figures/length_vs_distance2.png",
#       width = 180, height = 110, units = "mm", dpi = 600)

#option 3: no colour, no legend, with labels
ggplot(df_l, aes(length, distance)) +
  stat_smooth(method = "lm", col = "black") +
  geom_point() +
  labs(x = "Mean fish length (cm)", 
       y = "Maximum observed distance travelled (m)") +
  theme_paper_large() +
  scale_y_continuous(tran = "log10") +
  theme(legend.position = "none") +
  geom_text_repel(aes(label=species), size = 2, 
                  box.padding = unit(0.75, "lines"),
                  max.overlaps = 100, max.time = 10,
                  force_pull = 10)

#ggsave("../Figures/length_vs_distance3.png",
#       width = 180, height = 110, units = "mm", dpi = 600)



