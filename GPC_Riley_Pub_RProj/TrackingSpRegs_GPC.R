##############################################################################
################ Tracking spatial regimes as an early warning ################
################### for a species of conservation concern ####################
############################## Roberts et al. ################################
##############################################################################

#=============================================================================
## Preparations

# # Clear environment?
# rm(list=ls())

# List of packages necessary to run this script:
needed.packs <- list("data.table", "stringr", "reshape2", "ggplot2", 
                     "perm", "plyr", "ggfortify", "car", "here",
                     "vegan", "rgdal", "rgeos", "sp", "raster", "boot",
                     "lme4", "devtools", "parallel", "parallelsugar",
                     "AICcmodavg", "mgcv", "MASS", "boot", "scales",
                     "cowplot", "ncf", "cowplot", "nlme")

# Install/Load packages
package.check <- lapply(needed.packs, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# install_github('nathanvan/parallelsugar')

# How many cores
numCores <- 30

# Projection
projUTM <- "+proj=utm +zone=14 +datum=WGS84 +units=m"

### Load data:
dat <- fread(here("TrackingSpRegs_GPC_Data.csv"))

# Years of interest
YRS <- sort(unique(dat$Year))

#=============================================================================
## Model selection

# Fit the models
mods <- lapply(c(paste0("present ~ spCov_", 
                        c(3, 13, 33, 45, 57, 83, 107, 119, 133),
                        "_scale + (1|Year)"),
                 "present ~ (1|Year)"),
               as.formula)

# Name the models for the model selection table
names(mods) <- c(str_extract(sapply(mods, function(x) as.character(x)[3]), 
                             "spCov_[0-9]{1,3}_scale")[1:(length(mods) - 1)],
                 "Null")

# Fit the models
fits <- lapply(mods, 
               function(x) glmer(x, family = binomial, data = dat))

# Make a model selection table and save it
modSel_tab <- aictab(fits)
write.csv(as.data.frame(modSel_tab),
          here("GPC_Results/ModelSelection_Table.csv"),
          row.names = FALSE)

### Get the model coefficients for all candidate models.
# This is for a Supplementary Table S1.
tableS1 <- rbindlist(
  Map(x = fits,
      y = names(fits),
      function(x, y) {
        temp <- cbind(fixef(x))
        df <- data.frame(Model = y,
                         Parameters = row.names(temp),
                         Estimates = temp[,1])
        return(df)
      }))
tableS1[, "Parameters" := str_replace(str_remove(Parameters, "_[0-9]+_scale"),
                                      "spCov", "SpCov")]
tableS1[, "Estimates" := round(Estimates, 2)]
tableS1[ , "Model" := paste0("lek presence ~ ",
                             ifelse(str_detect(Model, "[0-9]+"),
                                    paste0("SpCov ",
                                           ((as.numeric(str_extract(Model, "[0-9]+")) * 30)^2 * 0.0001),
                                           " ha + (1|Year)"),
                                    "(1|Year)"))]
tableS1 <- dcast(tableS1,
                 Model ~ Parameters,
                 value.var = "Estimates")
tableS1 <- tableS1[c(7, 8, 9, 10, 3, 4, 6, 5, 2, 1),]
write.csv(tableS1,
          here("GPC_Results/RSFModel_Coefficients.csv"),
          row.names = FALSE)

#=============================================================================
## Plot model predictions

# Get range of spatial covariance 
spCov = range(na.omit(dat$spCov_45_scale))

# data.frame for 'new.data' argument in predict()
nd <- expand.grid(spCov_45_scale = seq(spCov[1], spCov[2], length = 100), 
                  Year = factor(sort(unique(dat$Year))))

# Helper functions for bootstrapping confidence limits
myFunc1 <- function(mm) {
  predict(mm, newdata = nd, re.form = NULL)
} # Random effects
myFunc2 <- function(mm) {
  predict(mm, newdata = nd, re.form = ~0)
} # Just fixed effects

### Make a cluster for bootstrapping confidence limits
cl <- makeCluster(30)

# The top model for passing to the cluster
fit_good <- fits$spCov_45_scale

# Export the 'lme4' package and some objects in the global environment to 
# the cluster
clusterEvalQ(cl, library("lme4"))
clusterExport(cl, list("spCov", "nd", "myFunc1", "myFunc2", "fit_good"))

# Bootstrapping confidence limits for fixed and random effects
bigBoot1 <- bootMer(x = fit_good, 
                    FUN = myFunc1, 
                    nsim = 10000, 
                    parallel = "snow", 
                    cl = cl)
bigBoot2 <- bootMer(x = fit_good, 
                    FUN = myFunc2, 
                    nsim = 10000, 
                    parallel = "snow", 
                    cl = cl)

# Close the cluster
stopCluster(cl)

# Get 95% CI for fixed and random effects
getCI <- function(boot, CI_low_name, CI_up_name) {
  CI_low <- apply(inv.logit(boot$t), 2, 
                  function(x) as.numeric(quantile(x, 
                                                  probs=.0275, 
                                                  na.rm=TRUE)))
  CI_up <- apply(inv.logit(boot$t), 2, 
                 function(x) as.numeric(quantile(x, 
                                                 probs=.975, 
                                                 na.rm=TRUE)))
  dt <- data.table(CI_low, CI_up)
  setnames(dt, names(dt), c(CI_low_name, CI_up_name))
  return(dt)
}

# Make predictions based on the random and fixed effects
nd$fit_random <- predict(fit_good, 
                         newdata = nd, 
                         type = "response")
nd$fit_fixed <- predict(fit_good, 
                        newdata = nd, 
                        type = "response", 
                        re.form = NA)

# Make a data.table with the predictions and bootstrapped confidence limits
nd <- data.table(nd,
                 # getCI(boot = bigBoot1, "CI_low_random", "CI_up_random"),
                 getCI(boot = bigBoot2, "CI_low_fixed", "CI_up_fixed"))

# Plot fixed effects
rsf_model_fixed <- ggplot(data = nd) +
  geom_ribbon(mapping = aes(x = spCov_45_scale, 
                            ymin = CI_low_fixed,
                            ymax = CI_up_fixed),
              fill = "grey70") +
  geom_line(mapping = aes(x = spCov_45_scale, 
                          y = fit_fixed,
                          group = Year),
            color = "green",
            size = 0.5) + 
  scale_x_continuous(limits = c(spCov[1], spCov[2])) +
  theme_minimal() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_blank()) +
  ylab("GPC Lek Occurence\nRelative Probability") + 
  xlab("Grass:Woody Boundary Strength (scaled)") +
  geom_rug(aes(x = spCov_45_scale),
           data = dat[present ==1],
           sides = "t", inherit.aes = FALSE)
rsf_model_fixed
ggsave(plot = rsf_model_fixed,
       filename = here("GPC_Results",
                       "RSFPrediction_Fixed_Figure 1.jpeg"),
       dpi = 600)

# Plot the Resource Selection Function model
rsf_model_random <- ggplot(data = nd) +
  geom_ribbon(mapping = aes(x = spCov_45_scale, 
                            ymin = CI_low_random,
                            ymax = CI_up_random,
                            group = Year),
              fill = "grey70") +
  geom_line(mapping = aes(x = spCov_45_scale, 
                          y = fit_random,
                          group = Year
                          ,
                          color = as.numeric(as.character(Year))
  ),
  # color = "green",
  size = 0.5) + 
  scale_x_continuous(limits = c(spCov[1], spCov[2])) +
  theme_minimal() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_blank()) +
  ylab("GPC Lek Occurence\nRelative Probability") + 
  xlab("Grass:Woody Boundary Strength (scaled)") +
  geom_rug(aes(x = spCov_45_scale),
           data = dat[present ==1],
           sides = "t", inherit.aes = FALSE)
# rsf_model_random
# ggsave(plot = rsf_model_random,
#        filename = here("GPC_Results",
#                        "RSFPrediction_Random.jpeg"),
#        dpi = 300)

# Calculate pseudo R-squared
r.squaredGLMM(fit_good)

#=============================================================================
## Make maps of spatial covariance and predicted lek occurrence

# Maneuver areas map of Ft Riley 
riley <- readOGR(dsn = here("Fort Riley Lek location maps"),
                 layer = "FtRiley_ManeuverAreas")
riley <- spTransform(riley,
                     CRS(proj4string(rast_list$rasters_107x107[[1]])))
riley_fort <- fortify(riley)
setnames(riley_fort, c("long", "lat"), c("x", "y"))

# Clip rasters to survey boundaries, scale rasters
rastALL_clip <- lapply(rast_list,
                       function(x){
                         lapply(x,
                                function(y) {
                                  scale(mask(crop(y, survey_bounds),
                                             survey_bounds),
                                        center = FALSE)
                                })
                       })
# fortify all rasters
rastALL_fort <- mclapply(rastALL_clip,
                         function(z) {
                           temp_ls <- lapply(z,
                                             function(y) {
                                               temp <- na.omit(as.data.frame(y, xy = TRUE))
                                               temp <- as.data.table(temp)
                                               return(temp)
                                             })
                           temp_join <- Reduce(join, temp_ls)
                           # return(temp_join)
                           temp_melt <- melt.data.table(temp_join,
                                                        id.vars = c("x", "y"))
                           return(temp_melt)
                         },
                         mc.cores = numCores)

# bind list
rastALL_bind <- rbindlist(rastALL_fort)

# Split columns
rastALL_bind[ , "Scale" := factor(str_extract(variable, "[0-9]{1,3}W"),
                                  levels = c("3W", "13W", "33W",
                                             "45W", "57W", "83W",
                                             "107W", "119W", "133W"))]
rastALL_bind[ , "Year" := str_extract(variable, "[0-9]{4}")]


#
# Clip rasters to survey boundaries
rasts_clip <- lapply(rast_list$rasters_45x45,
                     function(x) mask(crop(x, survey_bounds),
                                      survey_bounds))

# Create data.table of raster values by year
rast_dt <- rbindlist(Map(X = rasts_clip,
                         Y = names(rasts_clip),
                         function(X, Y) {
                           names(X) <- modSel_tab$Modnames[1]
                           X <- scale(X, center = FALSE)
                           temp <- na.omit(as.data.frame(X, xy = TRUE))
                           temp$Year <- Y
                           return(temp)
                         }))

# Which years to plot?
plot_years <- c(1994, 2006, 2017)

# Use the model to predict lek occurrence on the map
rast_pred <- predict(newdata = rast_dt,
                     object = fit_good,
                     type = "response",
                     re.form = NA)
rast_pred <- data.table(Probability = rast_pred,
                        rast_dt)

# Make a color ramp
probBreaks <- seq(min(rast_pred$Probability), 
                  max(rast_pred$Probability), 
                  0.01)
probRamp <- colorRampPalette(c("black", 
                               # "darkred",
                               # "red",
                               "orange", 
                               "palegoldenrod", 
                               "greenyellow",
                               "green",
                               "green3"))(length(probBreaks) - 1)


# Plot the model predctions on a map
pred_maps <- 
  ggplot(data = rast_pred[Year %in% plot_years],
         aes(x = x,
             y = y,
             fill = Probability)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               fortify(riley),
               fill= "grey80") +
  geom_tile() +
  north(x.min = min(rast_pred$x),
        x.max = max(rast_pred$x),
        y.min = min(rast_pred$y),
        y.max = max(rast_pred$y),
        symbol = 12,
        location = "bottomleft",
        scale = 0.2) +
  scale_fill_gradientn(colors = probRamp,
                       # values = probBreaks,
                       limits = c(min(rast_pred$Probability),
                                  max(rast_pred$Probability)),
                       name = "Relative\nProbability of\nGPC Lek\nOccurrence",
                       breaks = seq(0, 0.75, 0.25)) +
  facet_wrap(~ Year, ncol = 5) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 16),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Northing") +
  xlab("Easting")
# pred_maps

# ggsave(plot = pred_maps,
#        filename = here("GPC_Results",
#                        "LekOccurrenceProbability_45x45SpCov.jpeg"),
#        dpi = 300,
#        height = 4,
#        width = 11)


# ## Map spatial covariance
vegBreaks <- seq(-2.5, 
                 0.5,
                 0.01)
vegRamp <- colorRampPalette(c("darkred", 
                              "red", 
                              "orange", 
                              "palegoldenrod", 
                              "lightskyblue3"))(length(vegBreaks) - 1)
vegRamp <- c(rep(vegRamp[1], 
                 length(seq(min(rast_pred[Year %in% plot_years]$spCov_45_scale),
                            -2.51, 
                            0.01))),
             vegRamp,
             rep(vegRamp[length(vegRamp)], 
                 length(seq(0.5,
                            max(rast_pred[Year %in% plot_years]$spCov_45_scale),
                            0.01))))

spcov_maps <-
  ggplot() +
  geom_polygon(aes(x = x, y = y, group = group),
               riley_fort,
               fill= "grey80") +
  geom_tile(data = rast_pred[Year %in% plot_years],
            aes(x = x,
                y = y,
                fill = spCov_45_scale)) +
  scale_fill_gradientn(colors = vegRamp,
                       # values = rescale(vegBreaks),
                       limits = c(min(rast_pred[Year %in% plot_years]$spCov_45_scale),
                                  max(rast_pred[Year %in% plot_years]$spCov_45_scale)),
                       name = "Grass:Woody\nBoundary\nStrength") +
  facet_wrap(~ Year, ncol = 5) +
  coord_equal() +
  theme_bw() +
  ggsn::scalebar(x.min = min(rast_pred$x),
                 x.max = max(rast_pred$x),
                 y.min = min(rast_pred$y),
                 y.max = max(rast_pred$y),
                 dist = 2.5, st.size=3, height=0.02,
                 st.dist = .05,
                 dist_unit = "km",
                 transform = FALSE,
                 location="bottomleft" ) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 16),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  ylab("Northing") +
  xlab("Easting")
# spcov_maps

ggsave(plot = spcov_maps,
       filename = here("GPC_Results",
                       "VegTransitionsMap_45x45SpCov_Figure 3.jpeg"),
       dpi = 600,
       height = 4,
       width = 11)

# Combine prediction and spatial covariance plots
pred_cov <- plot_grid(spcov_maps + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                      pred_maps + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
                      ncol = 1,
                      rel_heights = c(0.5, 0.5),
                      align = "hv",
                      labels = c("A", "B"),
                      label_size = 18)
ggsave(plot = pred_cov,
       filename = here("GPC_Results/SpCov_LekProb_plot_Figure 2.jpeg"),
       dpi = 600,
       height = 8,
       width = 11)

#=============================================================================
## Temporal early warning

# Calculate median of spCov_57_scale over for each year. We won't
# use the mean because spCov is crazy skewed.
median45 <- rastALL_bind[Scale == "45W",
                         list(median = median(value, na.rm = T),
                              mean = mean(value, na.rm = T),
                              quant75 = quantile(value, 
                                                 probs = 0.75,
                                                 na.rm = T),
                              quant25 = quantile(value, 
                                                 probs = 0.25,
                                                 na.rm = T)), 
                         by = "Year"]
median45[ , "Year" := as.numeric(Year)]

# Fit mean spCov to GLS
sos_gls <- gls(median ~ Year,
               correlation = corAR1(),
               data = median45)

# Predict GLS
gls_fit <- as.data.frame(predictSE(sos_gls,
                                   newdata = data.frame(Year = 1994:2017),
                                   type = "response"))
gls_fit$Year <- 1994:2017

# Plot temporal early warning trend
safeOperate_plot <-
  ggplot() +
  geom_tile(data = expand.grid(spCov_45_scale = seq(-1.35, 0.1, 0.01),
                               Year = 1994:2017),
            mapping = aes(x = Year,
                          y = spCov_45_scale,
                          fill = spCov_45_scale)) +
  scale_fill_gradientn(values = rescale(vegBreaks),
                       colors = vegRamp) +
  geom_segment(aes(y = min(dat[present == 1]$spCov_45_scale, na.rm = T),
                   yend = min(dat[present == 1]$spCov_45_scale, na.rm = T),
                   x = 1993.5, xend = 2017.5),
               size = 2,
               color = "black") +
  geom_ribbon(data = gls_fit,
              aes(x = Year,
                  ymin = fit - se.fit * 1.96,
                  ymax = fit + se.fit * 1.96),
              fill = "grey70",
              alpha = 0.6) +
  geom_point(data = median45,
             aes(x = as.numeric(Year),
                 y = median),
             color = "grey40",
             size = 2) +
  geom_errorbar(data = median45,
                aes(x = as.numeric(Year),
                    ymin = quant25,
                    ymax = quant75),
                width = 0,
                color = "grey40") +
  geom_line(data = gls_fit,
            aes(x = Year,
                y = fit),
            color = "red") +
  scale_x_continuous(breaks = seq(1994, 2017, 3)) +
  scale_y_continuous(limits = c(-1.35, 0.1),
                     breaks = c(0, -0.25, -0.5, 
                                round(min(dat[present == 1]$spCov_45_scale, na.rm = T), 2), 
                                -1.25)) +
  annotate("text",
           x = 1999, y = 0,
           label = "Available Space",
           size = 4) +
  annotate("text",
           x = 2003, y = -1.27,
           label = "Unavailable Space",
           size = 4,
           color = "white") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 9),
        panel.grid = element_blank()) +
  ylab("Grass:Woody Boundary Strength (scaled)") +
  xlab("Year")
safeOperate_plot
ggsave(plot = safeOperate_plot,
       filename = here("GPC_Results/safeOperating_45W_Figure 3.jpeg"),
       dpi = 600,
       height = 4.25,
       width = 6)

# How much of the study area is within safe operating space?
nrow(rast_dt[Year == 1994 & spCov_45_scale >= -0.87])/nrow(rast_dt[Year == 1994])
nrow(rast_dt[Year == 2017 & spCov_45_scale >= -0.87])/nrow(rast_dt[Year == 2017])


#=============================================================================
## Supplementary Figures

# Plot the ALL YEARS spatial covariance on a map (supplementary figure)
spcov_Supplementary <-
  ggplot(data = na.omit(rast_pred),
         aes(x = x,
             y = y,
             fill = spCov_45_scale)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               fortify(riley),
               fill= "grey80") +
  geom_tile() +
  scale_fill_gradientn(colors = vegRamp,
                       # values = rescale(vegBreaks),
                       limits = c(min(rast_pred[Year %in% plot_years]$spCov_45_scale),
                                  max(rast_pred[Year %in% plot_years]$spCov_45_scale)),
                       name = "Grass:Woody\nBoundary\nStrength") +
  facet_wrap(~ Year, ncol = 8) +
  coord_equal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0,0,0,0), "cm"))
# ggsave(plot = spcov_Supplementary,
#        filename = here("GPC_Results",
#                        "VegTransitionsMap_45x45SpCov_Supplementary.jpeg"),
#        dpi = 300,
#        height = 7,
#        width = 11)

# Plot the ALL YEARS model predctions on a map (supplementary figure)
pred_maps_Supplementary <- 
  ggplot(data = na.omit(rast_pred),
         aes(x = x,
             y = y,
             fill = Probability)) +
  geom_polygon(aes(x = long, y = lat, group = group),
               fortify(riley),
               fill= "grey80") +
  geom_tile() +
  scale_fill_gradientn(colors = probRamp,
                       # values = probBreaks,
                       limits = c(min(rast_pred$Probability),
                                  max(rast_pred$Probability)),
                       name = "Relative\nProbability of\nGPC Lek\nOccurrence",
                       breaks = seq(0, 0.75, 0.25)) +
  facet_wrap(~ Year, ncol = 8) +
  coord_equal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 12),
        plot.margin = unit(c(0,0,0,0), "cm"))
# ggsave(plot = pred_maps_Supplementary,
#        filename = here("GPC_Results",
#                        "LekOccurrenceProbability_45x45SpCov_Supplementary.jpeg"),
#        dpi = 300,
#        height = 7,
#        width = 11)

routes <- spTransform(readOGR(here("Fort Riley Lek location maps"),
                              "FtRiley_Lek_SurveyRoutes"),
                      CRS(proj4string(riley)))
routes_fort <- fortify(routes) 

route_map <- 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = id),
               fortify(riley),
               fill= "grey80") +
  geom_path(aes(x = long, y = lat, group = group),
            routes_fort,
            color = "red",
            size = 0.75) +
  coord_equal() +
  theme_classic() +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 9)) +
  ylab("Northing") +
  xlab("Easting")
ggsave(plot = route_map,
       filename = here("GPC_Results", "Routes_Map.jpeg"),
       dpi = 300)
