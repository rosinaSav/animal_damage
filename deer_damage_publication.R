library(brms)
library(coda)
library(geosphere)
library(ggmap)
library(ggplot2)
library(gridExtra)
library(spdep)
library(stringr)
library(vioplot)

# Functions
######
# make a plot of the model-estimated number 
# of damaged quadrants as a function of some predictor variable
# keeping all other predictors fixed
# this is a wrapper function around means_by_predictor
make_epred_plot <- function(model, dataset, categories, variable, pretty_name, original_vector = NA, title = NA, plot_to_file = TRUE, show_legend = TRUE, enlarge = FALSE, panel_label = "") {
  file_name <- paste("figures_intensive/", title, "_", variable, ".pdf", sep = "")
  if (plot_to_file == TRUE) {
    if (variable == "edge_simple") {
      pdf(file_name, width = 4.5, height = 4)
    } else {
      pdf(file_name, width = 6, height = 4)
    }
  }
  margins <- par("mar")
  par("mar" = c(5.1, 4.6, 4.1, 2.1))
  if (enlarge == TRUE) {
    par("mar" = c(5.3, 5, 5.5, 2.4))
  }
  if (plot_to_file == TRUE) {
    par(mfrow = c(1,1))
  }
  means_by_predictor(model, dataset, categories, 
                       variable,
                       pretty_name = pretty_name, original_vector = original_vector,
                       show_legend = show_legend)
  if (plot_to_file == TRUE) {
    dev.off()
  }
  mtext(panel_label, side = 3, adj = 0, cex = 2.5)
  par("mar" = margins)
}

# make a plot of the model-estimated number 
# of damaged quadrants as a function of some predictor variable
# keeping all other predictors fixed
means_by_predictor <- function(model, dataset, categories_raw, predictor, original_vector = NA, dpar = NA, pretty_name = NA, show_legend = TRUE) {
  # scale the predictor using the mean and SD from the original data set
  if (length(original_vector) > 1) {
    original_mean <- mean(original_vector)
    original_sd <- sd(original_vector)
    categories <- (categories_raw - original_mean)/original_sd
    for_labelling <- categories_raw
  }
  else {
    categories <- categories_raw
    for_labelling <- categories
  }
  # create a synthetic data set for posterior_epred()
  new_data <- dataset[1:length(categories),]
  new_data[, predictor] <- categories
  # default values for predictors other than the focal
  # predictor
  all_cats <- list("edge_simple"="Mature vegetation", "Tree_height" = mean(dataset$Tree_height), "Distance_road" = mean(dataset$Distance_road), "Distance_water_fine" = mean(dataset$Distance_water_fine), "Distance_water_medium" = mean(dataset$Distance_water_medium), "Elevation" = mean(dataset$Elevation),
                   "Slope" = mean(dataset$Slope), "aspect_categorical" = "S")
  for (pred in names(all_cats)) {
    if (!pred == predictor) {
      new_data[, pred] <- all_cats[[pred]]
    }
  }
  ylab = "NA"
  # the conditions here check whether we're predicting based on 
  # the negative binomial component, the zero-inflated 
  # component or both
  if (is.na(dpar)) {
    pred_means <- posterior_epred(model,
                                  newdata = new_data,
                                  re_formula = NA) 
    y <- pred_means
    ylab <- "Est. # of quadrants damaged"
    ymax <- max(y) + 0.2
    main <- paste(pretty_name, "(% damaged)")
  }
  else if (dpar == "mu") {
    pred_means <- posterior_epred(model,
                                  newdata = new_data,
                                  re_formula = NA,
                                  dpar = dpar)   
    y <- pred_means * 100
    ylab <- "Estimated % damaged"
    ymax <- max(y) + 1
    main <- paste(pretty_name, "(% damaged)")
  } else if (dpar == "zi") {
    pred_means <- posterior_epred(model,
                                  newdata = new_data,
                                  re_formula = NA,
                                  dpar = dpar)
    y <- (1 - pred_means) * 100
    ylab <- "Est. probability of damage"
    ymax <- 100
    main <- paste(pretty_name, "(presence of damage)")
    
  }
  # make plot
  y_short <- y[sample(1:dim(y)[1], 500),]
  x <- seq(1, length(categories))
  trans <- rgb(173, 216, 230, alpha = 5, maxColorValue = 255)
  blue_trans <- rgb(173, 216, 230, alpha = 80, 
                    maxColorValue = 255)
  par("bty" = "n")
  plot(y[1,] ~ x,
       pch = "", xlab = "", ylab = "",
       ylim = c(0, ymax), xaxt = "n",
       main = "",
       xlim = c(0.5, length(categories) + 0.5),
       cex.axis = 1.2,
       cex.lab = 1.2,
       yaxt = "n")
  for_labelling[for_labelling == "Young pine"] <- "YP"
  for_labelling[for_labelling == "Mature vegetation"] <- "MV"
  axis(side = 1, at = x,
       labels = for_labelling, cex.axis = 2.5,
       tick = FALSE)
  axis(side = 2, at = seq(0, ymax), labels = seq(0, ymax),
       cex.axis = 2.5)
  mtext(ylab, side = 2, line = 3, cex = 1.5)
  for (column in 2:dim(y_short)[2]) {
    for (draw in 1:dim(y_short)[1]) {
      lines(c(y_short[draw, column - 1], y_short[draw, column]) ~ c(column - 1, column),
            col = blue_trans) 
    }
  }
  vioplot(y,
          names = rep("", length(for_labelling)),
          xlab = "",
          add = TRUE,
          col = trans)
  par("xpd" = TRUE)
  text(y = -3, x = seq(1, length(for_labelling)), srt = 45, 
       labels = for_labelling, cex = 1.5)
  mtext(pretty_name, side = 1, line = 3,
        cex = 1.5)
  # print out effect sizes and posterior probabilities
  counter1 <- 0
  for (cat in categories) {
    counter1 <- counter1 + 1
    counter2 <- 0
    for (cat2 in categories) {
      counter2 <- counter2 + 1
      if (!cat == cat2) {
        diff <- y[, categories == cat2] - y[, categories == cat]
        hpdi <- HPDinterval(as.mcmc(diff), 0.95)
        confidence <- mean(diff > 0) * 100
        effect_size <- round(median(diff), 3)
        # only print out effects with posterior probability > 0.9
        if (confidence > 90) {
          print(paste("Increase in ", for_labelling[counter2], " over ", for_labelling[counter1], " (", round(confidence, 1), "% prob.)", sep = ""))
          print("Posterior median:")
          print(effect_size)
          print("95% HPDI:")
          print(paste(round(hpdi[1], 3), " : ", round(hpdi[2], 3), sep = ""))
        }
      }
    }
  }
}

# make a graph of the observed distribution of quadrant damage
# overlaid with model posterior predictions
pp_check_custom <- function(model, original_vector, title, panel_label) {
  original_vector <- original_vector[!is.na(original_vector)]
  set.seed(99)
  # get posterior predictions
  ppd <- posterior_predict(model)
  trans_light_green <- rgb(144, 238, 144, alpha = 20, 
                           maxColorValue = 255)
  light_green <- rgb(144, 238, 144, 
                     maxColorValue = 255)
  # count the number of trees per damage category
  true_values <- table(original_vector)
  # sample 500 random draws from the PPD
  sampled <- sample(1:dim(ppd)[1], 500)
  # make plot
  plot(as.integer(true_values) ~ seq(0,4),
       main = title, xlab = "",
       ylim = c(0, 1500), pch = "", cex.lab = 2,
       cex.axis = 1.7, las = 2, ylab = "# of trees\n")
  mtext("# of damaged quadrants", side = 1, line = 5, cex = 1.7)
  for (line in sampled) {
    curr_table <- table(ppd[line,])
    lines(as.integer(curr_table) ~ seq(0,4),
          col = trans_light_green)
  }
  lines(as.integer(true_values) ~ seq(0,4),
        lwd = 2, col = "forestgreen") 
  legend("topright", lwd = c(2, 1), 
         col = c("forestgreen", light_green),
         legend = c("Observed", "Predicted"),
         cex = 1.5)
  mtext(panel_label, side = 3, adj = 0, cex = 1.5)
}

# check PPD for relationship with a numerical predictor
ppd_numerical_predictor <- function(model, data, panel_label, predictor, outcome, xlab, main, cuts, labels, mean = FALSE, show_legend = TRUE) {
  # cut numerical predictor into discrete categories
  categorical <- cut(data[, predictor], breaks = cuts,
                     labels = labels)
  par(xpd = TRUE)
  trans_light_green <- rgb(144, 238, 144, alpha = 20, 
                           maxColorValue = 255)
  light_green <- rgb(144, 238, 144, 
                     maxColorValue = 255)
  # get posterior predictions
  ppd <- posterior_predict(model)
  # can also be run with posterior means
  if (mean == TRUE) {
    ppd <- posterior_epred(model,
                           re_formula = NA)
  }
  # prepare plotting area
  maximum <- 3
  plot(ppd[1,] ~ as.integer(categorical),
       ylim = c(0, maximum), pch = "",
       xlab = "",
       ylab = "# of damaged quadrants", main = main,
       xaxt = "n", cex.lab = 2,
       cex.axis = 1.7, las = 2,
       yaxt = "n")
  axis(side = 2, at = seq(0,3), labels = seq(0,3), 
       cex.axis = 1.7, las = 2)
  if (predictor == "Tree_height") {
    line_number <- 5
  } else {
    line_number <- 9
  }
  mtext(xlab, side = 1, line = line_number, cex = 1.7)
  labels <- levels(categorical)
  axis(side = 1, at = seq(1, length(levels(categorical))), labels = levels(categorical),
       las = 2, cex.axis = 1.5)
  # get 500 random draws from PPD
  selected <- sample(1:nrow(ppd), 500)
  # add in lines from the PPD
  for (row in selected) {
    curr_values <- ppd[row, ]
    # summarise based on the bins as for the observed values
    agg <- aggregate(curr_values ~ categorical, FUN = "mean")
    lines(agg[,2] ~ seq(1, (length(cuts) - 1)),
          col = trans_light_green)
  }
  # summarised observed outcome based on predictor bins
  agg <- aggregate(data[, outcome] ~ categorical, FUN = "mean")
  # add line for observed data
  lines(agg[,2] ~ seq(1, (length(cuts) - 1)), lwd = 2, col = "forestgreen")
  # add legend and panel label
  if (show_legend == TRUE) {
    legend("topright", lwd = c(2, 1), 
           col = c("forestgreen", light_green),
           legend = c("Observed", "Predicted"),
           cex = 1.5)
  }
  mtext(panel_label, side = 3, adj = 0, cex = 1.5)
}

# prepare a data frame to use for building the matrix of neighbours
# for spatial correlations
prepare_df_for_sac <- function(curr_coupe, damage) {
  # filter to specific plot, coupe or location
  if (curr_coupe == "all") {
    filt <- seq(1:nrow(damage))
  } else {
    filt <- damage$Coupe == curr_coupe
  }
  # build data frame structure
  no_plots <- length(unique(damage[filt, "Plot_number"]))
  curr_df <- data.frame("plot_number" = rep(0, no_plots),
                        "longitude" = rep(0, no_plots),
                        "latitude" = rep(0, no_plots),
                        "perc_damage_below" = rep(0, no_plots),
                        "level_damage_below" = rep(0, no_plots),
                        "perc_damage_above" = rep(0, no_plots),
                        "level_damage_above" = rep(0, no_plots),
                        "coords" = rep(0, no_plots),
                        "number_damage_below" = rep(0, no_plots),
                        "number_damage_above" = rep(0, no_plots))
  # for each desired predictor, summarise it over the relevant spatial unit
  curr_df$plot_number <- aggregate(Longitude ~ Plot_number,
                                   FUN = mean,
                                   data = damage[filt,])[,1]
  curr_df$longitude <- aggregate(Longitude ~ Plot_number,
                                 FUN = mean,
                                 data = damage[filt,])[,2]
  curr_df$latitude <- aggregate(Latitude ~ Plot_number,
                                FUN = mean,
                                data = damage[filt,])[,2]
  curr_df$perc_damage_below <- aggregate(damage_below ~ Plot_number,
                                         data = damage[filt,],
                                         FUN = function(x) {mean(x > 0)})[,2]
  curr_df$level_damage_below <- aggregate(damage_below ~ Plot_number,
                                          FUN = mean,
                                          data = damage[filt,])[,2]
  curr_df$perc_damage_above <- aggregate(damage_above ~ Plot_number,
                                         FUN = function(x) {mean(x > 0)},
                                         data = damage[filt,])[,2]
  curr_df$level_damage_above <- aggregate(damage_above ~ Plot_number,
                                          FUN = mean,
                                          data = damage[filt,])[,2]
  curr_df$number_damage_below <- aggregate(damage_below ~ Plot_number,
                                           FUN = function(x) {sum(x > 0)},
                                           data = damage[filt,])[,2]
  curr_df$number_damage_above <- aggregate(damage_above ~ Plot_number,
                                           FUN = function(x) {sum(x > 0)},
                                           data = damage[filt,])[,2]
  coords <- lapply(seq(1,length(curr_df$longitude)), 
                   FUN = function(x) {c(curr_df$latitude[x], curr_df$longitude[x])})
  curr_df$coords <- coords
  curr_df <- st_as_sf(curr_df, coords = c("longitude", "latitude"),
                      crs = "+proj=longlat +ellps=WGS84")
  return(curr_df)
}
######

# Prepare data
######
damage <- read.csv("DamageAssessmentIntensiveFinal.csv")
# will be included in certain figure file names
string <- "_intensive"

# I will join categories 4 and 5 because that way, I'll be able to 
# treat this as a binomial problem
damage[, "Stem_damage_b_1m"][damage[, "Stem_damage_b_1m"] == 5] <- 4
damage[, "Stem_damage_a_1m"][damage[, "Stem_damage_a_1m"] == 5] <- 4
# make clean variables for damage
damage$damage_count <- damage[, "Stem_damage_b_1m"] + damage[, "Stem_damage_a_1m"]
damage$damage_below <- damage[, "Stem_damage_b_1m"]
damage$damage_above <- damage[, "Stem_damage_a_1m"]

# number of trees per coupe
aggregate(Tree_number ~ Coupe, FUN = length, data = damage)
######

# Data exploration
######
# relationship with # of deer/wallaby scats
# reencode scat data as presence/absence
deer_scats_categorical <- factor(damage$Deer_scats > 0)
levels(deer_scats_categorical) <- c("No scats", "Scats")
wallaby_scats_categorical <- factor(damage$Macropod_scats > 0)
levels(wallaby_scats_categorical) <- c("No scats", "Scats")
damage_below_plot <- aggregate(damage_below ~ Plot_number,
                               data = damage,
                               FUN = function(x) {mean(x > 0) * 100})[,2]
damage_above_plot <- aggregate(damage_above ~ Plot_number,
                               data = damage,
                               FUN = function(x) {mean(x > 0) * 100})[,2]
deer_scats_plot <- factor(aggregate(as.integer(deer_scats_categorical) - 1 ~ Plot_number,
                                    data = damage,
                                    FUN = max)[,2])
levels(deer_scats_plot) <- c("No scats", "Scats")
wallaby_scats_plot <- factor(aggregate(as.integer(wallaby_scats_categorical) - 1 ~ Plot_number,
                                       data = damage,
                                       FUN = max)[,2])
levels(wallaby_scats_plot) <- c("No scats", "Scats")

# plot damage based on whether scats are present
pdf(paste("figures_intensive/scats", string, ".pdf", sep = ""), height = 7.5, width = 8)
par("mfrow" = c(2,2))
par("xpd" = FALSE)
par("mar" = c(5.6, 5.1, 6, 2.1))
P <- round(wilcox.test(damage_below_plot ~ deer_scats_plot)$p.value, 2)
vioplot(damage_below_plot ~ deer_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2)
mtext(paste("Deer\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (<1m)",
      side = 2, line = 3.5, cex = 1.8)
mtext(expression(bold("A")), side = 3, adj = 0, cex = 2.2)
P <- round(wilcox.test(damage_below_plot ~ wallaby_scats_plot)$p.value, 2)
vioplot(damage_below_plot ~ wallaby_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2,
        cex.main = 1.8)
mtext(paste("Macropods\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (<1m)",
      side = 2, line = 3.5, cex = 1.8)
mtext(expression(bold("B")), side = 3, adj = 0, cex = 2.2)
P <- round(wilcox.test(damage_above_plot ~ deer_scats_plot)$p.value, 2)
vioplot(damage_above_plot ~ deer_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2,
        cex.main = 1.8)
mtext(paste("Deer\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (>1m)",
      side = 2, line = 3.5, cex = 1.8)
mtext(expression(bold("C")), side = 3, adj = 0, cex = 2.2)
P <- round(wilcox.test(damage_above_plot ~ wallaby_scats_plot)$p.value, 2)
vioplot(damage_above_plot ~ wallaby_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2,
        cex.main = 1.8)
mtext(paste("Macropods\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (>1m)",
      side = 2, line = 3.5, cex = 1.8)
mtext(expression(bold("D")), side = 3, adj = 0, cex = 2.2)
dev.off()
######

# Prepare data 2
######
# predictors that will be kept
to_keep <- c()

# simplify encoding for nearest edge
damage$edge_simple <- damage$Closest_edge
damage$edge_simple[damage$edge_simple %in% c("CL", "GR")] <- "Open"
damage$edge_simple[damage$edge_simple %in% c("NV", "MP")] <- "Mature vegetation"
damage$edge_simple[damage$edge_simple == "YP"] <- "Young pine"

to_keep <- c("edge_simple")
to_keep <- c(to_keep, "Tree_height")
damage$Wallowing <- factor(damage$Wallowing)
to_keep <- c(to_keep, "Wallowing")
to_keep <- c(to_keep, "Distance_road")
to_keep <- c(to_keep, "Distance_water_fine", "Distance_water_medium")
# encode aspect as a categorical variable
aspect_categorical <- rep("N", nrow(damage))
aspect_categorical[damage$Aspect >= 45 & damage$Aspect < 135] <- "E"
aspect_categorical[damage$Aspect >= 135  & damage$Aspect < 225] <- "S"
aspect_categorical[damage$Aspect >= 225  & damage$Aspect < 315] <- "W"
aspect_categorical <- factor(aspect_categorical, 
                             levels = c("N", "E", "S", "W"))
damage$aspect_categorical <- aspect_categorical
to_keep <- c(to_keep, "Slope")
to_keep <- c(to_keep, "Coupe", "Plot_number")

to_keep <- c(to_keep, "Aspect", "aspect_categorical", 
             "damage_count",
             "damage_below",
             "damage_above")
damage <- damage[, to_keep]
dim(damage)
# we lose 4 because of damage and because of tree height
damage <- damage[complete.cases(damage),]
# no west-facing plots
damage$aspect_categorical <- droplevels(damage$aspect_categorical)
dim(damage)
# number of observations per coupe
aggregate(edge_simple ~ Coupe, data = damage, FUN = length)

# add in GPS coordinates
damage$Latitude <- rep(0, dim(damage)[1])
damage$Longitude <- rep(0, dim(damage)[1])
damage$Elevation <- rep(0, dim(damage)[1])
gps <- read.csv("intensive_gps.tsv", sep = "\t")
gps$clean_plot_numbers <- sapply(gps$Plot_number, FUN = function(x) {y <- str_extract_all(x, "\\d+");y[[1]][2]})
for (obs in 1:nrow(damage)) {
  plot_name <- damage[obs, "Plot_number"]
  plot_number <- str_extract_all(plot_name, "\\d+")[[1]][2]
  coupe <- damage[obs, "Coupe"]
  coords <- gps[gps$clean_plot_numbers == plot_number & gps$Coupe == coupe,]
  damage[obs, "Longitude"] <- coords$Longitude
  damage[obs, "Latitude"] <- coords$Latitude
  damage[obs, "Elevation"] <- coords$Elevation
}

# the API key has been retracted for publication
API_key <- "XXX"
register_google(key = API_key)

# display the proportion of damaged trees per coupe
# on a satellite map
names_mapping <- list("coupe1" = "Coupe 1",
                      "coupe2" = "Coupe 2",
                      "coupe3" = "Coupe 3")
variable <- "damage_below"
variable_pretty <- "Prop. damaged"
pdf("figures_intensive/maps.pdf", height = 4, width = 12)
maps <- list()
for (curr_coupe in unique(damage$Coupe)) {
  damage_small <- aggregate(damage_below ~ Plot_number + Latitude + Longitude,
                            data = damage[damage$Coupe == curr_coupe,],
                            FUN = function(x) {mean(x > 0)})
  # adjust settings to each coupe
  if (curr_coupe == "coupe3") {
    zoom <- 16
    size <- c(640, 640)
    left <- 0.003
    top <- 0.002
    top2 <- top - 0.0001
    text_dist_scaling <- 0.7
  } else if (curr_coupe == "coupe2"){
    zoom <- 15
    size <- c(640, 640)
    left <- 0.004
    top <- 0.001
    top2 <- top - 0.0001
    text_dist_scaling <- 0.5
  } else {
    zoom <- 16
    size <- c(640, 640)
    left <- 0.002
    top <- 0.0008
    top2 <- top - 0.00007
    text_dist_scaling <- 0.5
  }
  # calculate the length of the scale bar by seeing what 250 m corresponds
  # to in GPS coordinates
  right <- destPoint(c((min(damage_small$Longitude) - left), (max(damage_small$Latitude) + top)), 90, 250)[1]
  colnames(damage_small)[4] <- variable_pretty
  # fetch map
  mymap <- get_map(location = c(longitude = mean(damage_small$Longitude), latitude = mean(damage_small$Latitude)), 
                   maptype = "hybrid",
                   zoom = zoom,
                   size = size)
  # add points for the observed proportion of damaged trees
  # at each plot
  mymap <- ggmap(mymap) +
    geom_point(data = damage_small, 
               aes(x = Longitude, y = Latitude, fill = !!sym(variable_pretty)), 
               size = 3, 
               shape = 21, alpha = 0.8) +
    scale_fill_gradient(low = "steelblue3", high = "salmon") +
    guides(alpha=FALSE, size=FALSE) + ggtitle(names_mapping[[curr_coupe]]) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=22),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15)) 
  
  # add scale bar
  scale_bar <- data.frame(longitude = c(min(damage_small$Longitude) - left,
                                        right,
                                        right,
                                        min(damage_small$Longitude) - left),
                          latitude = c(max(damage_small$Latitude) + top,
                                       max(damage_small$Latitude) + top,
                                       max(damage_small$Latitude) + top2,
                                       max(damage_small$Latitude) + top2))
  midpoint <- mean(c(min(damage_small$Longitude) - left, right))
  mymap <- mymap + geom_polygon(data = scale_bar, 
                                  aes(x = longitude, y = latitude),
                                fill = "white") +
    annotate("text", x = midpoint, 
             y = (max(damage_small$Latitude) + (top * text_dist_scaling)),
             label = "250 m", colour = "white")
  # add a legend to the last coupe
  if (!curr_coupe == "coupe3") {
    mymap <- mymap + theme(legend.position="none")
  }
  # adjust window size
  if (curr_coupe == "coupe2") {
    mymap <- mymap +
      scale_x_continuous(limits = c(146.667, 146.682)) +
      scale_y_continuous(limits = c(-36.638, -36.623))
  }
  maps[[curr_coupe]] <- mymap
}
grid.arrange(maps[["coupe1"]], maps[["coupe2"]], maps[["coupe3"]],
             ncol = 3)
dev.off()
######

# Model
######
# define priors
prior_b_weak <- set_prior("normal(0,5)", class = "b", dpar = "mu")
prior_b_zi_weak <- set_prior("normal(0,5)", class = "b", dpar = "zi")

prior_sd <- set_prior("normal(0,1)", class = "sd", dpar = "mu")
prior_sd_zi <- set_prior("normal(0,1)", class = "sd", dpar = "zi")
prior_sdcar <- set_prior("normal(0,1)", class = "sdcar")
prior_intercept <- set_prior("normal(0,6)", 
                             class = "Intercept",
                             dpar = "mu")
prior_intercept_zi <- set_prior("normal(0,6)", 
                                class = "Intercept",
                                dpar = "zi")

# build neighbour matrix for escar
curr_df <- prepare_df_for_sac("all", damage)
neighbour_mats <- list()
for (k in c(2, 4, 6)) {
  knn <- knearneigh(curr_df, k = k)
  plot_nb <- knn2nb(knn, row.names = curr_df$plot_number)
  names(plot_nb) <- curr_df$plot_number
  neighbour_mat <- matrix(0, nrow = length(curr_df$plot_number),
                          ncol = length(curr_df$plot_number))
  for (plot in 1:length(plot_nb)) {
    neighbour_mat[plot, plot_nb[[plot]]] <- 1
    neighbour_mat[plot_nb[[plot]], plot] <- 1
  }
  rownames(neighbour_mat) <- curr_df$plot_number
  colnames(neighbour_mat) <- curr_df$plot_number
  neighbour_mats[[k]] <- neighbour_mat
}

# overview of distribution for final data set
pdf(paste("figures_intensive/distribution_of_damage", string, ".pdf", sep = ""), width = 7, height = 3.5)
par(mfrow = c(1,2))
counts <- table(damage$damage_below)
bp <- barplot(counts ~ names(counts), ylab = "# of trees", 
              main = "<1m", xlab = "# of quadrants damaged", col = "lightblue", 
              cex.lab = 1.2, cex.main = 1.2,
              cex.axis = 1.2, las = 1, cex.names = 1.2,
              ylim = c(0, 1000))
percentages <- round(counts/sum(counts) * 100, 2)
par("xpd" = TRUE)
text(x = bp[,1],
     y = counts + 50,
     labels = paste(percentages, "%", sep = ""),
     cex = 0.8)
theta_below <- mean(damage$damage_below)/4
prob_below <- dbinom(as.numeric(names(counts)), size = 4, prob = theta_below)
counts <- table(damage$damage_above)
counts <- c(counts, 0)
names(counts)[length(names(counts))] <- "4"
bp <- barplot(counts ~ names(counts), ylab = "# of trees", 
              main = ">1m", xlab = "# of quadrants damaged", col = "lightblue",
              cex = 1, cex.lab = 1.2, cex.main = 1.2,
              cex.axis = 1.2, las = 1, cex.names = 1.2,
              ylim = c(0, 2000))
percentages <- round(counts/sum(counts) * 100, 2)
par("xpd" = TRUE)
text(x = bp[,1],
     y = counts + 80,
     labels = paste(percentages, "%", sep = ""),
     cex = 0.8)
theta_above <- mean(damage$damage_above)/4
prob_above <- dbinom(as.numeric(names(counts)), size = 4, prob = theta_above)
dev.off()

# fit model
pine_model_both_re_aspect_categorical_below <- brm(bf(damage_below|trials(4) ~ scale(Distance_water_fine) + scale(Distance_road) + edge_simple + scale(Tree_height) + scale(Slope) + aspect_categorical + scale(Elevation) + car(neighbour_mat, gr = Plot_number), zi ~ scale(Distance_water_fine) + scale(Distance_road) + edge_simple + scale(Tree_height) + scale(Slope) + aspect_categorical + scale(Elevation) + car(neighbour_mat, gr = Plot_number)),
                                                   family = zero_inflated_binomial(),
                                                   data = damage,
                                                   backend = "cmdstanr",
                                                   cores = 4,
                                                   save_pars = save_pars(all = TRUE),
                                                   seed = 5,
                                                   prior = c(prior_b_weak, prior_b_zi_weak,
                                                             prior_intercept, prior_intercept_zi,
                                                             prior_sdcar), 
                                                   data2 = list("neighbour_mat" = neighbour_mats[[6]]),
                                                   iter = 10000, warmup = 1000)

# use the medium-scale metric for distance to water
pine_model_both_re_aspect_categorical_below_medium_dist <- brm(bf(damage_below|trials(4) ~ scale(Distance_water_medium) + scale(Distance_road) + edge_simple + scale(Tree_height) + scale(Slope) + aspect_categorical + scale(Elevation) + car(neighbour_mat, gr = Plot_number), zi ~ scale(Distance_water_medium) + scale(Distance_road) + edge_simple + scale(Tree_height) + scale(Slope) + aspect_categorical + scale(Elevation) + car(neighbour_mat, gr = Plot_number)),
                                                               family = zero_inflated_binomial(),
                                                               data = damage,
                                                               backend = "cmdstanr",
                                                               cores = 4,
                                                               save_pars = save_pars(all = TRUE),
                                                               seed = 5,
                                                               prior = c(prior_b_weak, prior_b_zi_weak,
                                                                         prior_intercept, prior_intercept_zi,
                                                                         prior_sdcar), 
                                                               data2 = list("neighbour_mat" = neighbour_mats[[6]]),
                                                               iter = 10000, warmup = 1000)



# visualise the model-predicted effects of different predictors
pdf("figures_intensive/height_water_road_distance.pdf",
    width = 13.5, height = 4)
par(mfrow = c(1,3))
make_epred_plot(pine_model_both_re_aspect_categorical_below, damage, c(1, 3, 5, 7, 9, 11), "Tree_height", "Tree height (m)", original_vector = damage$Tree_height, title = paste("below", string, sep = ""),
                panel_label = "A", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below, damage, c(0, 450), "Distance_water_fine", "Distance to water (m) (fine)", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "B", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below, damage, c(35, 99), "Distance_road", "Distance to road (km)", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "C", plot_to_file = FALSE)
dev.off()

pdf("figures_intensive/height_water_road_distance_medium_scale.pdf",
    width = 13.5, height = 4)
par(mfrow = c(1,3))
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c(1, 3, 5, 7, 9, 11), "Tree_height", "Tree height (m)", original_vector = damage$Tree_height, title = paste("below", string, sep = ""),
                panel_label = "A", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c(0, 450), "Distance_water_medium", "Distance to water (medium)", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "B", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c(35, 99), "Distance_road", "Distance to road (km)", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "C", plot_to_file = FALSE)
dev.off()

pdf("figures_intensive/aspect.pdf",
    width = 9, height = 4)
par(mfrow = c(1,2))
make_epred_plot(pine_model_both_re_aspect_categorical_below, damage, c("S", "E"), "aspect_categorical", "Aspect", original_vector = NA, title = paste("below_south_west", string, sep = ""),
                panel_label = "A", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below, damage, c("N", "S"), "aspect_categorical", "Aspect", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "B", plot_to_file = FALSE)
dev.off()

# the same visualisations for the model using the medium-scale measure
# for distance to water
pdf("figures_intensive/aspect_medium.pdf",
    width = 9, height = 4)
par(mfrow = c(1,2))
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c("S", "E"), "aspect_categorical", "Aspect", original_vector = NA, title = paste("below_south_west", string, sep = ""),
                panel_label = "A", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c("N", "S"), "aspect_categorical", "Aspect", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "B", plot_to_file = FALSE)
dev.off()

pdf("figures_intensive/topographical_medium.pdf",
    width = 13.5, height = 4)
par(mfrow = c(1,3))
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c("Mature vegetation", "Young pine"), "edge_simple", "Closest edge", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "A", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c(250, 500), "Elevation", "Elevation", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "B", plot_to_file = FALSE)
make_epred_plot(pine_model_both_re_aspect_categorical_below_medium_dist, damage, c(0, 35), "Slope", "Slope", original_vector = NA, title = paste("below", string, sep = ""),
                panel_label = "C", plot_to_file = FALSE)
dev.off()

pdf("figures_intensive/PPD.pdf", height = 8, width = 8.5)
margins <- par("mar")
par("mar" = c(10.1, 6.5, 4.1, 2.1))
par(mfrow = c(2,2))
pp_check_custom(pine_model_both_re_aspect_categorical_below, damage$damage_below, "", "A")
cuts <- seq(0, 9)
labels <- paste(seq(0,8), seq(1,9), sep = "-")
ppd_numerical_predictor(pine_model_both_re_aspect_categorical_below, damage, "B", "Tree_height", "damage_below", "Tree height (m)", "", cuts = cuts, labels = labels,
                        show_legend = FALSE)
cuts <- seq(0, 180, by = 20)
labels <- paste(seq(0,160, by = 20), seq(20,180, by = 20), sep = "-")
ppd_numerical_predictor(pine_model_both_re_aspect_categorical_below, damage, "C", "Distance_road", "damage_below", "Distance to road (m)", "", cuts = cuts, labels = labels,
                        show_legend = FALSE)
cuts <- seq(0, 325, by = 25)
labels <- paste(seq(0,300, by = 25), seq(25,325, by = 25), sep = "-")
ppd_numerical_predictor(pine_model_both_re_aspect_categorical_below, damage, "D", "Distance_water_medium", "damage_below", "Distance to water (m) (fine)", "", cuts = cuts, labels = labels,
                        show_legend = FALSE)
dev.off()
######