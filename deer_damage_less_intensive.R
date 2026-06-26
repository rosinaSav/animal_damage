library(brms)
library(coda)
library(geosphere)
library(ggmap)
library(ggplot2)
library(spdep)
library(stringr)
library(vioplot)

# Functions
######
# return the lower boundary of a 95% HPDI
HPDIlow <- function(x) {
  return(HPDinterval(as.mcmc(x), prob = 0.95)[1])
}

# return the upper boundary of a 95% HPDI
HPDIhigh <- function(x) {
  return(HPDinterval(as.mcmc(x), prob = 0.95)[2])
}

# make a bar plot of the distirbution of trees
# into different categories of quadrant damage
make_distribution_plot <- function(data_in, title, to_add) {
  counts <- table(data_in)
  bp <- barplot(counts ~ names(counts), ylab = "# of trees\n", 
                main = title, xlab = "# of quadrants damaged", col = "lightblue", 
                cex.lab = 1.5, cex.main = 1.5,
                cex.axis = 1.5, las = 1, cex.names = 1.5,
                ylim = c(0, 1200))
  # calculate percentages to add to the graph
  percentages <- round(counts/sum(counts) * 100, 2)
  par("xpd" = TRUE)
  text(x = bp[,1],
  y = counts + to_add,
  labels = paste(percentages, "%", sep = ""),
  cex = 0.8)
  theta_below <- mean(damage$damage_below)/4
  prob_below <- dbinom(as.numeric(names(counts)), size = 4, prob = theta_below)
}

# make a plot of the model-estimated number 
# of damaged quadrants as a function of some predictor variable
# keeping all other predictors fixed
# this is a wrapper function around means_by_predictor/
# means_by_predictor_by_location_and_age/
# means_by_predictor_by_location
make_epred_plot <- function(model, dataset, categories, variable, pretty_name, original_vector = NA, title = NA, by_location = FALSE, by_location_age = FALSE, plot_to_file = TRUE, show_legend = TRUE, enlarge = FALSE, panel_label = "") {
  file_name <- paste("figures_less_intensive/", title, "_", variable, ".pdf", sep = "")
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
  if (by_location == TRUE) {
    means_by_predictor_by_location(model, dataset, categories, 
                       variable,
                       pretty_name = pretty_name, original_vector = original_vector,
                       show_legend = show_legend)
  } else if (by_location_age == TRUE) {
    par("xpd" = TRUE)
    means_by_predictor_by_location_and_age(model, dataset, categories, 
                                   variable,
                                   pretty_name = pretty_name, original_vector = original_vector,
                                   show_legend = show_legend, 
                                   enlarge = enlarge,
                                   panel_label = panel_label,
                                   plot_to_file = plot_to_file)
  }
    else {
    means_by_predictor(model, dataset, categories, 
                                   variable,
                                   pretty_name = pretty_name, original_vector = original_vector,
                       show_legend = show_legend)
    }
  if (plot_to_file == TRUE) {
    dev.off()
  }
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
  all_cats <- list("edge_simple"="Mature vegetation", "Tree_height" = mean(dataset$Tree_height), "Age" = "2", "Density_road" = mean(dataset$Density_road), "Distance_water_fine" = mean(dataset$Distance_water_fine), "Distance_water_medium" = mean(dataset$Distance_water_medium), "Elevation" = mean(dataset$Elevation),
                     "Slope" = mean(dataset$Slope), "aspect_categorical" = "S", "Location" = "Merriang")
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
       pch = "", xlab = "", ylab = ylab,
       ylim = c(0, ymax), xaxt = "n",
       main = "",
       xlim = c(0.5, length(categories) + 0.5),
       cex.axis = 1.2,
       cex.lab = 1.2)
  axis(side = 1, at = x,
       labels = for_labelling, cex.axis = 1.7)
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

# similar to means_by_predictor() but making separate predictions for
# the two locations
# see inline comments in means_by_predictor() for further information
means_by_predictor_by_location <- function(model, dataset, categories_raw, predictor, original_vector = NA, dpar = NA, pretty_name = NA, show_legend = TRUE) {
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
  new_data <- dataset[1:(length(categories) * 2),]
  new_data[, predictor] <- rep(categories, 2)
  new_data$Location <- c(rep("Merriang", length(categories)), rep("Bright", length(categories)))
  all_cats <- list("edge_simple"="Mature vegetation", "Tree_height" = mean(dataset$Tree_height), "Age" = "2", "Density_road" = mean(dataset$Density_road), "Distance_water_fine" = mean(dataset$Distance_water_fine), "Distance_water_medium" = mean(dataset$Distance_water_medium), "Elevation" = mean(dataset$Elevation),
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
  
  y_short <- y[sample(1:dim(y)[1], 500),]
  x <- seq(1, length(categories) * 2)
  trans <- rgb(173, 216, 230, alpha = 5, maxColorValue = 255)
  blue_trans <- rgb(173, 216, 230, alpha = 80, 
                    maxColorValue = 255)
  desert_trans <- rgb(250, 213, 165, alpha = 80, 
                      maxColorValue = 255)
  par("bty" = "n")
  plot(y[1,] ~ x,
       pch = "", xlab = "", ylab = ylab,
       ylim = c(0, ymax), xaxt = "n",
       main = "",
       xlim = c(0.5, (length(categories) * 2) + 0.5),
       cex.lab = 1.7,
       cex.axis = 1.7, las = 2,
       yaxt = "n")
  axis(side = 2, at = seq(0,4), labels = seq(0,4), 
       cex.axis = 1.7, las = 2)
  text(y = -5, x = c(1.5, 3.5), 
       labels = c("Merriang", "Bright"), cex = 1.7)
  for (column in 2:length(for_labelling)) {
    for (draw in 1:nrow(y_short)) {
      lines(c(y_short[draw, column - 1], y_short[draw, column]) ~ c(column - 1, column),
            col = blue_trans) 
    }
  }
  for (column in (length(for_labelling) + 2):(length(for_labelling) * 2)) {
    for (draw in 1:nrow(y_short)) {
      lines(c(y_short[draw, column - 1], y_short[draw, column]) ~ c(column - 1, column),
            col = desert_trans) 
    }
  }
  vioplot(y,
          names = rep("", length(for_labelling) * 2),
          xlab = "",
          add = TRUE,
          col = trans)
  par("xpd" = TRUE)
  mtext(side = 1, line = 1, text = rep(for_labelling, 2),
        cex = 1.7, at = seq(1, length(for_labelling) * 2))
  mtext(pretty_name, side = 1, line = 3,
        cex = 1.7)
  mtext(side = 3, line = 1, text = ~italic("Merriang"),
        cex = 1.8, at = c(1.5))
  mtext(side = 3, line = 1, text = ~italic("Bright"),
        cex = 1.8, at = c(3.5))
  counter1 <- 0
  for (cat in categories) {
    counter1 <- counter1 + 1
    counter2 <- 0
    for (cat2 in categories) {
      counter2 <- counter2 + 1
      if (!cat == cat2) {
        diff_merriang <- y[, new_data$Location == "Merriang" & new_data[,predictor] == cat2] - y[, new_data$Location == "Merriang" & new_data[,predictor] == cat]
        diff_bright <- y[, new_data$Location == "Bright" & new_data[,predictor] == cat2] - y[, new_data$Location == "Bright" & new_data[,predictor] == cat]
        diff <- diff_bright - diff_merriang
        prob <- mean(diff > 0) * 100
        prob_merriang <- mean(diff_merriang > 0) * 100
        prob_bright <- mean(diff_bright > 0) * 100
        
        if (prob_merriang > 50) {
          print(paste("Increase in ", for_labelling[counter2], " over ", for_labelling[counter1], " in Merriang (", round(prob_merriang, 1), "% prob.)", sep = ""))
          print("Posterior median:")
          effect_size <- round(median(diff_merriang), 3)
          print(effect_size)
          print("95% HPDI:")
          hpdi <- HPDinterval(as.mcmc(diff_merriang), 0.95)
          print(paste(round(hpdi[1], 3), " : ", round(hpdi[2], 3), sep = ""))
        }
        if (prob_bright > 50) {
          print(paste("Increase in ", for_labelling[counter2], " over ", for_labelling[counter1], " in Bright (", round(prob_bright, 1), "% prob.)", sep = ""))
          print("Posterior median:")
          effect_size <- round(median(diff_bright), 3)
          print(effect_size)
          print("95% HPDI:")
          hpdi <- HPDinterval(as.mcmc(diff_bright), 0.95)
          print(paste(round(hpdi[1], 3), " : ", round(hpdi[2], 3), sep = ""))
        }
        if (prob > 50) {
          print(paste("Increase in Bright over Merriang (", round(prob, 1), "% prob.)", sep = ""))
          print("Posterior median:")
          effect_size <- round(median(diff), 3)
          print(effect_size)
          print("95% HPDI:")
          hpdi <- HPDinterval(as.mcmc(diff), 0.95)
          print(paste(round(hpdi[1], 3), " : ", round(hpdi[2], 3), sep = ""))
        }
      }
    }
  }
}

# similar to means_by_predictor() but making separate predictions for
# the two locations and the two age classes
# see inline comments in means_by_predictor() for further information
means_by_predictor_by_location_and_age <- function(model, dataset, categories_raw, predictor, original_vector = NA, dpar = NA, pretty_name = NA, show_legend = TRUE, enlarge = FALSE, panel_label = "", plot_to_file = TRUE) {
  par("xpd" = TRUE)
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
  new_data <- dataset[1:(length(categories) * 4),]
  new_data[, predictor] <- rep(categories, 4)
  new_data$Location <- c(rep("Merriang", length(categories) * 2), rep("Bright", length(categories) * 2))
  new_data$Age <- rep(c(rep("2", length(categories)), rep("4", length(categories))), 2)
  all_cats <- list("edge_simple"="Mature vegetation", "Tree_height" = mean(dataset$Tree_height), "Density_road" = mean(dataset$Density_road), "Distance_water_fine" = mean(dataset$Distance_water_fine), "Distance_water_medium" = mean(dataset$Distance_water_medium), "Elevation" = mean(dataset$Elevation),
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
  
  y_short <- y[sample(1:dim(y)[1], 100),]
  x <- seq(1, length(categories) * 4)
  trans <- rgb(173, 216, 230, alpha = 5, maxColorValue = 255)
  dark_blue <- rgb(49, 75, 162, 
                   maxColorValue = 255)
  darkblue_trans <- rgb(49, 75, 162, alpha = 120, 
                    maxColorValue = 255)
  darkblue_lines <- rgb(11, 32, 102,
                        maxColorValue = 255)
  light_blue <- rgb(154, 191, 213, 
                    maxColorValue = 255)
  lightblue_trans <- rgb(154, 191, 213, alpha = 120, 
                         maxColorValue = 255)
  lightblue_lines <- rgb(46, 153, 214,
                        maxColorValue = 255)
  dark_desert <- rgb(207, 153, 10, 
                   maxColorValue = 255)
  darkdesert_trans <- rgb(245, 0, 0, alpha = 120, 
                         maxColorValue = 255)
  darkdesert_lines <- rgb(225, 20, 20,
                        maxColorValue = 255)
  light_desert <- rgb(247, 246, 185, 
                    maxColorValue = 255)
  lightdesert_trans <- rgb(255, 206, 27, alpha = 120, 
                         maxColorValue = 255)
  lightdesert_lines <- rgb(211, 165, 26,
                         maxColorValue = 255)
  par("bty" = "n")
  if (enlarge == TRUE) {
    axis_size <- 3
    text_size <- 2.1
    text_line <- 4
    final_ylab <- "# of quadrants"
    legend_bty <- "n"
  } else {
    axis_size <- 1.7
    text_size <- 1.5
    text_line <- 3
    final_ylab <- ylab
    legend_bty <- "o"
    if (plot_to_file == "FALSE") {
      legend_bty <- "n"
      axis_size <- 2
      text_size <- 1.7
      final_ylab <- "# of quadrants"
    }
  }
  plot(y[1,] ~ x,
       pch = "", xlab = "", ylab = final_ylab,
       ylim = c(0, ymax), xaxt = "n",
       main = "",
       xlim = c(0.5, (length(categories) * 2) + 0.5),
       cex.lab = axis_size,
       cex.axis = axis_size, las = 2,
       yaxt = "n")
  axis(side = 2, at = seq(0,4), labels = seq(0,4), 
       cex.axis = axis_size, las = 2)
  # text(y = -5, x = c(mean(c(1, length(categories) * 2)), mean(c(length(categories) * 2 + 1, length(categories * 4)))), 
  #      labels = c("Merriang", "Bright"), cex = 1.5)
  medians <- apply(y, MARGIN = 2, FUN = median)
  HPDI_low <- apply(y, MARGIN = 2, FUN = HPDIlow)
  HPDI_high <- apply(y, MARGIN = 2, FUN = HPDIhigh)
  polygon(x = c(seq(1:length(for_labelling)), rev(seq(1:length(for_labelling)))),  
       y = c(HPDI_low[1:length(for_labelling)], rev(HPDI_high[1:length(for_labelling)])), 
       border = NA, col = lightblue_trans)
  lines(medians[1:length(for_labelling)] ~ seq(1,length(for_labelling)), col = lightblue_lines, lwd = 3)
  curr_seq <- (length(for_labelling) + 1):(length(for_labelling) * 2)
  polygon(x = c(seq(1,length(for_labelling)), rev(seq(1,length(for_labelling)))),  
          y = c(HPDI_low[curr_seq], rev(HPDI_high[curr_seq])), 
          border = NA, col = darkblue_trans)
  lines(medians[curr_seq] ~ seq(1,length(for_labelling)), col = darkblue_lines, lwd = 3)
  curr_seq <- ((length(for_labelling) * 2) + 1):(length(for_labelling) * 3)
  polygon(x = c(seq((length(for_labelling)) + 1,(length(for_labelling) * 2)), rev(seq((length(for_labelling)) + 1,(length(for_labelling) * 2)))),  
          y = c(HPDI_low[curr_seq], rev(HPDI_high[curr_seq])), 
          border = NA, col = lightdesert_trans)
  lines(medians[((length(for_labelling) * 2) + 1):(length(for_labelling) * 3)] ~ seq((length(for_labelling)) + 1,(length(for_labelling) * 2)), col = lightdesert_lines, lwd = 3)
  curr_seq <- ((length(for_labelling) * 3) + 1):(length(for_labelling) * 4)
  polygon(x = c(seq(((length(for_labelling)) + 1),(length(for_labelling) * 2)), rev(seq(((length(for_labelling)) + 1),(length(for_labelling) * 2)))),  
          y = c(HPDI_low[curr_seq], rev(HPDI_high[curr_seq])), 
          border = NA, col = darkdesert_trans)
  lines(medians[((length(for_labelling) * 3) + 1):(length(for_labelling) * 4)] ~ seq(((length(for_labelling)) + 1),(length(for_labelling) * 2)), col = darkdesert_lines, lwd = 3)
  par("xpd" = TRUE)
  if (predictor == "edge_simple") {
    for_labelling <- c("MV", "YP", "O")
  }
  mtext(side = 1, line = 1, text = rep(for_labelling, 2),
        cex = text_size, at = seq(1, length(for_labelling) * 2))
  mtext(pretty_name, side = 1, line = text_line,
        cex = text_size)
  if (show_legend == TRUE) {
    if (enlarge == TRUE) {
      ypos <- 6
      xpos <- 1
      legend_cex <- 3
    } else {
      ypos <- 6
      xpos <- 0.5
      legend_cex <- 1.2
      if (plot_to_file == FALSE) {
        xpos <- 0.7
        legend_cex <- 1.8
      }
    }
    
    legend_text <- c("2 yrs Merriang", "4 yrs Merriang", "2 yrs Bright", "4 yrs Bright")
    if (enlarge == TRUE) {
      legend_text <- c("2 yrs M.", "4 yrs M.", "2 yrs B.", "4 yrs B.")
    }
    if (plot_to_file == FALSE) {
      legend_text <- c("2 yrs M.", "4 yrs M.", "2 yrs B.", "4 yrs B.")
    }
    legend(x = xpos, y = ypos, legend = legend_text,
           fill = c(lightblue_lines, darkblue_lines, lightdesert_lines, darkdesert_lines),
           ncol = 2, cex = legend_cex, bty = legend_bty) 
  }
  mtext(panel_label, side = 3, adj = 0, cex = 2)
  counter1 <- 0
  # print out comparisons
  if (predictor == "Tree_height") {
    categories <- c(min(categories), max(categories))
  }
  for (cat in categories) {
    counter1 <- counter1 + 1
    counter2 <- 0
    for (cat2 in categories) {
      counter2 <- counter2 + 1
      if (!cat == cat2) {
        print("*******************")
        print(paste(cat2, "minus", cat))
        print(cat)
        print(cat2)
        diff_merriang2 <- y[, new_data$Location == "Merriang" & (new_data[,predictor] == cat2 & new_data[, "Age"] == "2")] - y[, new_data$Location == "Merriang" & (new_data[,predictor] == cat & new_data[, "Age"] == "2")]
        diff_bright2 <- y[, new_data$Location == "Bright" & (new_data[,predictor] == cat2 & new_data[, "Age"] == "2")] - y[, new_data$Location == "Bright" & (new_data[,predictor] == cat & new_data[, "Age"] == "2")]
        diff_merriang4 <- y[, new_data$Location == "Merriang" & (new_data[,predictor] == cat2 & new_data[, "Age"] == "4")] - y[, new_data$Location == "Merriang" & (new_data[,predictor] == cat & new_data[, "Age"] == "4")]
        diff_bright4 <- y[, new_data$Location == "Bright" & (new_data[,predictor] == cat2 & new_data[, "Age"] == "4")] - y[, new_data$Location == "Bright" & (new_data[,predictor] == cat & new_data[, "Age"] == "4")]
        diff_location2 <- diff_bright2 - diff_merriang2
        diff_location4 <- diff_bright4 - diff_merriang4
        diff_agebright <- diff_bright4 - diff_bright2
        diff_agemerriang <- diff_merriang4 - diff_merriang2
        all_diffs <- list("merriang2" = diff_merriang2,
                          "bright2" = diff_bright2,
                          "merriang4" = diff_merriang4,
                          "bright4" = diff_bright4,
                          "location2" = diff_location2,
                          "location4" = diff_location4,
                          "agebright" = diff_agebright,
                          "agemerriang" = diff_agemerriang)
        for (diff_name in names(all_diffs)) {
          diff <- all_diffs[[diff_name]]
          prob <- round(mean(diff > 0) * 100, 2)
          prob_neg <- 100 - prob
          print(diff_name)
          print("Probability of positive/negative effect:")
          print(paste(prob, "%/", prob_neg, "%", sep = ""))
          effect_size <- round(median(diff), 3)
          print("Effect size:")
          print(effect_size)
          print("95% HPDI:")
          hpdi <- HPDinterval(as.mcmc(diff), 0.95)
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
prepare_df_for_sac <- function(curr_coupe, damage, parameter_large, parameter_narrow, parameter_narrowest) {
  # filter to specific plot, coupe or location
  if (curr_coupe == "all") {
    filt <- seq(1:nrow(damage))
  } else {
    filt <- damage[, parameter_large] == curr_coupe
  }
  # build data frame structure
  no_plots <- length(unique(damage[filt, parameter_narrow]))
  curr_df <- data.frame("coupe_number" = rep(0, no_plots),
                        "longitude" = rep(0, no_plots),
                        "latitude" = rep(0, no_plots),
                        "coords" = rep(0, no_plots),
                        "mean_damage_below" = rep(0, no_plots),
                        "mean_damage_above" = rep(0, no_plots))
  # for each desired predictor, summarise it over the relevant spatial unit
  if (parameter_narrowest == FALSE) {
    aggregate_by_plot_below <- aggregate(damage[filt, "damage_below"] ~ damage[filt, parameter_narrow],
                                         FUN = function(x) {mean(log(x + 1))})
    aggregate_by_plot_above <- aggregate(damage[filt, "damage_above"] ~ damage[filt, parameter_narrow],
                                         FUN = function(x) {mean(log(x + 1))})
    curr_df$mean_damage_below <- aggregate(aggregate_by_plot_below[,2] ~ aggregate_by_plot_below[,1],
                                           FUN = function(x) {exp(mean(x))})[,2]
    curr_df$mean_damage_above <- aggregate(aggregate_by_plot_above[,2] ~ aggregate_by_plot_above[,1],
                                           FUN = function(x) {exp(mean(x))})[,2]
  } else {
    aggregate_by_plot_below <- aggregate(damage[filt, "damage_below"] ~ damage[filt, parameter_narrow] + damage[filt, parameter_narrowest],
                                         FUN = function(x) {mean(log(x + 1))})
    aggregate_by_plot_above <- aggregate(damage[filt, "damage_above"] ~ damage[filt, parameter_narrow] + damage[filt, parameter_narrowest],
                                         FUN = function(x) {mean(log(x + 1))})
    curr_df$mean_damage_below <- aggregate(aggregate_by_plot_below[,3] ~ aggregate_by_plot_below[,1],
                                           FUN = function(x) {exp(mean(x))})[,2]
    curr_df$mean_damage_above <- aggregate(aggregate_by_plot_above[,3] ~ aggregate_by_plot_above[,1],
                                           FUN = function(x) {exp(mean(x))})[,2]
  }
  curr_df$plot_number <- aggregate(damage[filt, "Longitude"] ~ damage[filt, parameter_narrow],
                                   FUN = mean)[,1]
  curr_df$longitude <- aggregate(damage[filt, "Longitude"] ~ damage[filt, parameter_narrow],
                                 FUN = mean)[,2]
  curr_df$latitude <- aggregate(damage[filt, "Latitude"] ~ damage[filt, parameter_narrow],
                                FUN = mean)[,2]
  curr_df$Age <- aggregate(damage[filt, "Age"] ~ damage[filt, parameter_narrow],
                                FUN = unique)[,2]
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
damage <- read.csv("DamageAssessLess_IntensiveFinal.csv")
# encode age 5 as 4
damage$Age[damage$Age == 5] <- 4

# remove two coupes where height/DBH measurements were found to be unreliable
# the coupe names have been anonymised for publication
damage <- damage[!damage$HVP_Coupe %in% c("coupe1", "coupe2"),]
# count trees per coupe
aggregate(Coupe ~ HVP_Coupe, data = damage, FUN = length)

# reformatting of plot names
plot_number <- c()
for (plot in damage$Plot_number) {
  split <- str_split(plot, "P")[[1]]
  plot_number <- c(plot_number, split[length(split)])
}
damage$Plot_number <- plot_number
damage <- damage[damage$Plot_number != "",]
damage$Plot_number_specific <- paste(damage$Coupe, damage$Plot_number, sep = "_")

# I will join categories 4 and 5 (ring-barked) because that way, I'll be able to 
# treat this as a binomial problem
damage[, "Stem_damage._..1m"][damage[, "Stem_damage._..1m"] == 5] <- 4
damage[, "Stem_damage.1m"][damage[, "Stem_damage.1m"] == 5] <- 4
# make clean variables for damage
damage$damage_count <- damage[, "Stem_damage._..1m"] + damage[, "Stem_damage.1m"]
damage$damage_below <- damage[, "Stem_damage._..1m"]
damage$damage_above <- damage[, "Stem_damage.1m"]

# predictors that will be kept
to_keep <- c()

# simplify encoding for nearest edge
damage$edge_simple <- damage$Closest_edge
damage$edge_simple[damage$edge_simple %in% c("CL", "GR")] <- "Open"
damage$edge_simple[damage$edge_simple %in% c("NV", "MP")] <- "Mature vegetation"
damage$edge_simple[damage$edge_simple == "YP"] <- "Young pine"

damage$Age <- factor(damage$Age)
to_keep <- c(to_keep, "Age")
damage$Wallowing <- factor(damage$Wallowing)
to_keep <- c(to_keep, c("Wallowing", "edge_simple", "Tree_height"))
to_keep <- c(to_keep, "Density_road")
to_keep <- c(to_keep, "Distance_water_fine", "Distance_water_medium")
to_keep <- c(to_keep, "Elevation")
# reencode aspect as categorical
aspect_categorical <- rep("N", nrow(damage))
aspect_categorical[damage$Aspect >= 45 & damage$Aspect < 135] <- "E"
aspect_categorical[damage$Aspect >= 135  & damage$Aspect < 225] <- "S"
aspect_categorical[damage$Aspect >= 225  & damage$Aspect < 315] <- "W"
aspect_categorical <- factor(aspect_categorical, 
                             levels = c("N", "E", "S", "W"))
damage$aspect_categorical <- aspect_categorical

to_keep <- c(to_keep, "Slope")

to_keep <- c(to_keep, "Coupe", "Plot_number_specific", "Location", "Plot_number", "HVP_Coupe")

to_keep <- c(to_keep, "Aspect", "aspect_categorical",
             "damage_below",
             "damage_above")
######

# Data exploration 
######
# relationship with DBH
pdf("figures_less_intensive/tree_height_pattern.pdf",
    width = 6, height = 3)
par(mfrow = c(1,2))
vioplot(Tree_height ~ Age, data = damage,
     col = c("lightblue", "steelblue3"), main = "",
     las = 2,
     cex.lab = 1.2, cex.axis = 1.2,
     xlab = "", ylab = "")
mtext("A", side = 3, adj = 0, cex = 1.5)
mtext("Age (years)", side = 1, line = 3,
      cex = 1.2)
mtext("Height (m)", side = 2, line = 3,
      cex = 1.2)
grey_trans <- rgb(140, 140, 140, alpha = 40,
                  maxColorValue = 255)
plot(Tree_height ~ DBH, data = damage,
     col = grey_trans, main = "",
     ylab = "Height (m)", pch = 20,
     xlab = "DBH (cm)",
     cex.lab = 1.2, cex.axis = 1.2,
     las = 2)
mtext("B", side = 3, adj = 0, cex = 1.5)
dev.off()

corr <- cor.test(damage$Tree_height, damage$DBH)

# relationship between DBH and damage
par(mfrow = c(1,3))
vioplot(damage$DBH ~ damage$damage_below, col = "skyblue",
        xlab = "# of damaged quadrants",
        ylab = "DBH", main = "All trees",
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
vioplot(damage$DBH[damage$Age == "2"] ~ damage$damage_below[damage$Age == "2"], col = "skyblue",
        xlab = "# of damaged quadrants",
        ylab = "DBH", main = "Age 2 only",
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
vioplot(damage$DBH[damage$Age == "4"] ~ damage$damage_below[damage$Age == "4"], col = "skyblue",
        xlab = "# of damaged quadrants",
        ylab = "DBH", main = "Age 4 only",
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)

# relationship with # of deer/wallaby scats
# reencode scat data as presence/absence
deer_scats_categorical <- factor(damage$Deer_scats > 0)
levels(deer_scats_categorical) <- c("No scats", "Scats")
wallaby_scats_categorical <- factor(damage$Macropod_scats > 0)
levels(wallaby_scats_categorical) <- c("No scats", "Scats")

damage_below_plot <- aggregate(damage_below ~ Plot_number_specific,
                               data = damage,
                               FUN = function(x) {mean(x > 0) * 100})[,2]
damage_above_plot <- aggregate(damage_above ~ Plot_number_specific,
                               data = damage,
                               FUN = function(x) {mean(x > 0) * 100})[,2]
deer_scats_plot <- factor(aggregate(as.integer(deer_scats_categorical) - 1 ~ Plot_number_specific,
                                    data = damage,
                                    FUN = max)[,2])
levels(deer_scats_plot) <- c("No scats", "Scats")
wallaby_scats_plot <- factor(aggregate(as.integer(wallaby_scats_categorical) - 1 ~ Plot_number_specific,
                                       data = damage,
                                       FUN = max)[,2])
levels(wallaby_scats_plot) <- c("No scats", "Scats")

blue <- rgb(173, 216, 230, 
            maxColorValue = 255)
desert <- rgb(250, 213, 165, 
              maxColorValue = 255)

# are deer and wallaby scats present at the same plots?
deer_merriang_table <- aggregate(Deer_scats ~ Plot_number_specific,
          data = damage[damage$Location == "Merriang",],
          FUN = function(x) {mean(x > 0)})
wallabies_merriang_table <- aggregate(Macropod_scats ~ Plot_number_specific,
                                 data = damage[damage$Location == "Merriang",],
                                 FUN = function(x) {mean(x > 0)})
t1 <- table(deer_merriang_table[,2], wallabies_merriang_table[,2])
t1/sum(t1)

deer_bright_table <- aggregate(Deer_scats ~ Plot_number_specific,
                                 data = damage[damage$Location == "Bright",],
                                 FUN = function(x) {mean(x > 0)})
wallabies_bright_table <- aggregate(Macropod_scats ~ Plot_number_specific,
                                      data = damage[damage$Location == "Bright",],
                                      FUN = function(x) {mean(x > 0)})
t2 <- table(deer_bright_table[,2], wallabies_bright_table[,2])
t2/sum(t2)

deer_merriang <- mean(aggregate(Deer_scats ~ Plot_number_specific,
          data = damage[damage$Location == "Merriang",],
          FUN = function(x) {mean(x > 0)})[,2]) * 100
deer_bright <-mean(aggregate(Deer_scats ~ Plot_number_specific,
               data = damage[damage$Location == "Bright",],
               FUN = function(x) {mean(x > 0)})[,2]) * 100
wallabies_merriang <- mean(aggregate(Macropod_scats ~ Plot_number_specific,
               data = damage[damage$Location == "Merriang",],
               FUN = function(x) {mean(x > 0)})[,2]) * 100
wallabies_bright <- mean(aggregate(Macropod_scats ~ Plot_number_specific,
               data = damage[damage$Location == "Bright",],
               FUN = function(x) {mean(x > 0)})[,2]) * 100

# plot scat presence based on location
pdf("figures_less_intensive/scats_by_location.pdf",
    height = 4, width = 5)
par(mfrow = c(1,1))
margins <- par("mar")
par("mar" = c(6.1, 5.1, 4.1, 2.1))
bp <- barplot(c(deer_merriang, deer_bright, wallabies_merriang, wallabies_bright),
        col = c(blue, desert, blue, desert),
        xlab = "",
        ylab = "% of plots wih scats",
        space = c(0.2, 0.2, 1, 0.2), las = 2, 
        cex.lab = 1.7,
        cex.names = 1.7, cex.axis = 1.7)
mtext(text = "Deer", line = 1, side = 1, at = mean(bp[1:2]),
      cex = 1.7)
mtext(text = "Macropod", line = 1, side = 1, at = mean(bp[3:4]),
      cex = 1.7)
legend("topleft", fill = c(blue, desert),
       legend = c("Merriang", "Bright"),
       cex = 1.5, bty = "n")
dev.off()
par("mar" = margins)

# add age
deer_merriang2 <- mean(aggregate(Deer_scats ~ Plot_number_specific,
                                data = damage[damage$Location == "Merriang" & damage$Age == 2,],
                                FUN = function(x) {mean(x > 0)})[,2]) * 100
deer_merriang4 <- mean(aggregate(Deer_scats ~ Plot_number_specific,
                                 data = damage[damage$Location == "Merriang" & damage$Age == 4,],
                                 FUN = function(x) {mean(x > 0)})[,2]) * 100
deer_bright2 <- mean(aggregate(Deer_scats ~ Plot_number_specific,
                             data = damage[damage$Location == "Bright" & damage$Age == 2,],
                             FUN = function(x) {mean(x > 0)})[,2]) * 100
deer_bright4 <- mean(aggregate(Deer_scats ~ Plot_number_specific,
                               data = damage[damage$Location == "Bright" & damage$Age == 4,],
                               FUN = function(x) {mean(x > 0)})[,2]) * 100
wallabies_merriang2 <- mean(aggregate(Macropod_scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Merriang" & damage$Age == 2,],
                                     FUN = function(x) {mean(x > 0)})[,2]) * 100
wallabies_merriang4 <- mean(aggregate(Macropod_scats ~ Plot_number_specific,
                                      data = damage[damage$Location == "Merriang" & damage$Age == 4,],
                                      FUN = function(x) {mean(x > 0)})[,2]) * 100
wallabies_bright2 <- mean(aggregate(Macropod_scats ~ Plot_number_specific,
                                   data = damage[damage$Location == "Bright" & damage$Age == 2,],
                                   FUN = function(x) {mean(x > 0)})[,2]) * 100
wallabies_bright4 <- mean(aggregate(Macropod_scats ~ Plot_number_specific,
                                    data = damage[damage$Location == "Bright" & damage$Age == 4,],
                                    FUN = function(x) {mean(x > 0)})[,2]) * 100

pdf("figures_less_intensive/scats_by_location_and_age.pdf",
    height = 4, width = 6)
dark_blue <- rgb(49, 75, 162, alpha = 120, 
                      maxColorValue = 255)
light_blue <- rgb(154, 191, 213, alpha = 120, 
                       maxColorValue = 255)
dark_desert <- rgb(245, 0, 0, alpha = 120, 
                        maxColorValue = 255)
light_desert <- rgb(255, 206, 27, alpha = 120, 
                         maxColorValue = 255)

par(mfrow = c(1,1))
margins <- par("mar")
par("xpd" = TRUE)
par("mar" = c(6.1, 5.1, 4.1, 2.1))
bp <- barplot(c(deer_merriang2, deer_merriang4, deer_bright2, deer_bright4, wallabies_merriang2, wallabies_merriang4, wallabies_bright2, wallabies_bright4),
              col = c(light_blue, dark_blue, light_desert, dark_desert, light_blue, dark_blue, light_desert, dark_desert),
              xlab = "",
              ylab = "% of plots wih scats",
              space = c(0.2, 0.2, 0.5, 0.2, 1, 0.2, 0.5, 0.2), las = 2, 
              cex.lab = 1.7,
              cex.names = 1.7, cex.axis = 1.7)
mtext(text = ~italic("Merriang"), line = 1, side = 1, at = mean(bp[1:2]),
      cex = 1.7)
mtext(text = ~italic("Bright"), line = 1, side = 1, at = mean(bp[3:4]),
      cex = 1.7)
mtext(text = ~italic("Merriang"), line = 1, side = 1, at = mean(bp[5:6]),
      cex = 1.7)
mtext(text = ~italic("Bright"), line = 1, side = 1, at = mean(bp[7:8]),
      cex = 1.7)
mtext(text = "Deer", line = 3, side = 1, at = mean(bp[1:4]),
      cex = 1.7)
mtext(text = "Macropod", line = 3, side = 1, at = mean(bp[5:8]),
      cex = 1.7)
legend(x = -1, y = 85, legend = c("2 yrs Merriang", "4 yrs Merriang", "2 yrs Bright", "4 yrs Bright"),
       fill = c(light_blue, dark_blue, light_desert, dark_desert),
       ncol = 2, cex = 1.7, bty = "n")
dev.off()
par("mar" = margins)

deer_merriang_count <- sum(aggregate(Deer_scats ~ Plot_number_specific,
                                 data = damage[damage$Location == "Merriang",],
                                 FUN = function(x) {mean(x > 0)})[,2])
deer_merriang_total <- nrow(aggregate(Deer_scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Merriang",],
                                     FUN = function(x) {mean(x > 0)}))
wallabies_merriang_count <- sum(aggregate(Macropod_scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Merriang",],
                                     FUN = function(x) {mean(x > 0)})[,2])
wallabies_merriang_total <- nrow(aggregate(Macropod_scats ~ Plot_number_specific,
                                      data = damage[damage$Location == "Merriang",],
                                      FUN = function(x) {mean(x > 0)}))

deer_bright_count <- sum(aggregate(Deer_scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Bright",],
                                     FUN = function(x) {mean(x > 0)})[,2])
deer_bright_total <- nrow(aggregate(Deer_scats ~ Plot_number_specific,
                                      data = damage[damage$Location == "Bright",],
                                      FUN = function(x) {mean(x > 0)}))
wallabies_bright_count <- sum(aggregate(Macropod_scats ~ Plot_number_specific,
                                          data = damage[damage$Location == "Bright",],
                                          FUN = function(x) {mean(x > 0)})[,2])
wallabies_bright_total <- nrow(aggregate(Macropod_scats ~ Plot_number_specific,
                                           data = damage[damage$Location == "Bright",],
                                           FUN = function(x) {mean(x > 0)}))

damage$scats <- damage$Deer_scats + damage$Macropod_scats
two_merriang_count <- sum(aggregate(scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Merriang" & damage$Age == 2,],
                                     FUN = function(x) {mean(x > 0)})[,2])
two_merriang_total <- nrow(aggregate(scats ~ Plot_number_specific,
                                    data = damage[damage$Location == "Merriang" & damage$Age == 2,],
                                    FUN = function(x) {mean(x > 0)}))

four_merriang_count <- sum(aggregate(scats ~ Plot_number_specific,
                                    data = damage[damage$Location == "Merriang" & damage$Age == 4,],
                                    FUN = function(x) {mean(x > 0)})[,2])
four_merriang_total <- nrow(aggregate(scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Merriang" & damage$Age == 4,],
                                     FUN = function(x) {mean(x > 0)}))

two_bright_count <- sum(aggregate(scats ~ Plot_number_specific,
                                    data = damage[damage$Location == "Bright" & damage$Age == 2,],
                                    FUN = function(x) {mean(x > 0)})[,2])
two_bright_total <- nrow(aggregate(scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Bright" & damage$Age == 2,],
                                     FUN = function(x) {mean(x > 0)}))

four_bright_count <- sum(aggregate(scats ~ Plot_number_specific,
                                     data = damage[damage$Location == "Bright" & damage$Age == 4,],
                                     FUN = function(x) {mean(x > 0)})[,2])
four_bright_total <- nrow(aggregate(scats ~ Plot_number_specific,
                                      data = damage[damage$Location == "Bright" & damage$Age == 4,],
                                      FUN = function(x) {mean(x > 0)}))

# perform Fisher's Exact test for whether the proportion of plots with
# deer/wallaby scats differs between the locations
cont_table <- matrix(c(deer_merriang_count, wallabies_merriang_count, deer_merriang_total, wallabies_merriang_total), ncol = 2, byrow = TRUE)
fisher.test(cont_table)
cont_table <- matrix(c(deer_bright_count, wallabies_bright_count, deer_bright_total, wallabies_bright_total), ncol = 2, byrow = TRUE)
fisher.test(cont_table)

cont_table <- matrix(c(two_merriang_count, four_merriang_count, two_merriang_total, four_merriang_total), ncol = 2, byrow = TRUE)
fisher.test(cont_table)
cont_table <- matrix(c(two_bright_count, four_bright_count, two_bright_total, four_bright_total), ncol = 2, byrow = TRUE)
fisher.test(cont_table)

# plot damage levels based on whether scats are present
pdf(paste("figures_less_intensive/scats.pdf", sep = ""), height = 7.5, width = 8)
par("mfrow" = c(2,2))
par("xpd" = FALSE)
par("mar" = c(5.6, 5.1, 6, 2.1))
P <- signif(wilcox.test(damage_below_plot ~ deer_scats_plot)$p.value, 2)
vioplot(damage_below_plot ~ deer_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2)
mtext(paste("Deer\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (<1m)",
      side = 2, line = 3.5, cex = 1.8)
P <- signif(wilcox.test(damage_below_plot ~ wallaby_scats_plot)$p.value, 2)
vioplot(damage_below_plot ~ wallaby_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2,
        cex.main = 1.8)
mtext(paste("Macropods\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (<1m)",
      side = 2, line = 3.5, cex = 1.8)
P <- signif(wilcox.test(damage_above_plot ~ deer_scats_plot)$p.value, 2)
vioplot(damage_above_plot ~ deer_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2,
        cex.main = 1.8)
mtext(paste("Deer\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (>1m)",
      side = 2, line = 3.5, cex = 1.8)
P <- signif(wilcox.test(damage_above_plot ~ wallaby_scats_plot)$p.value, 2)
vioplot(damage_above_plot ~ wallaby_scats_plot, col = "lightblue",
        xlab = "", ylab = "",
        las = 1,
        cex.axis = 1.8, cex.names = 2,
        cex.main = 1.8)
mtext(paste("Macropods\n(P = ", P, ")", sep = ""), side = 3, cex = 1.8, line = 2)
mtext(text = "% damaged (>1m)",
      side = 2, line = 3.5, cex = 1.8)
dev.off()

# are tree height and age correlated?
stripchart(Tree_height ~ Age, data = damage,
           pch = 1, col = "lightblue",
           vertical = TRUE,
           method = "jitter",
           xlab = "Age", ylab = "Tree height")
cor.test(damage$Tree_height, damage$Age)

# yes, but I will keep both in the model: can be an important 
# covariate when evaluating the effect of age!
cor.test(damage$Tree_height, damage$damage_count)
cor.test(damage$Age, damage$damage_count)

# types of damage as a function of age
damage_only_complete <- damage[, c("Age", "damage_below", "damage_above", 
                                   "Browsing_damage", "Leader_damage")]
damage_only_complete <- damage_only_complete[complete.cases(damage_only_complete),]
below <- aggregate(damage_below ~ Age, data = damage_only_complete,
                   FUN = function(x) {mean(x > 0)})[,2]
above <- aggregate(damage_above ~ Age, data = damage_only_complete,
                   FUN = function(x) {mean(x > 0)})[,2]
leader <- aggregate(Leader_damage ~ Age, data = damage_only_complete,
                    FUN = function(x) {mean(x > 0)})[,2]
browsing <- aggregate(Browsing_damage ~ Age, data = damage_only_complete,
                    FUN = function(x) {mean(x > 0)})[,2]
# always more damage in older trees
# don't know if it's because the animals prefer them
# or if the older trees have had more time to be damaged
# browsing/leader damage are less present than stem damage for the younger trees
pdf("figures_less_intensive/type_of_damage.pdf",
    height = 4, width = 5)
all_data <- as.matrix(cbind(below, above, leader, browsing))
barplot(all_data, beside = TRUE,
        xlab = "Type of damage",
        names.arg = c("Stem <1m", "Stem >1m",
                      "Leader", "Browsing"),
        ylab = "Proportion of damaged trees",
        col = c("lightblue", "steelblue3"))
legend("top", legend = c("2 yrs", "4 yrs"),
       fill = c("lightblue", "steelblue3"))
dev.off()

damage_only_complete$damage_below_binary <- damage_only_complete$damage_below > 0
damage_only_complete$damage_above_binary <- damage_only_complete$damage_above > 0

# leader damage linked to damage below
# browsing damage linked to damage above
# leader and browsing damage correlate
damage_types <- c("Browsing_damage", "Leader_damage",
                  "damage_below_binary",
                  "damage_above_binary")
pretty_types <- c("browsing damage",
                  "leader damage",
                  "stem damage <1m",
                  "stem damage >1m")
for (age in c("2","4")) {
  for (predictor in damage_types) {
    for (predictor2 in damage_types) {
      if (predictor != predictor2) {
        curr_data <- damage_only_complete[damage_only_complete$Age == age,]
        agg <- aggregate(curr_data[, predictor] ~ curr_data[, predictor2],
                         FUN = mean)
        agg_counts <- aggregate(curr_data[, predictor] ~ curr_data[, predictor2],
                         FUN = sum)
        colnames(agg) <- c(predictor2, predictor)
        print(agg)
        perc <- round(agg[,2] * 100, 1)
        count <- agg_counts[,2]
        pdf(paste("figures_less_intensive/", predictor2, "_", predictor, "_Age", age, ".pdf", sep = ""),
            height = 4, width = 4)
        par("xpd" = TRUE)
        bp <- barplot(perc, col = "steelblue",
                      ylab = paste("% trees with", pretty_types[which(damage_types == predictor)]),
                      xlab = paste("Presence of", pretty_types[which(damage_types == predictor2)]),
                      las = 1,
                      names.arg = c("Absent", "Present"),
                      ylim = c(0,100),
                      main = paste("Age: ", age, " yrs", sep = ""))
        text(paste(perc, "% (n=", count,")", sep = ""),
             x = c(0.7, 1.9),
             y = perc + 10)
        dev.off()
      }
    }
  }
}

# tree height and DBH correlate strongly
plot(Tree_height ~ DBH, data = damage, ylab = "Tree height",
     col = "lightblue")
######

# Prepare data 2
######
# only keep the columns that are necessary for analysis
damage <- damage[, to_keep]
dim(damage)
# we lose some because of damage and because of tree height
damage <- damage[complete.cases(damage),]
damage$aspect_categorical <- droplevels(damage$aspect_categorical)
dim(damage)

# add in GPS coordinates
damage$Latitude <- rep(0, dim(damage)[1])
damage$Longitude <- rep(0, dim(damage)[1])
gps <- read.csv("LessIntensiveCoordinatesGPS.tsv", sep = "\t")
plot_number <- c()
for (plot in gps$Point_Name) {
  split <- str_split(plot, "P")[[1]]
  plot_number <- c(plot_number, split[length(split)])
}
gps$Point_name <- plot_number
gps$Plot_number_specific <- paste(gps$Coupe, gps$Point_name, sep = "_")
mean(damage$Plot_number_specific %in% gps$Plot_number_specific)
for (obs in 1:nrow(damage)) {
  plot_name <- damage[obs, "Plot_number_specific"]
  coupe <- damage[obs, "Coupe"]
  coords <- gps[gps$Plot_number_specific == plot_name & gps$Coupe == coupe,]
  damage[obs, "Longitude"] <- coords$Longitude
  damage[obs, "Latitude"] <- coords$Latitude
}
damage$Plot_number <- damage$Plot_number_specific
######

# Model
######

# overview of distribution for final data set
pdf(paste("figures_less_intensive/distribution_of_damage.pdf", sep = ""), width = 6, height = 6)
par(mfrow = c(2,2))
margins <- par("mar")
par("mar" = c(5.1, 5.6, 4.1, 2.1))
make_distribution_plot(damage[damage$Age == "2", "damage_below"],
                       title = "<1m (2 yrs)", to_add = 50)
damage_below2 <- damage[damage$Age == "2", "damage_below"]
mean(damage_below2[damage_below2 > 0])
make_distribution_plot(damage[damage$Age == "2", "damage_above"],
                       title = ">1m (2 yrs)", to_add = 80)
make_distribution_plot(damage[damage$Age == "4", "damage_below"],
                       title = "<1m (4 yrs)", to_add = 50)
damage_below4 <- damage[damage$Age == "4", "damage_below"]
mean(damage_below4[damage_below4 > 0])
make_distribution_plot(damage[damage$Age == "4", "damage_above"],
                       title = ">1m (4 yrs)", to_add = 80)
dev.off()

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
curr_df <- prepare_df_for_sac("all", damage, parameter_narrow = "Coupe",
                              parameter_large = "Location", parameter_narrowest = "Plot_number")
neighbour_mats <- list()
for (k in c(2,3,4)) {
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

# fit model
M5_below <- brm(bf(damage_below|trials(4) ~ scale(Distance_water_fine)*Location + scale(Distance_water_fine)*Age + scale(Density_road)*Location + scale(Density_road)*Age + edge_simple*Age + scale(Tree_height)*Location + scale(Tree_height)*Age + scale(Slope)*Location + scale(Slope)*Age + aspect_categorical*Location + aspect_categorical*Age + scale(Elevation)*Location + scale(Elevation)*Location + Age*Location + car(neighbour_mat, gr = Coupe) + (1|Plot_number), zi ~ scale(Distance_water_fine)*Location + scale(Distance_water_fine)*Age + scale(Density_road)*Location + scale(Density_road)*Age + edge_simple*Age + scale(Tree_height)*Location + scale(Tree_height)*Age + scale(Slope)*Location + scale(Slope)*Age + aspect_categorical*Location + aspect_categorical*Age + scale(Elevation)*Location + scale(Elevation)*Location + Age*Location + car(neighbour_mat, gr = Coupe) + (1|Plot_number)),
                             family = zero_inflated_binomial(),
                             data = damage,
                             backend = "cmdstanr",
                             cores = 4,
                             save_pars = save_pars(all = TRUE),
                             seed = 5,
                             prior = c(prior_b_weak, prior_b_zi_weak,
                                       prior_intercept, prior_intercept_zi,
                                       prior_sdcar,
                                       prior_sd, prior_sd_zi), 
                             data2 = list("neighbour_mat" = neighbour_mats[[3]]),
                             iter = 10000, warmup = 1000,
                             control = list("adapt_delta" = 0.99))

# use the medium-scale metric for distance to water
M5_below_medium <- brm(bf(damage_below|trials(4) ~ scale(Distance_water_medium)*Location + scale(Distance_water_medium)*Age + scale(Density_road)*Location + scale(Density_road)*Age + edge_simple*Age + scale(Tree_height)*Location + scale(Tree_height)*Age + scale(Slope)*Location + scale(Slope)*Age + aspect_categorical*Location + aspect_categorical*Age + scale(Elevation)*Location + scale(Elevation)*Location + Age*Location + car(neighbour_mat, gr = Coupe) + (1|Plot_number), zi ~ scale(Distance_water_medium)*Location + scale(Distance_water_medium)*Age + scale(Density_road)*Location + scale(Density_road)*Age + edge_simple*Age + scale(Tree_height)*Location + scale(Tree_height)*Age + scale(Slope)*Location + scale(Slope)*Age + aspect_categorical*Location + aspect_categorical*Age + scale(Elevation)*Location + scale(Elevation)*Location + Age*Location + car(neighbour_mat, gr = Coupe) + (1|Plot_number)),
                family = zero_inflated_binomial(),
                data = damage,
                backend = "cmdstanr",
                cores = 4,
                save_pars = save_pars(all = TRUE),
                seed = 5,
                prior = c(prior_b_weak, prior_b_zi_weak,
                          prior_intercept, prior_intercept_zi,
                          prior_sdcar,
                          prior_sd, prior_sd_zi), 
                data2 = list("neighbour_mat" = neighbour_mats[[3]]),
                iter = 10000, warmup = 1000,
                control = list("adapt_delta" = 0.99))

# visualise the model predictions for the effect of different predictors
pdf("figures_less_intensive/road_density_water_distance.pdf", width = 18, height = 4)
par(mfrow = c(1,3))
make_epred_plot(M5_below, damage, c(0.05, 0.25), "Density_road", "Road density", original_vector = NA, title = "<1m", by_location_age = TRUE, plot_to_file = FALSE, show_legend = TRUE, enlarge = TRUE, panel_label = "A")
make_epred_plot(M5_below, damage, c(0, 450), "Distance_water_fine", "Distance to water (fine)", original_vector = NA, title = "", by_location_age = TRUE, plot_to_file = FALSE, show_legend = FALSE, enlarge = TRUE, panel_label = "B")
make_epred_plot(M5_below_medium, damage, c(0, 450), "Distance_water_medium", "Distance to water (medium)", original_vector = NA, title = "", by_location_age = TRUE, plot_to_file = FALSE, show_legend = FALSE, enlarge = TRUE, panel_label = "C")
dev.off()

pdf("figures_less_intensive/aspect.pdf", width = 10, height = 4)
par(mfrow = c(1,2))
make_epred_plot(M5_below, damage, c("W", "E"), "aspect_categorical", "Aspect", original_vector = NA, title = "", by_location_age = TRUE, plot_to_file = FALSE, show_legend = TRUE, enlarge = FALSE, panel_label = "A")
make_epred_plot(M5_below, damage, c("N", "S"), "aspect_categorical", "Aspect", original_vector = NA, title = "", by_location_age = TRUE, plot_to_file = FALSE, show_legend = FALSE, enlarge = FALSE, panel_label = "B")
dev.off()

pdf("figures_less_intensive/elevation_slope.pdf", width = 10, height = 4)
par(mfrow = c(1,2))
make_epred_plot(M5_below, damage, c(250, 500), "Elevation", "Elevation", original_vector = NA, title = "", by_location_age = TRUE, plot_to_file = FALSE, show_legend = TRUE, enlarge = FALSE, panel_label = "A")
make_epred_plot(M5_below, damage, c(0, 35), "Slope", "Slope", original_vector = NA, title = "", by_location_age = TRUE, plot_to_file = FALSE, show_legend = FALSE, enlarge = FALSE, panel_label = "B")
dev.off()

# high- vs low-risk scenario
new_data <- damage[1:2,]
new_data$edge_simple <- "Mature vegetation"
new_data$aspect_categorical <- "W"
new_data$Slope <- median(damage$Slope)
new_data$Distance_water_fine <- median(damage$Distance_water_fine)
new_data$Age <- "4"
new_data$Location <- "Merriang"
new_data$Plot_number_specific <- "3_2"
new_data$Tree_height <- c(5.5, 5.5)
new_data$Density_road <- c(0.25, 0.05)
new_data$Elevation <- c(500, 250)
set.seed(5)
pred_means <- posterior_predict(M5_below, newdata = new_data, 
                              re_formula = NULL)
apply(pred_means, 2, FUN = mean)
apply(pred_means, 2, FUN = median)
apply(pred_means, 2, FUN = IQR)
apply(pred_means, 2, FUN = function(x) {mean(x == 0)})

HPDinterval(as.mcmc(pred_means[,1]))
diffs <- pred_means[,2] - pred_means[,1]
median(diffs)
HPDinterval(as.mcmc(diffs))
mean(diffs > 0)
hist(diffs)
pdf("figures_less_intensive/high_low_risk.pdf", width = 5, height = 3)
par(mfrow = c(1,1))
blue_trans <- rgb(173, 216, 230, alpha = 80, 
                  maxColorValue = 255)
desert_trans <- rgb(250, 213, 165, alpha = 80, 
                    maxColorValue = 255)
vioplot(pred_means[,1], pred_means[,2], damage$damage_below[damage$Location == "Merriang"], damage$damage_below[damage$Location == "Bright"], names = c("Low-risk", "High-risk", "Merriang", "Bright"),
        col = c("grey77", "grey35", blue_trans, desert_trans),
        ylab = "# of quadrants damaged\n",
        yaxt = "n", cex.lab = 1.5)
axis(side = 2, las = 2, at = seq(0,4),
     labels = seq(0,4), cex.axis = 1)
dev.off()

# posterior predictions
pdf("figures_less_intensive/PPD.pdf", height = 8, width = 8.5)
margins <- par("mar")
par("mar" = c(10.1, 6.5, 4.1, 2.1))
par(mfrow = c(2,2))
pp_check_custom(M5_below, damage$damage_below, "", "A")
cuts <- seq(0, 12)
labels <- paste(seq(0,11), seq(1,12), sep = "-")
ppd_numerical_predictor(M5_below, damage, "B", "Tree_height", "damage_below", 
                        "Tree height (m)", "", cuts = cuts, labels = labels,
                        show_legend = FALSE)
cuts <- seq(0.05, 0.20, by = 0.025)
labels <- paste(seq(0.05,0.175, by = 0.025), seq(0.075,0.20, by = 0.025), sep = "-")
ppd_numerical_predictor(M5_below, damage, "C", "Density_road", "damage_below", 
                        "Road density", "", cuts = cuts, labels = labels,
                        show_legend = FALSE)
cuts <- seq(0, 400, by = 100)
labels <- paste(seq(0,300, by = 100), seq(100,400, by = 100), sep = "-")
ppd_numerical_predictor(M5_below, damage, "D", "Distance_water_fine", 
                        "damage_below", "Distance to water (m)", "", 
                        cuts = cuts, labels = labels,
                        show_legend = FALSE)
dev.off()
par("mar" = margins)
######