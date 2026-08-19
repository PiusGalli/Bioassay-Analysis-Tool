read_setup_plate <- function(setup.path, plate.path){
  plate <- read.csv(plate.path, check.names = FALSE)
  setup <- read.csv(setup.path, check.names = FALSE)
  return(list(plate = plate, setup = setup))
}

read.DataSetup <- function(setup, plate) {
  data <- checkInput(plate)
  
  out <- ReplaceSampleNames(setup)
  setup <- out[[1]]
  lookup <- out[[2]]
  
  colnames(data)[1] <- "Well"
  
  checkTemp(data)
  
  data <- data[-2,] #Remove Temperature
  
  setup.list <- getSetup(setup)
  setup <- setup.list[[1]]
  c.setup <- setup.list[[2]]
  
  merged.data <- merge(setup, data, by = "Well")
  
  out <- list(merged.data, c.setup, lookup)
  return(out) 
}

ReplaceSampleNames <- function(setup) {
  original     <- as.character(setup$Sample)
  unique_names <- unique(original)
  
  width        <- max(2, nchar(length(unique_names)))
  placeholders <- paste0("Sample.", formatC(seq_along(unique_names),
                                            width = width, flag = "0"))
  
  setup$Sample <- placeholders[match(original, unique_names)]
  
  list(setup  = setup,
       lookup = setNames(unique_names, placeholders))
}

checkInput <- function(data) {
  if ("Device: infinite 200Pro" %in% data[, 1]) {
    idx <- which(data[, 1] == "Cycle Nr.")
    colnames(data) <- data[idx,]
    data <- data[idx+1:(nrow(data)-idx), ]
  }
  
  if ("End Time:" %in% data[, 1]) {
    data <- data[1:(nrow(data) - 5), ]
  }
  
  if (any(is.na(data))) {
    na <- which(is.na(data[3, ]))[1]
    data <- data[, 1:na-1]
  }
  
  return(data)
}

checkTemp <- function(data) {
  temp <- as.vector(t(data[2, 3:ncol(data)]))
  temp <- as.numeric(temp)
  
  if (any(temp > 31) || any(temp < 29)) {
    index.temp <- which(temp < 29 | temp > 31)
    temp.wrong <- temp[index.temp]
    return(warning(paste0("Index: ", index.temp, ". Temp: ", temp.wrong, collapse = " !!! ")))
  }
}

getSetup <- function(setup) {
  well.index <- which(str_detect(colnames(setup), "Well"))
  rest.index <- which(!str_detect(colnames(setup), "Well"))
  
  order.index <- append(well.index, rest.index)
  setup <- setup[, order.index]
  
  c.addition <- c("Time [s]")
  NA.addition <- rep(NA, ncol(setup) - 1)
  c.addition <- append(c.addition, NA.addition)
  
  setup <- rbind(c.addition, setup)
  setup.names <- colnames(setup)
  c.setup <- setup.names[!str_detect(setup.names, "Well")]
  
  setup.list <- list(setup, c.setup)
  
  return(setup.list)
}

#Function which reads data and checks the input file if it is ready for processing.

calculate.MeanSD <- function(out.read.DataSetup, time.unit = "min", df.mode = "Standard", starting.od = 0.1) {
  data <- out.read.DataSetup[[1]]
  setup <- out.read.DataSetup[[2]]
  lookup <- out.read.DataSetup[[3]]
  
  split.data <- SplitData(data, setup, lookup)
  
  length.indices <- length(which(colnames(data) %in% setup))
  
  time <- GetTime(data, length.indices, time.unit)
  blank <- GetBlank(split.data)
  
  df <- ProcessList(split.data, time, blank, df.mode, time.unit) #Looses Columns
  
  if (df.mode == "Standard") {
    df <- StandardFirstMeasurmentCorrection(df, time.unit, starting.od)
  } else {
    df <- ReplicateFirstMeasurmentCorrection(df, time.unit, starting.od)
  }
  
  return(df)
}

SplitData <- function(data, setup, lookup) {
  split.indices <- which(colnames(data) %in% setup)
  
  tmp.list <- list()
  for (i in 1:length(split.indices)) {
    tmp.list <- append(tmp.list, list(data[, split.indices[i]]))
  }
  
  split.data <- split(data[(length(split.indices) + 2):length(data)], tmp.list)
  
  split.data <- ifelse(lapply(split.data, nrow) != 0, split.data, NA)
  split.data <- split.data[!is.na(split.data)]
  
  split.data <- RetrieveSampleNames(split.data, lookup)
  
  return(split.data)
}

RetrieveSampleNames <- function(split.data, lookup) {
  New.list <- list()
  for (i in 1:length(lookup)) {
    old <- names(lookup)[i]
    new <- lookup[[i]]
    tmp.list <- split.data[str_detect(names(split.data), fixed(old))]
    names(tmp.list) <- gsub(old, new, names(tmp.list), fixed = TRUE)
    
    New.list <- append(New.list, tmp.list)
  }
  
  return(New.list)
}

GetTime <- function(data, length.indices, time.unit) {
  transposed.data <- as.data.frame(t(data))
  colnames(transposed.data) <- transposed.data[1, ]
  time <- as.numeric(transposed.data$`Time [s]`[-(1:(length.indices + 1))])
  
  if (time.unit == "min") {
    time <- time / 60
  } else if (time.unit == "h") {
    time <- time / 60 / 60
  }
  
  return(time)
}

GetBlank <- function(split.data) {
  blank.list <- split.data[str_detect(names(split.data), "blank|Blank")]
  blank.list <- blank.list[[1]]
  
  blank.list <- lapply(blank.list, as.numeric)
  blank.list <- sapply(blank.list, mean)
  
  blank <- mean(blank.list)
  
  return(blank)
}

ProcessList <- function(split.data, time, blank, df.mode, time.unit) {
  df <- data.frame(Time = time)
  split.data <- split.data[!str_detect(names(split.data), "blank|Blank")]
  list.names <- names(split.data)
  fun.getOD <- function(x) {(x - blank) * 10}
  
  for (i in 1:length(split.data)) {
    tmp.c <- split.data[[i]]
    name <- list.names[i]
    
    if (df.mode == "Standard") {
      tmp.df.mean <- ProcessListChunk(tmp.c, fun.getOD, time, blank, name, mode = "mean")
      tmp.df.sd <- ProcessListChunk(tmp.c, fun.getOD, time, blank, name, mode = "sd")
      
      df <- inner_join(df, tmp.df.mean, by = "Time")
      df <- inner_join(df, tmp.df.sd, by = "Time")
    } else if (df.mode == "Replicate") {
      tmp.df.rep <-  ProcessListChunk(tmp.c, fun.getOD, time, blank, name, mode = "replicates")
      
      df <- inner_join(df, tmp.df.rep, by = "Time")
    }
  }
  
  
  
  return(df)
}

StandardFirstMeasurmentCorrection <- function(df, time.unit, starting.od) {
  df.2to7 <- df[2:7, !(names(df) %in% "Time")]
  df.rest <- df[8:nrow(df),]
  
  df.means <- colMeans(df.2to7)
  df.offset <- df.means[seq(1, length(df.means), by = 2)]
  df.offset <- starting.od - df.offset
  
  if (time.unit == "min") {
    df.means <- c(14.5, df.means)
  } else if (time.unit == "h") {
    df.means <- c(0.242, df.means)
  } else {
    df.means <- c(870, df.means)
  }
  
  df.new <- rbind(df.means, df.rest)
  df.new[, seq(2, ncol(df), by = 2)] <- sweep(df.new[, seq(2, ncol(df), by = 2)], 2, df.offset, FUN = "+")
  
  print(paste("The Offsets for", names(df.offset), "is:", round(df.offset, 2)))
  
  return(df.new)
  
}

ReplicateFirstMeasurmentCorrection <- function(df, time.unit, starting.od) {
  df.2to7 <- df[2:7, !(names(df) %in% "Time")]
  df.rest <- df[8:nrow(df),]
  
  df.means <- colMeans(df.2to7)
  df.offset <- df.means
  df.offset <- starting.od - df.offset
  
  if (time.unit == "min") {
    df.means <- c(14.5, df.means)
  } else if (time.unit == "h") {
    df.means <- c(0.242, df.means)
  } else {
    df.means <- c(870, df.means)
  }
  
  df.new <- rbind(df.means, df.rest)
  df.new[, seq(2, ncol(df))] <- sweep(df.new[, seq(2, ncol(df))], 2, df.offset, FUN = "+")
  
  print(paste("The Offsets for", names(df.offset), "is:", round(df.offset, 2)))
  
  return(df.new)
  
}

ProcessListChunk <- function(tmp.c, fun.getOD, time, blank, name, mode = "mean") {
  tmp.c <- lapply(tmp.c, as.numeric)
  tmp.c <- lapply(tmp.c, fun.getOD)
  
  
  
  if (mode == "mean") {
    tmp.c <- lapply(tmp.c, mean)
  } else if (mode == "sd") {
    tmp.c <- lapply(tmp.c, sd)
    name <- paste(name, "sd", sep = "_")
  } else if (mode == "replicates") {
    tmp.rep.df <- as.data.frame(tmp.c)
    tmp.rep.df <- as.data.frame(t(tmp.rep.df))
    
    number_reps <- ncol(tmp.rep.df)
    
    rep <- paste0(".Rep", seq(1, number_reps))
    names <- paste(name, rep, sep = "")
    
    colnames(tmp.rep.df) <- names
    
    tmp.rep.df$Time <- time
    
    return(tmp.rep.df)
  }
  
  tmp.df <- as.data.frame(tmp.c)
  tmp.df <- as.data.frame(t(tmp.df))
  
  colnames(tmp.df) <- name
  
  tmp.df$Time <- time
  
  return(tmp.df)
}

################################ Fit functions

fit.logistic <- function(data, time) {
  d <- data.frame(y = data, t = time)
  
  k_init <- max(data) # carrying capacity is near the max
  n0_init <- min(data[data > 0]) # initial population size is near the min
  
  # Initial estimate for r
  glm_mod <- stats::glm(y / k_init ~ t,
                        family = stats::quasibinomial("logit"),
                        data = d)
  
  r_init <- stats::coef(glm_mod)[[2]] # slope which should only be positive
  if (r_init <= 0) {
    r_init <- 0.01
  }
  
  nls_mod <- nls(y ~ k / (1 + ( (k - n0) / n0) * exp(-r * t)),
                 start = list(k = k_init,
                              r = r_init,
                              n0 = n0_init),
                 control = list(maxiter = 500),
                 data = d)
  
  return(nls_mod)
}

safe.fit.logistic <- function(...){
  purrr:safely(fit.logistic, quiet = FALSE)
}

fit_auc_ci <- function(df, idx = 9, cut.at = 3000) {
  df <- na.omit(df)
  
  colnames_only <- names(df)[-1]
  base_of_col   <- sub("\\.Rep\\d+$", "", colnames_only)
  base_names    <- unique(base_of_col)
  
  auc_summary <- data.frame()
  d_fit <- data.frame(Time = df$Time)
  mean_sd_df <- data.frame()
  growth_params <- data.frame()
  auc_vals_df <- data.frame()
  
  for (name in base_names) {
    cols <- colnames_only[base_of_col == name]
    auc_vals <- numeric(length(cols))
    
    for (i in seq_along(cols)) {
      d_loop <- df[, c("Time", cols[i])]
      gc_fit <- SummarizeGrowth(data_t = d_loop$Time, data_n = d_loop[[2]])
      auc_vals[i] <- gc_fit$vals$auc_l
    }
    
    # Calculate mean + SD curve
    mean_vals <- rowMeans(df[cols], na.rm = TRUE)
    sd_vals   <- apply(df[cols], 1, sd, na.rm = TRUE)
    
    # Fit to mean curve to extract summary params
    mean_fit_successful <- tryCatch({
      mean_fit <- SummarizeGrowth(data_t = df$Time, data_n = mean_vals)
      
      # Store growth parameters based on mean curve
      growth_params <- rbind(growth_params, data.frame(
        Sample = name,
        k = mean_fit$vals$k,
        n0 = mean_fit$vals$n0,
        r = mean_fit$vals$r,
        t_mid = mean_fit$vals$t_mid,
        t_gen = mean_fit$vals$t_gen,
        auc_l = mean_fit$vals$auc_l,
        auc_e = mean_fit$vals$auc_e,
        sigma = mean_fit$vals$sigma
      ))
      
      # Predict curve for plotting
      pred_vals <- predict(mean_fit$model)
      d_fit_loop <- setNames(data.frame(pred_vals), name)
      d_fit <- cbind(d_fit, d_fit_loop)
      
      TRUE
    }, error = function(e) {
      message(paste("Mean fit failed for", name, ":", e$message))
      FALSE
    })
    
    auc_vals_df <- rbind(auc_vals_df, data.frame(
      Sample = name,
      Replicate = seq_along(auc_vals),
      auc = auc_vals
    ))
    
    auc_mean <- mean(auc_vals, na.rm = TRUE)
    auc_sd <- sd(auc_vals, na.rm = TRUE)
    auc_se <- auc_sd / sqrt(length(auc_vals))
    ci <- qt(0.975, df = length(auc_vals) - 1) * auc_se
    ci_lower <- auc_mean - ci
    ci_upper <- auc_mean + ci
    
    auc_summary <- rbind(auc_summary, data.frame(
      Sample = name,
      auc_l_mean = auc_mean,
      auc_l_sd = auc_sd,
      auc_l_se = auc_se,
      auc_l_ci_lower = ci_lower,
      auc_l_ci_upper = ci_upper
    ))
    
    mean_vals <- rowMeans(df[cols], na.rm = TRUE)
    sd_vals   <- apply(df[cols], 1, sd, na.rm = TRUE)
    
    mean_sd_df <- rbind(mean_sd_df,
                        data.frame(Time = df$Time,
                                   Sample = name,
                                   mean = mean_vals,
                                   sd = sd_vals,
                                   lower = mean_vals - sd_vals,
                                   upper = mean_vals + sd_vals))
    
    d_fit_successful <- tryCatch({
      mean_fit <- SummarizeGrowth(data_t = df$Time, data_n = mean_vals)
      pred_vals <- predict(mean_fit$model)
      d_fit_loop <- setNames(data.frame(pred_vals), name)
      d_fit <- cbind(d_fit, d_fit_loop)
      TRUE
    }, error = function(e) {
      message(paste("Failed to fit for", name, ":", e$message))
      FALSE
    })
  }
  
  d_fit_long <- if (ncol(d_fit) > 1) {
    d_fit %>%
      pivot_longer(-Time, names_to = "Sample", values_to = "Value")
  } else {
    data.frame(Time = numeric(0), Sample = character(0), Value = numeric(0))
  }
  
  mean_sd_df <- mean_sd_df %>%
    group_by(Sample) %>%
    filter(row_number() %% idx == 1) %>%
    ungroup()
  
  d_fit_long <- subset(d_fit_long, Time <= cut.at)
  mean_sd_df <- subset(mean_sd_df, Time <= cut.at)
  
  legend_labels <- auc_summary %>%
    mutate(
      legend = paste0(Sample, ": ",
                      round(auc_l_mean, 0), " ± ",
                      round(auc_l_sd, 0))
    ) %>%
    select(Sample, legend, auc_l_mean)
  
  ordered_levels <- legend_labels %>%
    arrange(desc(auc_l_mean)) %>%
    pull(legend)
  
  mean_sd_df <- mean_sd_df %>%
    left_join(legend_labels, by = "Sample") %>%
    mutate(legend = factor(legend, levels = ordered_levels))
  
  d_fit_long <- d_fit_long %>%
    left_join(legend_labels, by = "Sample") %>%
    mutate(legend = factor(legend, levels = ordered_levels))
  
  auc_summary <- auc_summary %>%
    left_join(legend_labels %>% select(Sample, legend), by = "Sample")
  
  anova_result <- aov(auc ~ Sample, data = auc_vals_df)
  anova_summary <- summary(anova_result)
  
  tukey_result <- TukeyHSD(anova_result)
  
  return(list(
    auc_summary = auc_summary,
    d_fit_long = d_fit_long,
    mean_sd_df = mean_sd_df,
    growth_params = growth_params,
    rep_auc_vals = auc_vals_df,
    anova_summary = anova_summary,
    tukey_result = tukey_result
  ))
}

################################ Plot functions
# Plot Functions which either plot all curves alone in one window or plots all curves in one df in one plot

plot.separateGrowthcurves <- function(data, nr.row, nr.col, time.unit = "min", lab.size = 1.3, 
                              plot.point = 4, ylim = NA, pch = 16, cut.min = 2000, cut.h = 30, error = TRUE, log.mode = FALSE) {
  
  if (log.mode == TRUE) {
    log.mode <- "y"
    y.start <- 0.05
  } else {
    log.mode <- ""
    y.start <- 0
  }
  
  if (time.unit == "min") {
    data <- subset(data, data$Time < cut.min)
  } else {
    data <- subset(data, data$Time < cut.h)
  }

  if (is.na(ylim)) {
    max.OD <- max(data[-1])
  } else {
    max.OD <- ylim
  }
  
  idx <- seq(1, nrow(data), plot.point)
  oldparams <- par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))
  par(mfrow = c(nr.row, nr.col), oma = c(2, 2, 0, 0))
  
  for (i in 2:ncol(data)) {
    if (!str_detect(colnames(data)[i], ".sd")) {
      tmp.fit <- safe.fit.logistic(data[ , i], data$Time)
      r <- 0
      sd <- 0
      
      if (is.null(tmp.fit$error)) {  
        r <- summary(tmp.fit$result)$coefficients[2]
        sd <- summary(tmp.fit$result)$coefficients[2, 2]
      }
      
      tryCatch(
        plot(data$Time[idx], data[idx, i], 
             pch = 16,
             xlab = "", ylab = "", 
             main = paste0(colnames(data)[i], "\n ", sprintf("r = (%.4f ± %.5f)/min", r, sd)),
             ylim = c(y.start, max.OD), log = log.mode),
        
        error = function(e) {
          par(oldparams)
          stop(e)
        }
      )
      
      if (error) {
        data.sd <- data[idx, i + 1]
        arrows(x0 = data$Time[idx], y0 = data[idx, i] - data.sd, x1 = data$Time[idx], y1 = data[idx, i] + data.sd, 
               code = 3, angle = 90, length = 0.05)
      }
      
      if (is.null(tmp.fit$error)) {
        lines(data$Time, predict(tmp.fit[[1]]), col = "red")
        #legend(x = 0, y = max.OD + 0.2, legend = c(sprintf("(%.4f ± %.5f)/min", r, sd)))
      }
    }
  }
  
  mtext(paste0("Time [", time.unit, "]"), side = 1, line = 0, outer = TRUE, cex = lab.size)
  mtext(expression("OD"[600]), side = 2, line = 0, outer = TRUE, cex= lab.size, las = 0)
  
  par(oldparams)
}

# All, split in motifs (ie. different concentrations but several compounds and vice-versa)
# Cut data early

DFtoDFLong <- function(data, plot.point) {
  # Create both data frames directly by filtering the column names
  data.noSD <- data[seq(1, nrow(data), by = plot.point), grep("_sd$", colnames(data), invert = TRUE)]
  data.SD <- data[seq(1, nrow(data), by = plot.point), grep("_sd$", colnames(data))]
  
  # Rename the columns in data.SD and add Time directly
  colnames(data.SD) <- sub("_sd$", "", colnames(data.SD))
  data.SD$Time <- data.noSD$Time
  
  long_data.noSD <- data.noSD %>%
    pivot_longer(cols = -Time, names_to = "Variable", values_to = "OD")
  
  long_data.SD <- data.SD %>%
    pivot_longer(cols = -Time, names_to = "Variable", values_to = "sd")
  
  merged_df <- long_data.noSD %>%
    inner_join(long_data.SD, by = c("Time", "Variable"))
  
  return(merged_df)
}

PrepToPlot <- function(data) {
  variables <- unique(data$Variable)
  
  # Fit logistic models for each variable and collect results
  results <- map_dfr(variables, function(variable) {
    # Filter the data for the current variable
    variable_data <- data %>% 
      filter(Variable == variable)
    
    # Ensure there are enough data points for fitting
    if (nrow(variable_data) < 4) {
      return(tibble(Variable = variable, Fitted = NA, error = "Not enough data points"))
    }
    
    # Fit the logistic model using safe.fit.logistic
    fit_result <- safe.fit.logistic(variable_data$OD, variable_data$Time)
    
    # Check if the fitting was successful
    if (!is.null(fit_result$result)) {
      r <- summary(fit_result$result)$coefficients[2]
      sd <- summary(fit_result$result)$coefficients[2, 2]
      #coeffs <- get.Coefficients(fit_result) #c("k", "k_se", "k_p", "r", "r_se", "r_p", "n0", "n0_se", "n0_p", "t_mid", "t_double")
      #print(coeffs)
      
      # Use the fitted model to predict values
      predictions <- data.frame(Time = variable_data$Time,
                                Variable = variable,
                                Fitted = predict(fit_result$result, newdata = variable_data),
                                error = NA,
                                growthrate = r,
                                growthrate_sd = sd)
                                #coeff = coeffs)  # No error here
      
      return(predictions)
    } else {
      return(tibble(Variable = variable, Fitted = NA, error = as.character(fit_result$error)))
    }
  })
  

  # Clean results by filtering out errors
  fitted_results <- results %>% filter(is.na(error))  # Keep rows where there is no error
  error_results <- results %>% filter(!is.na(error))  # Store the error results if needed
  
  # Merging the original data with fitted values
  result <- try({
    final_data <- data %>% 
      left_join(fitted_results, by = c("Time", "Variable"))
  }, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    final_data <- data
  }
  
  return(final_data)
  
}

plot.allGrowthcurves <- function(data, time.unit = "min", plot.point = 6, cut.min = 2000, cut.h = 30, ylim = NULL,
                                 error = TRUE, log.mode = FALSE, main = "", col = NULL, custom.colors = NULL, 
                                 ...) {
  
  if (time.unit == "min") {
    new.data <- subset(data, data$Time < cut.min)
  } else {
    new.data <- subset(data, data$Time < cut.h)
  }
  
  merged_df <- DFtoDFLong(new.data, plot.point)
  final_data <- PrepToPlot(merged_df)
  unique.variable <- unique(final_data$Variable)
  df_legend <- subset(final_data, final_data$Time == final_data$Time[1])
  constructs <- as.vector(df_legend$Variable)
  #constructs <- paste(df_legend$Variable, ": ", round(df_legend$growthrate, 5), " ± ", round(df_legend$growthrate_sd, 5), "/min", sep = "")
  
  len.data <- length(constructs)
  color_map <- c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
  
  # Determine color scheme based on `col` argument
  if (!is.null(custom.colors)) {
    if (length(custom.colors) != len.data) {
      stop(paste("The length of custom.colors must be equal to the number of constructs (", len.data, ").", sep = ""))
    }
    cols <- custom.colors
  } else if (is.null(col)) {
    cols <- viridis(len.data, option = "turbo", ...)
  } else if (col %in% color_map) {
    cols <- viridis(len.data, option = col, ...)
  } else {
    stop('Invalid color option. Choose from: "magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo".')
  }
  
  if ("Fitted" %in% colnames(final_data)) {
    plot <- ggplot(final_data, aes(x = Time)) +
      geom_point(aes(y = OD, color = Variable)) +
      geom_line(aes(y = Fitted, color = Variable)) +
      labs(title = main, x = paste("Time [", time.unit, "]", sep = ""), y = expression("OD"[600]), col = "Construct") +
      scale_color_manual(values = cols) + #, labels = constructs
      theme_classic()
  } else {
    print("Any fit did work, plotting only points.")
    plot <- ggplot(final_data, aes(x = Time)) +
      geom_point(aes(y = OD, color = Variable)) +
      labs(title = main, x = paste("Time [", time.unit, "]", sep = ""), y = expression("OD"[600]), col = "Construct") +
      scale_color_manual(values = cols) + #, labels = constructs
      theme_classic()
  }

  # Conditionally add ylim if provided
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim)
  }
  
  if (error & log.mode == FALSE) {
    plot <- plot + geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd, color = Variable))
  } else if (log.mode & error == FALSE) {
    plot <- plot + scale_y_log10()
  } else if (log.mode & error) {
    plot <- plot + scale_y_log10() + geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd, color = Variable))
  }
  print(plot)
  
  #results <- SplitCoeffs(final_data)
  
  #return(results)
}

SplitCoeffs <- function(data) {
  unique_data <- subset(data, data$Time == data$Time[1])
  
  result_df <- unique_data %>%
    filter(!is.na(coeff)) %>%
    separate(coeff, into = c("k", "k_se", "k_p", "r", "r_se", "r_p", "n0", "n0_se", "n0_p", "t_mid", "t_double"), sep = " ") %>%
    mutate(across(starts_with("k"), as.numeric)) %>%
    select(Variable, k, k_se, k_p, r, r_se, r_p, n0, n0_se, n0_p, t_mid, t_double)
  
  return(result_df)
}

SubsetByColnames <- function(df, string, negate = FALSE) {
  string <- paste(string, collapse = "|")

  if (negate) {
    return(df[!str_detect(colnames(df), string)])
  } else {
    return(df[str_detect(colnames(df), paste(string, "|Time", sep = ""))])
  }
}






