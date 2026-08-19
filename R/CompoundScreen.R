read_plate <- function(plate.path) {
  tryCatch(
    suppressMessages(as.data.frame(read_excel(plate.path))),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("cannot be opened|zip file", msg)) {
        stop(sprintf("File %s cannot be accessed because it is open", plate.path),
             call. = FALSE)
      }
      stop(e) 
    }
  )
}

prepare.compound.map <- function(map.path) {
  map <- read_plate(map.path)
  
  map <- as.data.frame(map)
  map$row.col <- paste(map$Row, map$Column, sep = "/")
  map$Name <- ifelse(is.na(map$Name), "DMSO", map$Name)
  map$CatalogID <- ifelse(is.na(map$CatalogID), "DMSO", map$CatalogID)

  return(map)
}

get_plate <- function(map, plate.paths, barcode = NULL) {
  plate <- read_plate(plate.paths)
  
  file_barcode <- plate[13, c("...5")]        
  if (is.na(file_barcode)) stop("Missing barcode in file: ", plate.paths)
  
  if (!is.null(barcode) && file_barcode != barcode) {
    stop("Barcode mismatch")                   
  }
  barcode <- file_barcode
  
  values <- plate[34:50, ]
  colnames(values) <- values[1,]
  values <- values[2:nrow(values), 2:ncol(values)]
  rownames(values) <- 1:16
  
  if (!barcode %in% map$PlateID) {
    stop(barcode, " not found in plate map. ", plate.paths)
  }
  
  plate_map <- subset(map, PlateID == barcode)
  
  plate_long <- values %>%
    rownames_to_column("Row") %>%
    pivot_longer(
      cols = -Row,
      names_to = "Column",
      values_to = "Value"
    )
  
  plate_long$row.col <- paste(plate_long$Row, plate_long$Column, sep = "/")
  
  final <- merge(plate_map, plate_long[, c("Value", "row.col")], by = "row.col")
  final$Value <- as.numeric(final$Value)
  
  return(final)
}


process.plates <- function(plates, dmso = 2, pos = 23){
  
  data <- do.call(rbind, lapply(plates, function(p) {
    prepare.plates(plate = p, dmso = dmso, pos = pos)
  }))
  
  data[c("row.col", "BatchID", "Supplier", "Type", "Volume", "Quantity")] <- NULL
  
  return(data)
}

prepare.plates <- function(plate, dmso = 2, pos = 23) {
  df <- plate
  
  df$robust_z_dmso <- calculate_z_score(df, dmso)
  df$robust_z_pos <- calculate_z_score(df, pos)
  
  df$plate_zprime <- check_plate(df)
  
  pos_values <- df$Value[df$Column == pos]
  dmso_values <- df$Value[df$Column == dmso]
  
  min.signal <- 500 # Ignore every DMSO control below 500 counts/s to calculate the mean
  
  valid_pos  <- pos_values[pos_values >= min.signal]
  valid_dmso <- dmso_values[dmso_values >= min.signal]
  
  positive_mean <- mean(valid_pos)
  dmso_mean     <- mean(valid_dmso)

  df <- subset(df, df$Type != "negative")
  df <- subset(df, df$Type != "positive")
  
  df$rel.luminescence <- df$Value / dmso_mean
  df$to.positive <- df$Value / positive_mean
  
  return(df)
}

run_plate_qc <- function(matched_data) {
  tapply(matched_data$plate_zprime, matched_data$PlateID, function(x) x[1])
}

read_all_plates <- function(map.paths, plate.paths){
  map <- do.call(rbind, lapply(map.paths, prepare.compound.map))
  plates <- lapply(plate.paths, function(p) get_plate(map, plate.paths = p))
  names(plates) <- sapply(plates, function(x) unique(x$PlateID))
  plates
}

plot.plate.bar <- function(df.comp, df.dmso.pos, df.plate.pos) {
  max.value <- max(c(max(df.comp$Value), max(df.dmso.pos$Value), max(df.plate.pos$Value)))
  
  pos_values <- df.dmso.pos$Value[df.dmso.pos$Type == "positive"]
  dmso_values <- df.dmso.pos$Value[df.dmso.pos$Type == "negative"]
  
  min.signal <- 500
  
  valid_pos  <- pos_values[pos_values >= min.signal]
  valid_dmso <- dmso_values[dmso_values >= min.signal]
  
  mean_pos <- mean(valid_pos)
  mean_dmso     <- mean(valid_dmso)
  
  all_pos <- rbind(df.dmso.pos, df.plate.pos)
  
  ctrl_compound <- df.plate.pos$Name[df.plate.pos$Type == "positive"][1]
  new_compound  <- df.dmso.pos$Name[df.dmso.pos$Type == "positive"][1]
  
  all_pos$binder <- ifelse(all_pos$Column %in% c(1, 24),
                      paste0("Plate control\n(", ctrl_compound, ")"),
                      paste0("New binder\n(", new_compound, ")"))
  all_pos$treatment <- ifelse(all_pos$Column %in% c(1, 2), "DMSO", "Positive compound")
  
  
  
  p1 <- ggplot() +
    geom_point(data = df.comp, aes(x = CatalogID, y = Value), col = "grey50") +
    xlab("Compounds") +
    ylab("Luminescence") +
    ylim(c(0, max.value)) +
    geom_hline(yintercept= mean_dmso, linetype="dashed", color = "red") +
    geom_hline(yintercept= mean_pos, linetype="dashed", color = "red") +
    theme_classic() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  p2 <- ggplot(all_pos, aes(x = binder, y = Value, fill = treatment)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                alpha = 0.6, size = 1) +
    scale_fill_manual(values = c("DMSO" = "grey80", "Positive compound" = "coral1")) +
    ylim(c(0, max.value)) +
    labs(x = NULL, y = "Luminescence", fill = NULL) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      #legend.position = c(0.02, 0.98),
      #legend.justification = c(0, 1),
      legend.background = element_rect(fill = NA)
    )
  
  p1 + p2
}

dot.plot <- function(df, y_var = "rel.luminescence", strongest = 0, threshold_z_dmso = 5, threshold_z_pos = 5, threshold_lum = 2, pos) {
  
  allowed <- c("rel.luminescence", "to.positive")
  
  if (!y_var %in% allowed) {
    stop("y_var must be one of: ", paste(allowed, collapse = ", "))
  }
  
  if (y_var == "rel.luminescence") {
    text <- "DMSO"
  } else {
    text <- pos
  }

  if (strongest == 0) {
    if (y_var == "rel.luminescence") {
      hit_idx <- df$robust_z_dmso > threshold_z_dmso & df$rel.luminescence > threshold_lum
    } else {
      hit_idx <- df$robust_z_pos > threshold_z_pos & df$to.positive > threshold_lum
    }
    layer1 <- df[hit_idx, ]
    layer2 <- df[!hit_idx, ]
  } else {
    df <- df[order(df[[y_var]], decreasing = TRUE), ]
    layer1 <- df[1:strongest, ]
    layer2 <- df[(strongest + 1):nrow(df), ]
  }
  
  number_comp <- nrow(df)
  offset_y <- 0.02 * max(df[[y_var]])
  offset_x <- 0.01 * number_comp
  
  ggplot() +
    geom_point(data = layer2, aes(x = CatalogID, y = .data[[y_var]]), col = "grey50") +
    geom_point(data = layer1, aes(x = CatalogID, y = .data[[y_var]]), col = "coral1") +
    xlab(label = "Compounds") +
    ylab(label = "rel. Luminescence") +
    geom_hline(yintercept=1, linetype="dashed", color = "red") +
    annotate("text", x=number_comp - offset_x, y= 1 + offset_y, hjust = 1, vjust = 0, label= text, col = "red", size = 6) +
    theme_classic() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank()) 
}

get.export.df <- function(df, hits = "positive", strongest = 0, threshold_z_dmso = 5, threshold_z_pos = 5, threshold_lum = 2) {
  if (hits == "positive") {
    subset(df, df$robust_z_pos > threshold_z_pos & df$to.positive > threshold_lum)
  }
  else if (hits == "dmso") {
    subset(df, df$robust_z_dmso > threshold_z_dmso & df$rel.luminescence > threshold_lum)
  }
  else if (hits == "all") {
    df
  } else if (hits == "best") {
    df <- df[order(df$rel.luminescence, decreasing = TRUE),]
    df <- df[1:strongest, ]
  }
}


resolve_smiles <- function(name) {
  res <- webchem::opsin_query(name)
  
  if (!is.na(res$smiles) && res$status == "SUCCESS") {
    res$smiles
  } else {
    NA_character_
  }
}

#draw_smiles <- function(smile) {
#  mol <- Chem$MolFromSmiles(smile)
#  d <- Draw$rdMolDraw2D$MolDraw2DSVG(300L, 200L)
#  d$DrawMolecule(mol)
#  d$FinishDrawing()
#  return(d$GetDrawingText())
#}

calculate_z_score <- function(df, control.col) {
  control_values <- df$Value[df$Column == control.col]
  control_value_median <- median(control_values, na.rm = TRUE)
  control_mad_sd <- mad(control_values, constant = 1.4826, na.rm = TRUE)
  
  if (control_mad_sd == 0 || is.na(control_mad_sd)) {
    qc_warning(paste("Dead control on plate", unique(df$PlateID)), "error")
    return(rep(NA_real_, nrow(df)))         
  }
  
  return((df$Value - control_value_median) / control_mad_sd)
}

#############
# QC Checks #
#############

# Helper to identify severity of warnings
qc_warning <- function(msg, severity = c("error", "warning")) {
  severity <- match.arg(severity)
  cond <- structure(
    class = c(paste0("qc_", severity), "qc_condition", "warning", "condition"),
    list(message = msg, call = NULL)
  )
  warning(cond)
}

check_plate <- function(df, plate.neg.col = 1, plate.pos.col = 24) {
  pos_values <- df$Value[df$Column == plate.pos.col]
  neg_values <- df$Value[df$Column == plate.neg.col]
  
  pos_mean <- mean(pos_values, na.rm = TRUE)
  neg_mean <- mean(neg_values, na.rm = TRUE)
  pos_sd   <- sd(pos_values, na.rm = TRUE)
  neg_sd   <- sd(neg_values, na.rm = TRUE)
  
  window <- abs(pos_mean - neg_mean)
  
  z_prime_score <- if (is.na(window) || window == 0) {
    NA_real_
  } else {
    1 - (3 * (pos_sd + neg_sd)) / window
  }
  
  plate_id <- unique(df$PlateID)
  
  if (is.na(z_prime_score)) {
    qc_warning(paste("Plate", plate_id, "- could not compute Z' (missing or flat controls)."), "error")
  } else if (z_prime_score < 0) {
    qc_warning(paste("Plate", plate_id, "- plate controls overlap; no usable assay window."), "error")
  } else if (z_prime_score < 0.5) {
    qc_warning(paste("Plate", plate_id, "- marginal assay window; interpret hits with caution."), "warning")
  }
  
  return(z_prime_score)
}

robust_Zprime_check_plate <- function(df, plate.neg.col = 1, plate.pos.col = 24) {
  pos_values <- df$Value[df$Column == plate.pos.col]
  neg_values <- df$Value[df$Column == plate.neg.col]

  pos_mean <- mean(pos_values, na.rm = TRUE)
  neg_mean <- mean(neg_values, na.rm = TRUE)
  pos_sd   <- sd(pos_values, na.rm = TRUE)
  neg_sd   <- sd(neg_values, na.rm = TRUE)

  window <- abs(pos_mean - neg_mean)

  z_prime_score <- if (is.na(window) || window == 0) {
    NA_real_
  } else {
    1 - (3 * (pos_sd + neg_sd)) / window
  }

  plate_id <- unique(df$PlateID)

  if (is.na(z_prime_score)) {
    warning(paste("Plate", plate_id, "- could not compute Z' (missing or flat controls)."))
  } else if (z_prime_score < 0) {
    warning(paste("Plate", plate_id, "- plate controls overlap; no usable assay window."))
  } else if (z_prime_score < 0.5) {
    warning(paste("Plate", plate_id, "- marginal assay window; interpret hits with caution."))
  }

  return(z_prime_score)
}

check_wells <- function(df) {
  df_plate_pos <- subset(df, df$Column == 1 | df$Column == 24)
  df_plate_pos$z_plate_pos <- calculate_z_score(df_plate_pos, 1)
}


