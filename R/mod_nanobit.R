
# NanoBiT screening UI
mod_nanobit_ui <- function(id, label = "NanoBiT") {
  ns <- NS(id)
  tagList(
    tags$style(HTML("
      .sticky-table {
        max-height: 300px;
        overflow-y: auto;
      }
      .sticky-table thead th {
        position: sticky;
        top: 0;
        background-color: var(--bs-body-bg);
        z-index: 1;
      }
    ")),
    uiOutput(ns("qc_banner")),
    sidebarLayout(
      sidebarPanel(
        card(
          card_header("Upload Files:"),
          fileInput(ns("map"), label = "Compound map", buttonLabel = "Browse...", multiple = TRUE),
          fileInput(ns("plate"), label = "Plate Map", buttonLabel = "Browse...", multiple = TRUE)
        ),
        card(
          card_header("Column"),
          layout_columns(
            col_widths = c(5,5),
            numericInput(ns("dmso_col"), label = "DMSO", value = 2, min = 1, max = 24, step = 1),
            numericInput(ns("pos_col"), label = "Positive", value = 23, min = 1, max = 24, step = 1)
          )
        ),
        actionButton(ns("match"), label = "Match Compounds with Values", class = "btn-primary", width = "100%"),
        br(""),
        actionButton(ns("send_to_console"), label = "Send data to R Studio", class = "btn-primary", width = "100%"),
        helpText("Available in your R console as 'nanobit_data' after closing the app."),
        br(""),
        card(
          card_header("Download Files:"),
          textInput(ns("fname"), label = "File name", value = "nanobit_hits.csv"),
          radioButtons(ns("df"), label = NULL, inline = FALSE,
                       choices = c("Everything" = "all", "Best" = "best"),
                       selected = "all"),
          shinyDirButton(ns("folder"), label = "Choose folder", title = "Select where to save", class = "w-100"),
          actionButton(ns("save"), label = "Save CSV", width = "100%")
        )
      ),
      mainPanel(
        navset_card_underline(
          title = "Plot",
          nav_spacer(),
          nav_panel("Processed",
                    radioButtons(ns("plotVar"), label = NULL, inline = TRUE,
                                 choices = list("relative DMSO" = 1, 
                                                "relative positive Compound" = 2),
                                 selected = 1
                    ),
                    plotOutput(ns("plot"), click = ns("plot_click")),
                    layout_columns(
                      col_widths = c(4,4,4),
                      sliderInput(ns("rel.lum"), label = "rel. Luminescence", min = 0, max = 5, value = 2, step = 0.25),
                      div(
                        conditionalPanel(
                          condition = "input.plotVar == 1",
                          ns = ns,                                  # <- important in a module
                          sliderInput(ns("z.DMSO"), "Z DMSO", min = 0, max = 5, value = 2, step = 0.25)
                        ),
                        conditionalPanel(
                          condition = "input.plotVar == 2",
                          ns = ns,
                          sliderInput(ns("z.POS"), "Z control Compound", min = 0, max = 5, value = 2, step = 0.25)
                        )
                      ),
                      sliderInput(ns("number.hits"), label = "Number of Hits", min = 0, max = 200, value = 0, step = 5)
                    ),
                    card(
                      helpText("Z-score/Z threshold: ≥2 very permissive (many false positives), ≥3 common/moderate, ≥4–5 stringent (good for validation), ≥6–7 very stringent (only the strongest hits).")
                    )
          ),
          nav_panel("Plates",
                    layout_columns(
                      col_widths = c(8, 2, 2),
                      selectInput(ns("barcode_select"), label = "Plate", choices = ""),
                      actionButton(ns("prev_plate"), "←", class = "mt-4"),
                      actionButton(ns("next_plate"), "→", class = "mt-4")
                    ),
                    textOutput(ns("z_prime_display")),                    
                    helpText(HTML(
                      "Z' measures assay quality — separation between plate controls relative to their variability.<br>
                      <b>0.5–1.0</b> robust separation, high-confidence screening<br>
                      <b>0–0.5</b> marginal; usable with caution, more false positives<br>
                      <b>below 0</b> controls overlap, plate unreliable"
                    )),
                    plotOutput(ns("platePlot"))
          ),
          
        ),
        card(
          card_header("Compound"),
          uiOutput(ns("structure"))
        ),
        card(
          card_header("Table"),
          height = 450,
          card_body(
            min_height = 350,
            layout_columns(
              col_widths = c(6,6),
              actionButton(ns("clear"), label = "Clear"),
              selectInput(ns("table_barcode_select"), label = NULL, choices = ""),
            ),
            div(class = "sticky-table", tableOutput(ns("table")))
          )
        )
      )
    )
  )
}
  

# NanoBiT screening server
mod_nanobit_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    #################
    ### reactive Vals
    #################
    
    all_plates <- reactiveVal(NULL)
    plate_z_primes <- reactiveVal(NULL)
    current_table <- reactiveVal(NULL)
    structure_url <- reactiveVal(NULL)
    smiles_val <- reactiveVal(NULL)
    qc_messages <- reactiveVal(list())
    
    var <- reactive({
      if (input$plotVar == 1) {
        "rel.luminescence"
      } else if (input$plotVar == 2) {
        "to.positive"
      }
    })
    
    #################
    ### Data Processing
    #################
    
    # Read Data
    observeEvent(input$match, {
      req(input$map, input$plate)
      all_plates(read_all_plates(input$map$datapath, input$plate$datapath))
    })
    
    matched_data <- reactive({
      req(all_plates())
      
      acc <- new.env(parent = emptyenv())
      acc$msgs <- list()
      
      result <- withCallingHandlers(
        process.plates(all_plates(), dmso = input$dmso_col, pos = input$pos_col),
        warning = function(w) {
          if (!inherits(w, "qc_condition")) return()
          acc$msgs[[length(acc$msgs) + 1]] <- list(
            severity = if (inherits(w, "qc_error")) "error" else "warning",
            text     = conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }
      )
      
      qc_messages(acc$msgs)
      result
    })
    
    # Get Plate with specific Barcode
    fetch_barcode_plate <- eventReactive(input$barcode_select, {
      req(input$barcode_select, all_plates())          
      
      plate <- all_plates()[[input$barcode_select]]
      
      list(
        dmso_pos  = subset(plate, Column %in% c(input$dmso_col, input$pos_col)),
        plate_pos = subset(plate, Column %in% c(1, 24)),
        compounds = subset(plate, !(Column %in% c(1, input$dmso_col, input$pos_col, 24)) & !(Type %in% c("negative", "positive")))
      )
    })
    
    #################
    ### Table
    #################
    
    # Get full Table
    observeEvent(matched_data(), {
      current_table(matched_data())
    })
    
    # Update Table with clicked Point
    observeEvent(input$plot_click, {
      subset <- nearPoints(matched_data(),
                           coordinfo = input$plot_click,
                           xvar = "CatalogID",
                           yvar = var(),
                           maxpoints = 1)
      
      current_table(subset)
    })
    
    # Update Table with clicked Point
    observeEvent(input$table_barcode_select, {
      req(matched_data()) 
      
      subset <- subset(matched_data(), matched_data()$PlateID == input$table_barcode_select)
      
      current_table(subset)
    })
    
    # Clear table selection and update to full Table
    observeEvent(input$clear, {
      req(current_table())
      current_table(matched_data())
    })
    
    # Output current table
    output$table <- renderTable({
      req(current_table())
      current_table()
    })
    
    # Get the name of the positive compound
    get_pos <- reactive({
      req(all_plates())
      plate <- all_plates()[[1]]
      nm <- as.character(plate$Name[plate$Column == input$pos_col][1])
      if (is.na(nm)) "positive Compound" else nm
    })
    
    #################
    ### Plot Data
    #################
    
    # All compounds plot
    dot_plot <- reactive({
      
      z_dmso <- input$z.DMSO
      z_pos <- input$z.POS
      rel_lum <- input$rel.lum
      number_hits <- input$number.hits
      
      df <- matched_data()
      pos <- get_pos()
      
      dot.plot(df, y_var = var(), strongest = number_hits, 
               threshold_z_dmso = z_dmso, threshold_z_pos = z_pos, threshold_lum = rel_lum, pos)
    })
    
    output$plot <- renderPlot({
      dot_plot()
    })
    
    # Plots for each plate by Barcode
    plate_plot <- reactive({
      
      ls <- fetch_barcode_plate()
      dmso_pos <- ls$dmso_pos
      compounds <- ls$compounds
      plate_pos <- ls$plate_pos
      
      plot.plate.bar(compounds, dmso_pos, plate_pos)
    })
    
    output$platePlot <- renderPlot({
      plate_plot()
    })
    
    
    #################
    ### QC
    #################
    
    # Get Z'
    observeEvent(matched_data(), {
      plate_z_primes(run_plate_qc(matched_data()))
    })
    
    output$z_prime_display <- renderText({
      req(input$barcode_select)
      zp <- plate_z_primes()[[input$barcode_select]]    
      paste("Z' for this plate:", round(zp, 2))
    })

    
    #################
    ### Update UI
    #################
    
    # Reset Number Hits slider when other sliders move
    observeEvent(
      list(input$rel.lum, input$z.DMSO, input$z.POS),   
      {
        updateSliderInput(session, "number.hits", value = 0)   
      },
      ignoreInit = TRUE
    )
    
    # Get Barcode and update dropdown menu in the Plates tab
    observeEvent(matched_data(), {
      barcodes <- sort(unique(matched_data()$PlateID))
      
      updateSelectInput(session, "barcode_select", choices = barcodes)
    })
    
    # Get Barcode and update dropdown menu in the Table
    observeEvent(matched_data(), {
      barcodes <- sort(unique(matched_data()$PlateID))
      
      updateSelectInput(session, "table_barcode_select", choices = barcodes)
    })
    
    # Step between plates
    step_plate <- function(direction) {
      barcodes <- sort(unique(matched_data()$PlateID))
      req(length(barcodes) > 0, input$barcode_select)
      i <- match(input$barcode_select, barcodes)
      req(!is.na(i))
      j <- i + direction
      if (j >= 1 && j <= length(barcodes)) {
        updateSelectInput(session, "barcode_select", selected = barcodes[j])
      }
    }
    
    observeEvent(input$prev_plate, step_plate(-1))
    observeEvent(input$next_plate, step_plate(1))
    
    
    #################
    ### Compound structure/SMILE
    #################
    
    # Fetch Compound SMILE and the Structure PNG
    observeEvent(input$plot_click, {
      row <- current_table()
      req(nrow(row) == 1)
      name <- row$Name
      smiles <- resolve_smiles(name)
      
      req(!is.na(smiles))
      #structure_url(
      #  paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", URLencode(smiles, reserved = TRUE), "/PNG?image_size=500x500&bgcolor=white")
      #)
      
      smiles_val(smiles)
    })
    
    # Clear Compound and update to showing nothing
    observeEvent(input$clear, {
      structure_url(NULL)
    })
    
    # Output picture of Structure
    output$structure <- renderUI({
      #req(structure_url())
      req(smiles_val())
      #div(
      #  style = "text-align: center;",
      #  tags$img(src = structure_url(), width = "50%"),
      #  tags$p(smiles_val(), style = "font-family: monospace; word-break: break-all; font-size: 10px;")
      #)
      #source_python(file.path("Python", "cheminformatics.py"))
      svg <- draw_smiles(smiles_val())
      div(
        style = "text-align: center;", 
        HTML(svg),
        tags$p(smiles_val(), style = "font-family: monospace; word-break: break-all; font-size: 10px;")
      )
    })
    
    #################
    ### Download
    #################
    
    # Get directory
    volumes <- c(Home = fs::path_home(), getVolumes()())
    shinyDirChoose(input, "folder", roots = volumes, session = session)
    
    chosen_dir <- reactive({
      parseDirPath(volumes, input$folder)
    })
    
    # Save Data
    observeEvent(input$save, {
      req(chosen_dir(), input$fname)         
      
      full_path <- file.path(chosen_dir(), input$fname)
      
      df_to_write <- get.export.df(
        matched_data(),
        hits            = input$df,          
        strongest       = input$number.hits,
        threshold_z_pos = input$z.POS,         
        threshold_lum   = input$rel.lum
      )
      
      write.csv(df_to_write, file = full_path, row.names = FALSE)
      
      showNotification(paste("Saved to", full_path))
    })
    
    #################
    ### Send data to RStudio
    #################

    observeEvent(input$send_to_console, {
      assign("nanobit_data", matched_data(), envir = globalenv())
      assign("control", fetch_barcode_plate(), envir = globalenv())
      showNotification("Sent to console as 'nanobit_data'")
    })
    
    #################
    ### Update UI
    #################

    # Updates the Plot radioButtons
    observeEvent(get_pos(), {
      pos <- get_pos()
      
      choices = setNames(
        c(1, 2),                                          
        c("relative DMSO", paste("relative", pos))        
      )
      
      updateRadioButtons(session, "plotVar", choices = choices, selected = 1)
    })
    
    # Updates the Z-score slider for the positive compound
    observeEvent(get_pos(), {
      pos <- get_pos()
  
      updateSliderInput(session, "z.POS", paste("Z", pos), value = 2)
    })
    
    # Updates which download radioButton is present - either realtive DMSO or positive Compound
    observeEvent(list(input$plotVar, get_pos()), {
      first_label <- if (input$plotVar == 1) {
        c("relative DMSO" = "dmso")
      } else {
        setNames("positive", paste("relative", get_pos()))
      }
      choices <- c(first_label, "Everything" = "all", "Best" = "best")
      updateRadioButtons(session, "df", choices = choices, selected = "all")
    }, ignoreInit = FALSE)
    
    #################
    ### Output warnings/errors
    #################
    
    output$qc_banner <- renderUI({
      msgs <- qc_messages()
      req(length(msgs) > 0)
      
      errors <- Filter(function(m) m$severity == "error", msgs)
      warns  <- Filter(function(m) m$severity == "warning", msgs)
      
      tagList(
        if (length(errors) > 0)
          div(class = "alert alert-danger",
              tags$b("Critical QC failures"),
              tags$ul(lapply(errors, function(m) tags$li(m$text)))),
        if (length(warns) > 0)
          div(class = "alert alert-warning",
              tags$b("QC warnings"),
              tags$ul(lapply(warns, function(m) tags$li(m$text))))
      )
    })
})}


