# server.r
landing_server <- function(input, output, session) {
  
  current_tool <- reactiveVal(NULL)        # NULL = show landing page
  
  # buttons on the landing page set the choice
  observeEvent(input$go_nanobit, current_tool("nanobit"))
  observeEvent(input$go_growth,  current_tool("growth"))
  
  # back button returns to the menu
  observeEvent(input$back, current_tool(NULL))
  
  # fill the slot: landing page when NULL, module when chosen
  output$main_ui <- renderUI({
    if (is.null(current_tool())) {
      # ---- the landing page content (tagList, NOT page_fluid) ----
      tagList(
        titlePanel("Lab Analysis Tools"),
        #br(""),
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("NanoBiT Compound Screen"),
            card_body(
              p("Match compound plates to measurements, explore hits, and export results."),
              actionButton("go_nanobit", "Open", class = "btn-primary")
            )
          ),
          card(
            card_header("Growth Curve Analysis"),
            card_body(
              p("Analyse growth curves from your plate reader output."),
              actionButton("go_growth", "Open", class = "btn-primary")
            )
          )
        )
      )
    } else if (current_tool() == "nanobit") {
      tagList(
        actionButton("back", "← Back to menu", class = "btn-secondary"),
        mod_nanobit_ui("nanobit")
      )
    } else if (current_tool() == "growth") {
      tagList(
        actionButton("back", "← Back to menu", class = "btn-secondary"),
        mod_growthcurver_ui("growth")
      )
    }
    
  })
  
  # module server — called once, sits ready
  mod_nanobit_server("nanobit")
  mod_growthcurver_server("growth")
}