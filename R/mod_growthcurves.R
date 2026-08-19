# Growthcurver screening UI
mod_growthcurver_ui <- function(id, label = "Grwothcurver") {
  gs <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        card(
          card_header("Upload Files:"),
          fileInput(gs("setup"), label = "Setup File", buttonLabel = "Browse...", 
                    multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
          fileInput(gs("plate"), label = "Measurements", buttonLabel = "Browse...", 
                    multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
        ),
        actionButton(gs("match"), label = "Match & Normalize", class = "btn-primary", width = "100%"),
        br("")
      ),  
      mainPanel(
        "More to come here!"
      )
    )
  )
}

mod_growthcurver_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
  })
}