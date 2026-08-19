# global.R
required_packages <- c("shiny", "bslib", "readxl", "tibble", 
                       "tidyr", "dplyr", "ggplot2", "writexl",
                       "webchem", "shinyFiles", "patchwork",
                       "purrr", "stringr", "scales", "stats", 
                       "viridis", "growthcurver", "reticulate") 

for (pkg in required_packages) suppressMessages(library(pkg, character.only = TRUE))

py_require("rdkit")
source_python(file.path("Python", "cheminformatics.py"))
