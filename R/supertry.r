library(shiny)
library(ggplot2)
library(Cairo)
library(Seurat)
library(shinycssloaders)
library(plotly)
library(shinyjs)
library(DT)
library(tidyverse)

func_supertry = function(){
    message("hello 3")
}

runInitJS <- function(){
    shinyjs::runjs('document.getElementById("sel.panel.mode").parentElement.style.display="none"')
}

ui <- function(seurat.object) {
    fluidPage(
        useShinyjs(),
        shiny::includeScript(file.path("..","js","script.js")),
        
        # titlePanel(UITitlePanel()),
        UITitlePanel(),
        fluidRow(
            column(
                width = 8,
                uiOutput("sel.panel"),
            ),
            column(
                width = 4,
                UIVisPanel(),
                UIListSelection(),
                style="background-color:rgb(240,240,240)"
            ),
        ),
        UIVisSelection(),
        conditionalPanel(
            "false", # always hide the download button
            downloadButton("downloadData")
        )
    )
}
