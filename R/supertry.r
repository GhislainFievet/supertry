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
    message("hello")
}

runInitJS <- function(){
    runjs('document.getElementById("sel.panel.mode").parentElement.style.display="none"')
}