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
    message("hello 4")
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

server <- function(input, output, session) {
    runInitJS()
    
    # updateSelectizeInput(session, "vp.sel.genes", choices = c_genes, server = TRUE)
    # updateSelectizeInput(session, "dp.sel.genes", choices = c_genes, server = TRUE)
    
    
    reactSelList <- reactiveValues(df_lists=l_init_selection$df_lists,
                              c_cell_selections=l_init_selection$c_cell_selections)

    observeEvent(input$quit.button.1,{stopApp(returnValue = list(
        "names"=reactSelList$df_lists[,c("selection","description","cell_number")],
                                      "cells"=reactSelList$c_cell_selections))})
    observeEvent(input$quit.button.2,{stopApp(returnValue = list(
        "names"=reactSelList$df_lists[,c("selection","description","cell_number")],
        "cells"=reactSelList$c_cell_selections))})

    
    observeEvent(input$create.sel.button,{
        updateTextInput(session,"sel.panel.mode", value = "creation")
    })
    
    output$sel.panel <- renderUI({
        if (input$sel.panel.mode == "main"){
            UIActionPanel()
        } else {
            UISelPanel()
        }
    })
    
    output$vis.plot <- renderPlot({
        message("start output$vis.plot")
        req(input$sel.panel.mode)
        req(input$vis.red.algo)
        df_2plot <- as.data.frame(seurat.object@reductions[[input$vis.red.algo]]@cell.embeddings)
        
        c_cell_list = c()
        
        if (input$sel.panel.mode == "main"){
            c_cell_list <- reactCellList()
        } else {
            click_data <- event_data("plotly_selected")
            click_data = NULL
            if(is.null(click_data)){
                c_cell_list <- c()
            } else {
                # curveNumber
                # pointNumber
                # x
                # y
                c_cell_list <- click_data$key
                if ( ! c_cell_list[1] %in% colnames(seurat.object)){
                # if ( can.be.numeric(c_cell_list[1])){
                    c_cell_list <- colnames(seurat.object)[as.numeric(c_cell_list)]
                }
            }
        }
        
        names(df_2plot)[1] <- 'x.seurselect'
        names(df_2plot)[2] <- 'y.seurselect'
        c_new_order <- sample(1:nrow(df_2plot))
        df_2plot <- df_2plot[c_new_order, ]
        c_colors <- seurat.object@meta.data[[input$vis.meta.data]][c_new_order]
        c_alpha <- unlist(lapply(rownames(df_2plot), function(x) if (x %in% c_cell_list){1} else {0.01}))
        message("end output$vis.plot")

        ggplot(df_2plot, aes(x = x.seurselect, y = y.seurselect, color=c_colors)) + geom_point(alpha=c_alpha) + theme_bw()
    })
    
    selDimPlot <- eventReactive(
        input$dp.sel.support.valid,{
            DimPlot(seurat.object, group.by=input$dp.sel.meta.data, reduction=input$dp.sel.red.algo)+
                aes(key=colnames(seurat.object))
    })
    selFeatureScatter <- eventReactive(
        input$fs.sel.support.valid,{
        FeatureScatter(seurat.object, feature1 = input$fs.sel.genes1, feature2 = input$fs.sel.genes2, group.by=input$fs.sel.meta.data)+
            # aes(key=colnames(seurat.object))
            aes(key=1:length(colnames(seurat.object)))
    })
    selVlnPlot <- eventReactive(
        
        input$vp.sel.support.valid,{
            VlnPlot(seurat.object, features = input$vp.sel.genes) + aes(key=1:length(colnames(seurat.object)))
    })
    selFeaturePlot <- eventReactive(
        input$fp.sel.support.valid,{
        FeaturePlot(seurat.object, features = input$fp.sel.genes, min.cutoff=input$fp.sel.cutoff.slider[1], max.cutoff=input$fp.sel.cutoff.slider[2], reduction = input$fp.sel.red.algo) + aes(key=colnames(seurat.object))
    })
        
    
    output$sel.plot <- renderPlotly({
        myplot <- switch(input$tp.sel.panel,
                         "DimPlot"=selDimPlot(),
                         "FeatureScatter"=selFeatureScatter(),
                         "VlnPlot"=selVlnPlot(),
                         "FeaturePlot"=selFeaturePlot()
                        )
        ggplotly(myplot) %>% 
          layout(dragmode = "select")
      })
                          
    ## Handle server side large gene set
    observeEvent(input$tp.sel.panel,
        {switch(input$tp.sel.panel,
                 "DimPlot"={},
                 "FeatureScatter"={updateSelectizeInput(session, "fs.sel.genes1", choices = c_genes, server = TRUE)
                                  updateSelectizeInput(session, "fs.sel.genes2", choices = c_genes, server = TRUE)},
                 "VlnPlot"={updateSelectizeInput(session, "vp.sel.genes", choices = c_genes, server = TRUE)},
                 "FeaturePlot"=updateSelectizeInput(session, "fp.sel.genes", choices = c_genes, server = TRUE)
                )
        })
                          
      ## returns the data related to data points selected by the user
      output$vis.sel.table <- renderPrint({
        c_cell_list = c()
        
        if (input$sel.panel.mode == "main"){
            c_cell_list <- reactCellList()
        } else {
            click_data <- event_data("plotly_selected")
            click_data = NULL
            if(is.null(click_data)){
                c_cell_list <- c()
            } else {
                c_cell_list <- click_data$key
                if ( ! c_cell_list[1] %in% colnames(seurat.object)){
                    c_cell_list <- colnames(seurat.object)[as.numeric(c_cell_list)]
                }
            }
        }
        c_cell_list
      })

    # observeEvent(input$sel.support,{
    #     message("observeEvent(input$sel.support")
    #     updateTabsetPanel(session, "tp.sel.panel", selected=input$sel.support)
    # })

    output$list.sel.table <- DT::renderDT({
        reactSelList$df_lists
    },
      escape = FALSE,
      selection = list(mode = "single", target = "row")
    )
                    
    reactCellList <- eventReactive(
        input$list.sel.table_rows_selected,{
        if (is.null(input$list.sel.table_rows_selected)){
            c()
        } else {
            selection_id <- isolate(reactSelList$df_lists$selection[input$list.sel.table_rows_selected])
            isolate(reactSelList$c_cell_selections[[selection_id]])
        }
    }, ignoreNULL = FALSE)
                          
    observeEvent(input$save.cancel,{
        updateTextInput(session,"sel.panel.mode", value = "main")
    })

    # Observe save selection
    observeEvent(input$dp.save.selection,{
        req(input$dp.save.selection)
        showModal(dataModal())
    })
    observeEvent(input$fs.save.selection,{
        req(input$fs.save.selection)
        showModal(dataModal())
    })
    observeEvent(input$vp.save.selection,{
        req(input$vp.save.selection)
        showModal(dataModal())
    })
    observeEvent(input$fp.save.selection,{
        req(input$fp.save.selection)
        showModal(dataModal())
    })


    observeEvent(input$list.sel.table_rows_selected, {
        req(input$list.sel.table_rows_selected)
        updateTextInput(session,"sel.panel.mode", value = "main")
        # showModal(modalDialog(
        #   title = "message",
        #   paste("This is a somewhat important message:", 
        #         df_lists$selection[input$list.sel.table_rows_selected]),
        #   easyClose = TRUE,
        #   footer = NULL))
      })
                          
   observeEvent(input$current_id, {
           c_str_split = unlist(strsplit(input$current_id,"_"))
           action_type = c_str_split[1]
           action_id = c_str_split[2]
       
          showModal(modalDialog(
          title = "message",
          str(session),
          easyClose = TRUE,
          footer = NULL))
       })
                          
    observeEvent(input$in_edit_selection,{
        showModal(editModal())
        str_title = reactSelList$df_lists[as.numeric(input$in_edit_selection), "selection"][1]
        str_description = reactSelList$df_lists[as.numeric(input$in_edit_selection), "description"][1]
        updateTextInput(session,"edit.name", value = str_title)
        updateTextInput(session,"edit.description", value = str_description)
    })
                          
    observeEvent(input$edit.cancel,{
        removeModal()
    })
                          
    observeEvent(input$edit.form.submit,{
        reactSelList$df_lists[as.numeric(input$in_edit_selection),"selection"] = input$edit.name
        reactSelList$df_lists[as.numeric(input$in_edit_selection),"description"] = input$edit.description
        removeModal()
    })

   observeEvent(
        input$save.selection.form.submit,{

        if (is.null(input$save.selection.form.submit)){
            reactSelList$df_lists <- isolate(reactSelList$df_lists)
            reactSelList$c_cell_selections <- isolate(reactSelList$c_cell_selections)
        } else {
            c_temp <- event_data("plotly_selected")$key
            if ( ! c_temp[1] %in% colnames(seurat.object)){
                # if ( can.be.numeric(c_cell_list[1])){
                    c_temp <- colnames(seurat.object)[as.numeric(c_temp)]
                }
          #   showModal(modalDialog(
          # title = "message",
          # paste("This is a somewhat important message:", 
          #       c_temp[1]),
          # easyClose = TRUE,
          # footer = NULL))
            i_new_id <- nrow(reactSelList$df_lists) + 1
            reactSelList$c_cell_selections[[input$selection.name]] <- c(c_temp)
            reactSelList$df_lists[i_new_id,] <- c(input$selection.name,
                                                 input$selection.description,
                                                 length(c_temp),
                                                  # onclick=$('#downloadData')[0].click()
                                                 paste0("<button id='iddl_",i_new_id,"' onclick=doDLFile(this.id)><i class='fa fa-download'></i></button> ",
                                     "<button id='idedit_",i_new_id,"' onclick=editSelection(this.id)><i class='fa fa-edit'></i></button> ",
                                     "<button id='iddel_",i_new_id,"' onclick=removeSelection(this.id)><i class='fa fa-trash'></i></button> ")
                                                                        )            
            # reactSelList$df_lists <- isolate(reactSelList$df_lists)
            # c_res$c_cell_selections <- isolate(c_cell_selections)
        }
        removeModal()
    })
                          
    observeEvent(input$selection_to_remove, {
        str_selection_id = reactSelList$df_lists$selection[as.numeric(input$selection_to_remove)]
        
        c_mask = 1:nrow(reactSelList$df_lists)
        c_mask = c_mask[c_mask != as.numeric(input$selection_to_remove)]
        reactSelList$df_lists = reactSelList$df_lists[c_mask,]
        
        
        c_cellselect_mask = names(reactSelList$c_cell_selections)
        c_cellselect_mask = c_cellselect_mask[c_cellselect_mask != str_selection_id]

        reactSelList$c_cell_selections = reactSelList$c_cell_selections[c_cellselect_mask]
    })
                          
      # observeEvent(input$save.selection.form.submit,{
      #     c_sel_ids <- event_data("plotly_selected")$key +
      #     showModal(modalDialog(
      #     title = "message",
      #     paste(input$selection.name, input$selection.description,c_sel_ids[1],c_sel_ids[2],c_sel_ids[3]),
      #     easyClose = TRUE,
      #     footer = NULL))
      # })
    

    # Download handlers
    output$dp.export.selection <- downloadHandler(
      filename = function() {
        "selection.txt"
      },
      content = function(file) {
          c_cell_list = c()
        if (input$sel.panel.mode == "main"){
            c_cell_list <- reactCellList()
        } else {
            click_data <- event_data("plotly_selected")
            if(is.null(click_data)){
                c_cell_list <- c()
            } else {
                c_cell_list <- click_data$key
                if ( ! c_cell_list[1] %in% colnames(seurat.object)){
                    c_cell_list <- colnames(seurat.object)[as.numeric(c_cell_list)]
                }
            }
        }
        write.table(c_cell_list, file, sep="\t", quote=F, row.names=F, col.names=F)
      }
    )
    output$fs.export.selection <- downloadHandler(
      filename = function() {
        "selection.txt"
      },
      content = function(file) {
          c_cell_list = c()
        if (input$sel.panel.mode == "main"){
            c_cell_list <- reactCellList()
        } else {
            click_data <- event_data("plotly_selected")
            if(is.null(click_data)){
                c_cell_list <- c()
            } else {
                c_cell_list <- click_data$key
                if ( ! c_cell_list[1] %in% colnames(seurat.object)){
                    c_cell_list <- colnames(seurat.object)[as.numeric(c_cell_list)]
                }
            }
        }
        write.table(c_cell_list, file, sep="\t", quote=F, row.names=F, col.names=F)
      }
    )
    output$vp.export.selection <- downloadHandler(
      filename = function() {
        "selection.txt"
      },
      content = function(file) {
          c_cell_list = c()
        if (input$sel.panel.mode == "main"){
            c_cell_list <- reactCellList()
        } else {
            click_data <- event_data("plotly_selected")
            if(is.null(click_data)){
                c_cell_list <- c()
            } else {
                c_cell_list <- click_data$key
                if ( ! c_cell_list[1] %in% colnames(seurat.object)){
                    c_cell_list <- colnames(seurat.object)[as.numeric(c_cell_list)]
                }
            }
        }
        write.table(c_cell_list, file, sep="\t", quote=F, row.names=F, col.names=F)
      }
    )
    output$fp.export.selection <- downloadHandler(
      filename = function() {
        "selection.txt"
      },
      content = function(file) {
          c_cell_list = c()
        if (input$sel.panel.mode == "main"){
            c_cell_list <- reactCellList()
        } else {
            click_data <- event_data("plotly_selected")
            if(is.null(click_data)){
                c_cell_list <- c()
            } else {
                c_cell_list <- click_data$key
                if ( ! c_cell_list[1] %in% colnames(seurat.object)){
                    c_cell_list <- colnames(seurat.object)[as.numeric(c_cell_list)]
                }
            }
        }
        write.table(c_cell_list, file, sep="\t", quote=F, row.names=F, col.names=F)
      }
    )
     
    output$downloadData <- downloadHandler(
        filename = function() {
          base_name = reactSelList$df_lists$selection[as.numeric(input$file_to_dl)]
          paste0(base_name,".tsv")
        },
        content = function(file) {
            base_name = reactSelList$df_lists$selection[as.numeric(input$file_to_dl)]
            df = reactSelList$c_cell_selections[[base_name]]
            write.csv(df, file, sep="\t", quote=F, row.names=F, col.names=F)
        }
    )                         
}
