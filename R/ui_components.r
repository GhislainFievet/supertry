UITitlePanel <- function(){
    fluidRow(
        column(
            width=8,
            h1("SeurSelect", align = "center"),
            textInput("sel.panel.mode", "", value="main")
        ),
        column(
            width=4,
            actionButton("quit.button.1","Back to console", icon = icon("arrow-left", style="margin-right:1em"), align="right"),
            style="display:flex;justify-content:flex-end;margin-top:1em"
        ),
        style="background-color:rgb(240,240,240); style:'margin-top:0'"
    )
}

UIVisPanel = function(){
    list(
        fluidRow(
            # column(
            #     width = 3,
            #     selectInput("vis.assays","Assays",names(seurat.object@assays))
            # ),
            column(
                width = 4,
                selectInput("vis.red.algo","Dim Reduc Algo ",c_reducs)
            ),
            column(
                width = 4,
                selectInput("vis.meta.data","Metadata",c_meta_data)
            )
        ),
        withSpinner(plotOutput("vis.plot"))
    )
}

UISelPanel = function(){ 
    list(
        fluidRow(
            column(
                width = 6,
                h2("Selection")
            ),
        ),
        tabsetPanel(
            id="tp.sel.panel",
            tabPanel(
                "DimPlot",
                "",
                fluidRow(                    
                    column(
                        width = 2,
                        selectInput("dp.sel.red.algo","Dim Reduc Algo ", c_reducs)
                    ),
                    column(
                        width = 2,
                        selectInput("dp.sel.meta.data","Metadata", c_meta_data)
                    ),
                    column(
                        width = 1,
                        actionButton("dp.sel.support.valid","OK"),
                        style="display:flex;align-items:center"
                    ),
                    column(
                        width = 5,
                        downloadButton("dp.export.selection","export current selection"),
                        style="display:flex;align-items:center; flex-direction: row-reverse;"
                    ),
                    column(
                        width = 2,
                        actionButton("dp.save.selection","save selection", icon=icon("save"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        style="display:flex;align-items:center"
                    ),
                    style="display:flex"
                )
            ),
            tabPanel(
                "FeatureScatter",
                "",
                fluidRow(
                    column(
                        width = 2,
                        selectInput("fs.sel.meta.data","Metadata", c_meta_data)
                    ),
                    column(
                        width = 2,
                        selectizeInput("fs.sel.genes1", "Gene1", choices = NULL, options = list(placeholder = 'Select gene 1')),
                        # selectInput("fs.sel.genes1","Gene 1", c_genes)
                    ),
                    column(
                        width = 2,
                        selectizeInput("fs.sel.genes2", "Gene2", choices = NULL, options = list(placeholder = 'Select gene 2')),
                        # selectInput("fs.sel.genes2","Gene 2", c_genes)
                    ),
                    column(
                        width = 1,
                        actionButton("fs.sel.support.valid","OK"),
                        style="display:flex;align-items:center"
                    ),
                    column(
                        width = 3,
                        downloadButton("fs.export.selection","export current selection"),
                        style="display:flex;align-items:center; flex-direction: row-reverse;"
                    ),
                    column(
                        width = 2,
                        actionButton("fs.save.selection","save selection", icon=icon("save"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        style="display:flex;align-items:center"
                    ),
                    style="display:flex"
                )
            ),
            tabPanel(
                "VlnPlot",
                "",
                fluidRow(
                    column(
                        width = 2,
                        selectizeInput("vp.sel.genes", "Gene", choices = NULL, options = list(placeholder = 'Select a gene')),
                        # selectizeInput("vp.sel.genes", "Gene", choices = NULL)
                        # selectInput("vp.sel.genes","Gene", c_genes)
                    ),
                    column(
                        width = 1,
                        actionButton("vp.sel.support.valid","OK"),
                        style="display:flex;align-items:center"
                    ),
                    column(
                        width = 7,
                        downloadButton("vp.export.selection","export current selection"),
                        style="display:flex;align-items:center; flex-direction: row-reverse;"
                    ),
                    column(
                        width = 2,
                        actionButton("vp.save.selection","save selection", icon=icon("save"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        style="display:flex;align-items:center"
                    ),
                    style="display:flex"
                )
            ),
            tabPanel(
                "FeaturePlot",
                "",
                fluidRow(
                    column(
                        width = 2,
                        selectInput("fp.sel.red.algo","Dim Reduc Algo ", c_reducs)
                    ),
                    column(
                        width = 2,
                        selectizeInput("fp.sel.genes", "Gene", choices = NULL, options = list(placeholder = 'Select a gene'))
                        # selectInput("fp.sel.genes","Gene", c_genes)
                    ),
                     column(4,
  
                      # Copy the line below to make a slider range 
                      sliderInput("fp.sel.cutoff.slider", label = "Cutoff", min = -1, 
                        max = 5, value = c(0, 4))
                    ),
                    column(
                        width = 1,
                        actionButton("fp.sel.support.valid","OK"),
                        style="display:flex;align-items:center"
                    ),
                    column(
                        width = 3,
                        downloadButton("fp.export.selection","export current selection"),
                        style="display:flex;align-items:center; flex-direction: row-reverse;"
                    ),
                    column(
                        width = 2,
                        actionButton("fp.save.selection","save selection", icon=icon("save"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                        style="display:flex;align-items:center"
                    ),
                    style="display:flex"
                )
            )
        ),
        withSpinner(plotlyOutput("sel.plot", width = "100%", height="100%"))
    )
}

UIVisSelection <- function(){ 
    list(
        h2("Cells"),
        verbatimTextOutput('vis.sel.table')
        # dataTableOutput('vis.sel.table')
        # verbatimTextOutput("click"),
    )
}

UIListSelection <- function(){ 
    list(
        h2("Selection list"),
        div(
        dataTableOutput('list.sel.table'),id="list.sel.table.parent", style="cursor:pointer")
    )
}

UIActionPanel <- function(){ 
    # list(
    #     h2("What do you want to do?"),
    #     actionButton('create.sel.button', "create a selection"),
    #     actionButton('quit.button.2', "close and go back to console")
    # )
    tags$div(        h1("What do you want to do?", style="margin:1em"),
        actionButton('create.sel.button', "Create a selection", style="margin:2em;padding:1em;color: #fff;font-size: x-large; background-color: #337ab7;"),
        actionButton('quit.button.2', "Close and go back to console", style="margin:2em;padding:1em;font-size: x-large;"),
            style="display:flex;flex-direction:column;height:500px")
}