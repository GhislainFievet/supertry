visPlot <- function(seurat.object, reduction, metadata) {
    
    df_2plot <- as.data.frame(seurat.object@reductions[[reduction]]@cell.embeddings)

    c_cell_list <- tail(rownames(df_2plot), 1000)

    names(df_2plot)[1] <- 'x.seurselect'
    names(df_2plot)[2] <- 'y.seurselect'
    c_new_order <- sample(1:nrow(df_2plot))
    df_2plot <- df_2plot[c_new_order, ]
    c_colors <- seurat.object@meta.data[[metadata]][c_new_order]

    c_alpha <- unlist(lapply(rownames(df_2plot), function(x) if (x %in% c_cell_list){1} else {0.01}))

    ggplot(df_2plot, aes(x = x.seurselect, y = y.seurselect, color=c_colors)) + geom_point(alpha=c_alpha) + theme_bw()
}
                             
dataModal <- function() {
      modalDialog(
        textInput("selection.name", "Choose a name for your selection"),
        textInput("selection.description", "Short description of the selection"),

        footer = tagList(
          modalButton("Cancel"),
          actionButton("save.selection.form.submit", "OK")
        )
      )
}
                             
editModal <- function() {
      modalDialog(
        textInput("edit.name", "Choose a name for your selection"),
        textInput("edit.description", "Short description of the selection"),

        footer = tagList(
          modalButton("edit.cancel"),
          actionButton("edit.form.submit", "save modifications")
        )
      )
}