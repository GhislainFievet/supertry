function editSelection(clicked_id) {
    let int_clicked_id = clicked_id.split("_")[1]
    Shiny.setInputValue("in_edit_selection", int_clicked_id, {priority: "event"});
}

function removeSelection(clicked_id) {
    let int_clicked_id = clicked_id.split("_")[1]
    Shiny.setInputValue("selection_to_remove", int_clicked_id, {priority: "event"});
}

function setFileToDL(clicked_id) {
    Shiny.setInputValue("file_to_dl", clicked_id, {priority: "event"});
}

function doDLFile(clicked_id){
    let int_clicked_id = clicked_id.split("_")[1]
    Shiny.setInputValue("file_to_dl", int_clicked_id, {priority: "event"});
    $('#downloadData')[0].click()
}
