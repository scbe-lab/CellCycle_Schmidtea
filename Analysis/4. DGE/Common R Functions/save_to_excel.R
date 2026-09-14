require(openxlsx)
save_to_excel <- function(list_of_objects, output_file_name, rowNames = TRUE, ... ){
  if(is.null(names(list_of_objects))) stop("E: objects are not named. Cannot place on excel sheets. Please revise your input.")
  
  library(openxlsx)
  
  # Create a blank workbook
  OUT <- createWorkbook()
  
  # Add sheets to the workbook and write the data to the sheets
  for( i in names(list_of_objects)){
    message("adding sheet ", i, " to workbook for file ", output_file_name)
    addWorksheet(OUT, i)
    writeData(OUT, sheet = i, x = list_of_objects[[i]], rowNames = rowNames,...)
  }
  
  # Export the file
  saveWorkbook(OUT, output_file_name)
  
  message("saved to ", output_file_name)
  
}