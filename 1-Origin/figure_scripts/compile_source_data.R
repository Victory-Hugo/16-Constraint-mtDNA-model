library(readr)
library(openxlsx)

# directory containing tsv files
tsv_dir <- "final_figures_source_data/"

# get a list of TSV files in the directory
tsv_files <- list.files(path = tsv_dir, pattern = "\\.tsv$", full.names = TRUE)

# create a new Excel workbook
wb <- createWorkbook()

# Loop through each TSV file
for (file in tsv_files) {
  # Read the TSV file into a data frame
  df <- read_tsv(file)
  
  # Get the base name of the file (without directory and extension)
  sheet_name <- tools::file_path_sans_ext(basename(file))
  
  # Add a new sheet to the workbook with the file name as the sheet name
  addWorksheet(wb, sheet_name)
  
  # Write the data frame to the sheet
  writeData(wb, sheet = sheet_name, x = df)
}

# Save the workbook to an XLSX file
output_file <- "final_figures_source_data/mito_constraint_source_data.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("Combined Excel file saved as", output_file)