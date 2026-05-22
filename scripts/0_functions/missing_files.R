missing_files <- files_to_check[!file.exists(files_to_check)]
if (length(missing_files) > 0) {
  stop(
    "Faltan los siguientes archivos:\n",
    paste(missing_files, collapse = "\n")
  )
}