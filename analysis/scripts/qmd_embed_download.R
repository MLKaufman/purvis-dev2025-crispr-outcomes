qmd_embed_download <- function() {
    current_file <- knitr::current_input()
    # strip the extension and replace with .qmd
    if (!is.null(current_file)) {
        current_file <- sub("\\.[^.]+$", ".qmd", current_file)
        # Clean up any extra underscores or special characters from the filename
        current_file <- basename(current_file)  # Get just the filename, not the path
        current_file <- gsub("^_+|_+$", "", current_file)  # Remove leading/trailing underscores
    }

    download_filename <- current_file  # Don't use URLencode for the filename attribute
    qmd_content <- readLines(current_file, warn = FALSE)
    qmd_text <- paste(qmd_content, collapse = "\n")
    qmd_encoded <- base64enc::base64encode(charToRaw(qmd_text))
    download_link <- paste0("data:application/octet-stream;base64,", qmd_encoded)

    cat('<div style="display: flex; justify-content: space-between; align-items: center; margin-top: 10px; margin-bottom: 20px;">')
    cat(
        '<a href="', download_link, '" download="', download_filename,
        '" style="display: inline-block; padding: 10px 20px; background-color: #007bff; color: white; text-decoration: none; border-radius: 5px;">', download_filename, '</a>'
    )
    cat("</div>")
}
