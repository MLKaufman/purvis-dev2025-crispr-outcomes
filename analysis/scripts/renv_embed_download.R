renv_embed_download <- function(lock_file = "renv.lock") {
    if (!file.exists(lock_file)) {
        cat('<div style="margin: 10px 0; padding: 10px; border: 1px solid #ff6b6b; border-radius: 5px; background-color: #ffe0e0; color: #d63031;">')
        cat("<strong>⚠️ Notice:</strong> renv.lock file not found. Run <code>renv::snapshot()</code> to create it.")
        cat("</div>")
    } else {
        download_filename <- URLencode(lock_file, reserved = TRUE)
        lock_content <- readLines(lock_file, warn = FALSE)
        lock_text <- paste(lock_content, collapse = "\n")
        lock_encoded <- base64enc::base64encode(charToRaw(lock_text))
        download_link <- paste0("data:application/octet-stream;base64,", lock_encoded)

        cat('<div style="margin-top: 10px; margin-bottom: 20px;">')
        cat(
            '<a href="', download_link, '" download="', download_filename,
            '" style="display: inline-block; padding: 10px 20px; background-color: #28a745; color: white; text-decoration: none; border-radius: 5px;">renv.lock</a>'
        )
        cat("</div>")
    }
}
