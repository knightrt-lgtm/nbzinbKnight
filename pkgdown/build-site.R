# Build the pkgdown site into docs/ and leave CLAUDE.md out of it.
#
# pkgdown renders every top-level .md file except README/LICENSE/NEWS, and
# has no setting to skip one. CLAUDE.md (the Claude Code project rules)
# belongs in the repository root, so instead we build the site normally and
# then remove the CLAUDE page from docs/: the HTML page, its .md copy, and
# its entries in the sitemap, the search index, and llms.txt.
#
# Run from the package root:  Rscript pkgdown/build-site.R
# (On Windows without a system pandoc, first set RSTUDIO_PANDOC to the
# folder containing pandoc.exe, e.g. RStudio's bundled Quarto tools folder.)

pkgdown::build_site(preview = FALSE)

unlink(file.path("docs", c("CLAUDE.html", "CLAUDE.md")))

drop_lines <- function(file) {
  if (file.exists(file)) {
    x <- readLines(file, warn = FALSE)
    writeLines(x[!grepl("CLAUDE", x, fixed = TRUE)], file)
  }
}
drop_lines(file.path("docs", "sitemap.xml"))
drop_lines(file.path("docs", "llms.txt"))

search <- file.path("docs", "search.json")
if (file.exists(search)) {
  s <- jsonlite::read_json(search)
  entry_path <- function(e) {
    p <- e[["path"]]
    if (length(p) == 0) "" else as.character(p[[1]])
  }
  keep <- !grepl("CLAUDE", vapply(s, entry_path, character(1)), fixed = TRUE)
  jsonlite::write_json(s[keep], search, auto_unbox = TRUE)
}

message("Site built in docs/ (CLAUDE.md page removed).")
