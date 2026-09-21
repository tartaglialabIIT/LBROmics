args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: workflow outfile pkglist.csv pkg1 pkg2 ...")
workflow <- args[[1]]
outfile <- args[[2]]
pkglist_file <- args[[3]]
pkgs <- args[-(1:3)]

hdr <- c(
  "# LIVE sessionInfo() dump — library() imports only (analysis scripts not fully executed)",
  paste0("# Workflow: ", workflow),
  paste0("# Captured: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("# Host: ", Sys.info()[["nodename"]], " / ", Sys.info()[["sysname"]], " ", Sys.info()[["release"]], " ", Sys.info()[["machine"]]),
  paste0("# R.home: ", R.home()),
  paste0("# .libPaths: ", paste(.libPaths(), collapse = " | ")),
  paste0("# Bioconductor: ", if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "BiocManager not installed"),
  "#",
  "# CAUTION: Versions below are from THIS machine at capture time.",
  "# They will often differ from manuscript Methods (e.g. DESeq2 1.30.1, Seurat 4.1.0).",
  "# Use for documenting a working modern stack / reviewer install aid; not bit-for-bit paper reproduction.",
  "#"
)

ok <- character(); fail <- character(); fail_msg <- character()
for (p in pkgs) {
  res <- tryCatch({
    suppressPackageStartupMessages(library(p, character.only = TRUE))
    list(ok = TRUE, msg = NA_character_)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (isTRUE(res$ok)) ok <- c(ok, p) else {
    fail <- c(fail, p); fail_msg <- c(fail_msg, res$msg)
  }
}

lines <- c(
  hdr,
  paste0("# Requested packages: ", paste(pkgs, collapse = ", ")),
  paste0("# Loaded OK (", length(ok), "): ", paste(ok, collapse = ", ")),
  paste0("# Failed (", length(fail), "): ",
         if (length(fail)) paste(sprintf("%s [%s]", fail, fail_msg), collapse = "; ") else "none"),
  "#",
  "----- sessionInfo() -----",
  capture.output(sessionInfo())
)
writeLines(lines, outfile)

ip <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
want <- unique(c(pkgs, ok))
want <- want[want %in% ip$Package]
out <- ip[match(sort(want), ip$Package), , drop = FALSE]
write.csv(out, pkglist_file, row.names = FALSE)
message(workflow, ": OK=", length(ok), " FAIL=", length(fail), " -> ", outfile)
quit(status = if (length(ok) == 0) 1 else 0)
