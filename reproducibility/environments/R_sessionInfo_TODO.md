# R environments — still need author export

No `sessionInfo()` / `renv.lock` was found under `Cerase_single_cell`.

Preprint + scripts imply **at least two** R stacks:

| Stack | Hint | Typical packages |
|-------|------|------------------|
| scRNA-seq Seurat / destiny / tradeSeq | Comment “runs in R4.2.2” in trajectory script | Seurat 4.1.0, destiny 3.4.0, slingshot 2.8.0, tradeSeq 1.8.0, condiments 1.6.0, scry 1.6.0 |
| bulk DESeq2 / WebGestalt / karyoploteR | DESeq2 1.30.1 (older Bioconductor) | DESeq2, WebGestaltR 0.4.5, karyoploteR 1.20.3, pheatmap 1.0.12 |
| SAMMY visualization | R 4.4.2 in Methods | Gviz 1.50.0, rtracklayer 1.66.0 |

On each analysis machine, please run:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
if (requireNamespace("BiocManager", quietly=TRUE)) {
  writeLines(as.character(BiocManager::version()), "Bioconductor_version.txt")
}
write.csv(installed.packages()[,c("Package","Version")], "installed_R_packages.csv")
```
