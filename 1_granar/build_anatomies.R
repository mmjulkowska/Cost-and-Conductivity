# Build v5 GRANAR anatomies from 1_granar/granar_params/ into <outputs>/granar/anatomies/ (run from anywhere).
# Output folder = $COSTCOND_OUT, default <repo>/outputs (same rule as paths.py).
# Seed = 10 x section + cell, plus an offset from seed_overrides.txt if listed; if GRANAR fails, retries with +100, +200, ... (up to +1000)
# and records the offset used. Optional args: ids to (re)build, e.g.  Rscript 1_granar/build_anatomies.R MP8_cell3
suppressMessages(library(granar))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])))
OUTROOT <- Sys.getenv("COSTCOND_OUT", file.path(dirname(HERE), "outputs"))
IN <- file.path(HERE, "granar_params"); OUT <- file.path(OUTROOT, "granar", "anatomies"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
args <- commandArgs(TRUE)
ovf <- file.path(HERE, "seed_overrides.txt")
ov <- if (file.exists(ovf) && file.info(ovf)$size>0) read.table(ovf, col.names=c("id","off"), stringsAsFactors=FALSE) else data.frame(id=character(), off=integer())
ids <- if (length(args)) args else sub("\\.xml$","", sort(list.files(IN, pattern="^MP\\d+_cell\\d+\\.xml$")))
res <- list()
for (id in ids) {
  n <- as.integer(sub("MP(\\d+)_.*","\\1",id)); k <- as.integer(sub(".*cell","",id))
  p <- read_param_xml(file.path(IN, paste0(id,".xml")))
  off0 <- if (id %in% ov$id) ov$off[ov$id==id][1] else 0; sim <- NULL; used <- NA
  for (off in unique(c(off0, seq(0,1000,100)))) {
    set.seed(n*10 + k + off)
    sim <- tryCatch(suppressWarnings(create_anatomy(parameters=p, verbatim=FALSE)), error=function(e) NULL)
    if (!is.null(sim)) { used <- off; break }
  }
  if (is.null(sim)) { cat(id, "FAILED all seeds\n"); res[[id]] <- data.frame(id=id, ok=FALSE, seed_offset=NA); next }
  invisible(capture.output(write_anatomy_xml(sim, file.path(OUT, paste0(id,"_root.xml")))))
  o <- sim$output; g <- function(nm,tp){ v <- o$value[o$name==nm & o$type==tp]; if (length(v)) v[1] else NA }
  res[[id]] <- data.frame(id=id, ok=TRUE, seed_offset=used, area_all=g("all","layer_area"), area_cortex=g("cortex","layer_area"),
     area_exo=g("exodermis","layer_area"), area_epi=g("epidermis","layer_area"), area_stelar=g("stelar","layer_area"),
     n_cortex_cells=g("cortex","n_cells"), n_exo_cells=g("exodermis","n_cells"), cortex_cell_area=g("cortex","cell_area"))
}
S <- do.call(rbind, lapply(res, function(d){ for (c in c("seed_offset","area_all","area_cortex","area_exo","area_epi","area_stelar","n_cortex_cells","n_exo_cells","cortex_cell_area")) if (!c %in% names(d)) d[[c]] <- NA; d }))
f <- file.path(OUTROOT, "granar", "granar_anatomy_summary.csv")
if (length(args) && file.exists(f)) { old <- read.csv(f); S <- rbind(old[!old$id %in% S$id,], S) }
write.csv(S[order(S$id),], f, row.names=FALSE)
new_ov <- S[S$ok & !is.na(S$seed_offset) & S$seed_offset>0, c("id","seed_offset")]
write.table(new_ov, ovf, row.names=FALSE, col.names=FALSE, quote=FALSE)
cat(sum(S$ok), "of", nrow(S), "anatomies built;", nrow(new_ov), "needed a seed offset\n")
