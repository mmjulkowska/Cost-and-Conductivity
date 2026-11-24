# Does Krs follow genetic structure more under salt than in control?
# Reads the 5 PCs and the INT phenotypes of the 232-accession run; reports
#  (1) variance in each trait explained by the 5 PCs, (2) per-PC correlations,
#  (3) trait means by genetic group (k-means on PCs), (4) EMMAX pseudo-h2 from the .reml files.
W  <- path.expand(Sys.getenv("WORK", "~/pareto_gwas_krs_v2"))   # the main Krs GWAS folder on the server
pc <- read.table(path.expand("~/pimp_gwas_for_ruth/genotypes/ALLCHR.filt.prune.PCA.txt"))
pc <- pc[, c(1, 2, which(apply(pc[, -(1:2)], 2, var) > 0) + 2)]          # drop the intercept column
names(pc) <- c("FID", "IID", paste0("PC", 1:(ncol(pc) - 2)))
traits <- c("logKrs.C", "logKrs.S", "logKrs.SminusC", "KrsResid.C", "KrsResid.S")
for (t in traits) {
  p <- read.table(file.path(W, "pheno/emmax", paste0("trait.", t, ".INT.pheno")), na.strings = "NA")
  pc[[t]] <- p$V3[match(pc$IID, p$V2)]
}
pcs <- grep("^PC", names(pc), value = TRUE)
cat("\n== (1) variance explained by the 5 PCs (R2, and P of the joint F-test)\n")
for (t in traits) {
  f <- summary(lm(reformulate(pcs, t), data = pc)); fs <- f$fstatistic
  cat(sprintf("%-16s n=%3d  R2=%.3f  P=%.2g\n", t, sum(!is.na(pc[[t]])), f$r.squared,
              pf(fs[1], fs[2], fs[3], lower.tail = FALSE)))
}
cat("\n== (2) Pearson r of each trait with each PC\n")
print(round(sapply(traits, function(t) sapply(pcs, function(k) cor(pc[[k]], pc[[t]], use = "complete.obs"))), 2))
cat("\n== (3) genetic groups: k-means (k = 3) on PC1-3, trait means (INT units) and ANOVA P\n")
set.seed(1); pc$group <- factor(kmeans(pc[, pcs[1:3]], 3, nstart = 50)$cluster)
print(table(group = pc$group))
for (t in traits) {
  m <- tapply(pc[[t]], pc$group, mean, na.rm = TRUE)
  a <- anova(lm(pc[[t]] ~ pc$group))$`Pr(>F)`[1]
  cat(sprintf("%-16s %s   P=%.2g\n", t, paste(sprintf("g%s=%+.2f", names(m), m), collapse = "  "), a))
}
cat("\n== (4) EMMAX pseudo-heritability Vg/(Vg+Ve), from the CHR01 .reml\n")
for (t in traits) {
  r <- scan(file.path(W, "out", paste0(t, ".CHR01.reml")), quiet = TRUE)
  cat(sprintf("%-16s h2 = %.2f\n", t, r[4] / (r[4] + r[5])))
}
write.table(pc, file.path(W, "structure", "pcs_with_krs_traits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
