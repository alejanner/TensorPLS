test_that("synthetic ORA preserves universes, boundaries, tails and BH", {
  skip_if(Sys.getenv("TENSORPLS_ENRICHMENT_INTEGRATION") != "true", "Opt-in database integration")
  skip_if_not_installed("clusterProfiler")
  local_mocked_bindings(.enrichment_annotation = function(ids, ...) {
    list(pairs = data.frame(input = ids, entrez_id = ids), recognized = ids, metadata = data.frame())
  }, .package = "TensorPLS")
  universe <- as.character(1:600)
  genes <- as.character(1:12)
  sets <- list(A = 1:10, B = 5:24, C = 1:500, small = 1:9,
               large = 1:501, zero = 100:120, outside = 601:620)
  h <- do.call(rbind, lapply(names(sets), function(nm) data.frame(term = nm, gene = as.character(sets[[nm]]))))
  z <- tensorpls_enrichment(genes, universe, databases = "Hallmark", hallmark_term2gene = h)
  t <- z$results$Hallmark$full
  expect_setequal(t$ID, c("A", "B", "C"))
  expect_setequal(z$results$Hallmark$effective_universe, as.character(1:501))
  expect_identical(z$results$Hallmark$effective_genes, genes)
  p <- vapply(t$ID, function(nm) phyper(length(intersect(genes, as.character(sets[[nm]]))) - 1,
    length(sets[[nm]]), 501 - length(sets[[nm]]), 12, lower.tail = FALSE), numeric(1))
  expect_equal(t$pvalue, unname(p), tolerance = 1e-12)
  expect_equal(t$p.adjust, unname(p.adjust(p, "BH")), tolerance = 1e-12)
  expect_identical(z$results$Hallmark$significant$ID, t$ID[t$p.adjust < .05])
  expect_identical(z$parameters$pAdjustMethod, "BH")
  expect_true(all(t$GeneRatio == paste0(t$Count, "/12")))
  expect_true(all(t$BgRatio == paste0(lengths(sets[t$ID]), "/501")))
  unannotated <- tensorpls_enrichment(c(genes, "600"), universe, databases = "Hallmark", hallmark_term2gene = h)
  expect_identical(unannotated$results$Hallmark$effective_genes, genes)
  expect_equal(unannotated$results$Hallmark$full$pvalue, t$pvalue)
  same <- tensorpls_enrichment(c(genes, genes), c(universe, universe), databases = "Hallmark", hallmark_term2gene = h)
  expect_identical(z$results$Hallmark$full, same$results$Hallmark$full)
  null <- tensorpls_enrichment(as.character(1:501), universe, databases = "Hallmark", hallmark_term2gene = h)
  expect_identical(null$results$Hallmark$status, "tested")
  expect_equal(nrow(null$results$Hallmark$significant), 0)
  expect_true(all(null$results$Hallmark$full$pvalue == 1))
})

test_that("installed human annotation maps exact symbols without alias rescue", {
  skip_if(Sys.getenv("TENSORPLS_ENRICHMENT_INTEGRATION") != "true", "Opt-in database integration")
  skip_if_not_installed("org.Hs.eg.db")
  a <- .enrichment_annotation(c("STAT1", "MX1", "ISG15", "IRF7", "stat1", "NOT_A_GENE"), "SYMBOL")
  expect_setequal(a$pairs$input, c("STAT1", "MX1", "ISG15", "IRF7"))
  expect_setequal(a$pairs$entrez_id, c("6772", "4599", "9636", "3665"))
  org <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  aliases <- AnnotationDbi::select(org, "7157", "ALIAS", "ENTREZID")
  expect_true("P53" %in% aliases$ALIAS)
  expect_equal(nrow(.enrichment_annotation("P53", "SYMBOL")$pairs), 0)
})

test_that("database universes require exact membership, not matching counts", {
  skip_if(Sys.getenv("TENSORPLS_ENRICHMENT_INTEGRATION") != "true", "Opt-in database integration")
  skip_if_not_installed("clusterProfiler")
  local_mocked_bindings(.enrichment_annotation = function(ids, ...) {
    list(pairs = data.frame(input = ids, entrez_id = ids), recognized = ids, metadata = data.frame())
  }, .package = "TensorPLS")
  h <- data.frame(term = "A", gene = as.character(1:20))
  z <- tensorpls_enrichment(as.character(1:5), as.character(1:30), databases = "Hallmark", hallmark_term2gene = h)
  bad <- z$results$Hallmark$native
  bad@universe[1] <- "99999"
  local_mocked_bindings(.enrichment_run = function(...) bad, .package = "TensorPLS")
  expect_error(tensorpls_enrichment(as.character(1:5), as.character(1:30), databases = "Hallmark",
    hallmark_term2gene = h), "Backend universe")
})

test_that("real ambiguous symbols keep every candidate and are excluded", {
  skip_if(Sys.getenv("TENSORPLS_ENRICHMENT_INTEGRATION") != "true", "Opt-in database integration")
  skip_if_not_installed("org.Hs.eg.db")
  org <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  pairs <- AnnotationDbi::select(org, AnnotationDbi::keys(org, keytype = "SYMBOL"),
    "ENTREZID", "SYMBOL")
  counts <- table(unique(pairs)$SYMBOL)
  symbol <- names(counts[counts > 1L])[1L]
  expect_false(is.na(symbol))
  a <- .enrichment_annotation(symbol, "SYMBOL")
  m <- .enrichment_map(symbol, "genes", a, "SYMBOL")
  expect_true(nrow(m) > 1)
  expect_setequal(m$entrez_id, unique(pairs$ENTREZID[pairs$SYMBOL == symbol]))
  expect_true(all(m$status == "ambiguous_symbol"))
})
