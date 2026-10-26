.enrichment_fixture <- function() {
  list(pairs = data.frame(input = c("STAT1", "MX1", "ISG15", "IRF7", "AMB", "AMB", "OLD", "OTHER"),
    entrez_id = c("1", "2", "3", "4", "5", "6", "7", "1")),
    recognized = as.character(1:6), metadata = data.frame())
}

test_that("malformed inputs fail before annotation or fitting", {
  local_mocked_bindings(.enrichment_annotation = function(...) stop("annotation reached"),
    .package = "TensorPLS")
  bad <- c("01", "0", " 1", "1 ", "'1'", '"1"', "1.0", "1e3", "1/2", "1;2", "-1", "+1")
  for (id in bad) expect_error(tensorpls_enrichment(id, id, databases = "GO_BP"), "Malformed")
  for (ids in list(1, NA_character_, ""))
    expect_error(tensorpls_enrichment(ids, ids, databases = "GO_BP"), "character vectors")
  expect_error(tensorpls_enrichment("1", "2", databases = "GO_BP"), "before conversion")
  expect_error(tensorpls_enrichment("1", "1", organism = "mouse"), "human")
  expect_error(tensorpls_enrichment("1", "1", id_type = "alias"), "id_type")
  for (db in list("unknown", character(), c("GO_BP", "GO_BP"), NA_character_))
    expect_error(tensorpls_enrichment("1", "1", databases = db), "databases")
})

test_that("symbol mapping retains all candidates and input occurrences", {
  local_mocked_bindings(.enrichment_annotation = function(...) .enrichment_fixture(),
    .enrichment_run = function(...) NULL, .package = "TensorPLS")
  ids <- c("STAT1", "MX1", "ISG15", "IRF7", "AMB", "OLD", "UNKNOWN", "stat1", "STAT1", "OTHER")
  z <- tensorpls_enrichment(ids, ids, id_type = "SYMBOL", databases = "GO_BP")
  expect_identical(z$original$genes, ids)
  expect_identical(z$converted$genes, as.character(1:4))
  expect_identical(z$converted$universe, as.character(1:4))
  a <- subset(z$mapping, set == "genes")
  expect_identical(a$entrez_id[a$input == "AMB"], c("5", "6"))
  expect_true(all(a$status[a$input == "AMB"] == "ambiguous_symbol"))
  expect_identical(a$status[a$input == "OLD"], "unrecognized_entrez")
  expect_true(all(a$status[a$input %in% c("UNKNOWN", "stat1")] == "unknown_symbol"))
  expect_equal(sum(a$duplicate_input), 1)
  expect_equal(sum(a$duplicate_entrez), 2)
  expect_setequal(unique(a$index), seq_along(ids))
  expect_equal(nrow(subset(z$exclusions, set == "genes")), 5)
  expect_identical(z$results$GO_BP$status, "no_testable_terms")
  expect_null(z$results$GO_BP$effective_universe)
})

test_that("Entrez recognition, duplicates and empty selections are explicit", {
  local_mocked_bindings(.enrichment_annotation = function(ids, id_type) {
    list(pairs = data.frame(input = c("1", "2"), entrez_id = c("1", "2")),
      recognized = c("1", "2"), metadata = data.frame())
  }, .enrichment_run = function(...) NULL, .package = "TensorPLS")
  ids <- c("1", "1", "999", paste(rep("9", 200), collapse = ""))
  z <- tensorpls_enrichment(ids, c(ids, "2"), databases = "GO_BP")
  expect_identical(z$converted$genes, "1")
  expect_identical(z$converted$universe, c("1", "2"))
  expect_identical(z$mapping$input[z$mapping$set == "genes"], ids)
  expect_true(all(z$exclusions$status == "unrecognized_entrez"))
  expect_equal(nrow(z$exclusions), 4)
  local_mocked_bindings(.enrichment_run = function(...) stop("backend called"), .package = "TensorPLS")
  for (ids in list(character(), "999")) {
    z <- tensorpls_enrichment(ids, ids, databases = "GO_BP")
    expect_identical(z$results$GO_BP$status, "no_recognized_genes")
    expect_equal(nrow(z$results$GO_BP$full), 0)
  }
})

test_that("converted subset consistency is independently checked", {
  calls <- 0L
  local_mocked_bindings(.enrichment_annotation = function(...) .enrichment_fixture(),
    .enrichment_map = function(ids, set, ...) {
      data.frame(set, index = 1L, input = "STAT1", entrez_id = if(set == "genes") "1" else "2",
        status = "recognized", duplicate_input = FALSE, duplicate_entrez = FALSE)
    }, .enrichment_run = function(...) stop("backend called"), .package = "TensorPLS")
  expect_error(tensorpls_enrichment("STAT1", "STAT1", "SYMBOL", "GO_BP"), "converted universe")
})

test_that("frozen Hallmark annotation is required and validated", {
  h <- data.frame(term = "A", gene = "1")
  for (bad in list(NULL, data.frame(), data.frame(x = "A", y = "1"),
    data.frame(term = "A", gene = 1), data.frame(term = NA_character_, gene = "1"),
    data.frame(term = "", gene = "1"), data.frame(term = "A", gene = "01"), rbind(h, h))) {
    expect_error(tensorpls_enrichment("1", "1", databases = "Hallmark", hallmark_term2gene = bad), "TERM2GENE")
  }
  expect_error(tensorpls_enrichment("1", "1", databases = "GO_BP", hallmark_version = NA_character_), "hallmark_version")
})

test_that("backend errors are not reported as nonsignificant results", {
  local_mocked_bindings(.enrichment_annotation = function(...) .enrichment_fixture(),
    .enrichment_run = function(...) stop("database failure"), .package = "TensorPLS")
  expect_error(tensorpls_enrichment("STAT1", "STAT1", "SYMBOL", "GO_BP"), "database failure")
})
