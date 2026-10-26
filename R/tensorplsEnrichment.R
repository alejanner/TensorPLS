#' Gene-based functional enrichment
#'
#' Run over-representation analysis for human GO Biological Process, Reactome,
#' and Hallmark gene sets.
#'
#' @param genes Selected gene identifiers.
#' @param universe Gene identifiers eligible for selection, including `genes`.
#' @param id_type `"ENTREZID"` or `"SYMBOL"`.
#' @param databases One or more of `"GO_BP"`, `"Reactome"`, and `"Hallmark"`.
#' @param organism Currently `"human"`.
#' @param hallmark_term2gene Data frame with character columns `term` and
#'   `gene` (Entrez IDs), required for Hallmark.
#' @param hallmark_version Optional annotation release label.
#' @return A list containing enrichment results, gene mappings, exclusions,
#'   parameters, and provenance.
#' @export
tensorpls_enrichment <- function(genes, universe, id_type = "ENTREZID",
  databases = c("GO_BP", "Reactome", "Hallmark"), organism = "human",
  hallmark_term2gene = NULL, hallmark_version = NULL) {
  if (!identical(organism, "human")) stop("Only organism='human' is supported.")
  if (length(id_type) != 1L || is.na(id_type) || !id_type %in% c("ENTREZID", "SYMBOL"))
    stop("id_type must be ENTREZID or SYMBOL.")
  for (x in list(genes = genes, universe = universe)) {
    if (!is.character(x) || anyNA(x) || any(!nzchar(x)))
      stop("genes and universe must be character vectors without NA or empty identifiers.")
    if (id_type == "ENTREZID" && any(!grepl("^[1-9][0-9]*$", x)))
      stop("Malformed ENTREZID: use positive digit strings without leading zeros or whitespace.")
  }
  if (!all(genes %in% universe)) stop("Selected genes must belong to the supplied universe before conversion.")
  allowed <- c("GO_BP", "Reactome", "Hallmark")
  if (
    !is.character(databases) ||
    !length(databases) ||
    anyNA(databases) ||
    anyDuplicated(databases) ||
    !all(databases %in% allowed)
  )
    stop("databases must be a nonempty, unique subset of GO_BP, Reactome, Hallmark.")
  if (!is.null(hallmark_version) && (!is.character(hallmark_version) ||
      length(hallmark_version) != 1L || is.na(hallmark_version) || !nzchar(hallmark_version)))
    stop("hallmark_version must be NULL or one nonempty release label.")
  if ("Hallmark" %in% databases) {
    h <- hallmark_term2gene
    if (
      !is.data.frame(h) ||
      !identical(names(h), c("term", "gene")) ||
      !nrow(h) ||
      !is.character(h$term) ||
      !is.character(h$gene) ||
      anyNA(h) ||
      any(!nzchar(h$term)) ||
      any(!grepl("^[1-9][0-9]*$", h$gene)) ||
      anyDuplicated(h)
    )
      stop("Hallmark requires a nonempty TERM2GENE data frame with unique character term/gene pairs and valid Entrez syntax.")
  }
  annotation <- .enrichment_annotation(unique(c(genes, universe)), id_type)
  mapping <- rbind(.enrichment_map(genes, "genes", annotation, id_type),
                   .enrichment_map(universe, "universe", annotation, id_type))
  converted <- lapply(c("genes", "universe"), function(set) {
    unique(mapping$entrez_id[mapping$set == set & mapping$status == "recognized"])
  })
  names(converted) <- c("genes", "universe")
  if (!all(converted$genes %in% converted$universe)) stop("Selected genes are outside the converted universe.")
  parameters <- list(id_type = id_type, databases = databases, organism = organism,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 500, significance = "p.adjust < 0.05")
  results <- lapply(databases, function(database) {
    native <- if (length(converted$genes)) {
      .enrichment_run(database, converted$genes, converted$universe, hallmark_term2gene)
    } else NULL
    if (is.null(native)) {
      return(list(status = if (length(converted$genes)) "no_testable_terms" else "no_recognized_genes",
        full = data.frame(), significant = data.frame(), effective_universe = NULL,
        effective_genes = NULL, native = NULL))
    }
    full <- native@result
    effective <- native@universe
    selected <- intersect(converted$genes, effective)
    if (!all(effective %in% converted$universe)) stop("Backend universe exceeds the submitted universe.")
    significant <- full[!is.na(full$p.adjust) & full$p.adjust < 0.05, , drop = FALSE]
    list(status = if (nrow(full)) "tested" else "no_testable_terms", full = full,
      significant = significant, effective_universe = effective,
      effective_genes = selected, native = native)
  })
  names(results) <- databases
  packages <- c("TensorPLS", "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler", "DOSE",
    if ("GO_BP" %in% databases) "GO.db", if ("Reactome" %in% databases) c("ReactomePA", "reactome.db"))
  versions <- vapply(packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
  list(results = results, original = list(genes = genes, universe = universe),
    converted = converted, mapping = mapping,
    exclusions = mapping[mapping$status != "recognized", , drop = FALSE],
    parameters = parameters, provenance = list(packages = versions,
      orgdb_metadata = annotation$metadata, R = R.version.string,
      hallmark_version = hallmark_version, hallmark_term2gene = hallmark_term2gene))
}

.enrichment_annotation <- function(ids, id_type) {
  for (p in c("AnnotationDbi", "org.Hs.eg.db")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Install optional enrichment dependency: ", p)
  }
  org <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  recognized <- AnnotationDbi::keys(org, keytype = "ENTREZID")
  pairs <- data.frame(input = character(), entrez_id = character())
  if (id_type == "SYMBOL") {
    keys <- intersect(ids, AnnotationDbi::keys(org, keytype = "SYMBOL"))
    if (length(keys)) {
      selected <- AnnotationDbi::select(org, keys = keys, keytype = "SYMBOL", columns = "ENTREZID")
      pairs <- unique(data.frame(input = selected$SYMBOL, entrez_id = selected$ENTREZID))
    }
  } else {
    keys <- intersect(ids, recognized)
    pairs <- data.frame(input = keys, entrez_id = keys)
  }
  list(pairs = pairs, recognized = recognized, metadata = AnnotationDbi::metadata(org))
}

.enrichment_map <- function(ids, set, annotation, id_type) {
  empty <- data.frame(set = character(), index = integer(), input = character(),
    entrez_id = character(), status = character(), duplicate_input = logical(),
    duplicate_entrez = logical())
  if (!length(ids)) return(empty)
  duplicates <- duplicated(ids)
  lookup <- split(annotation$pairs$entrez_id, annotation$pairs$input)
  rows <- lapply(seq_along(ids), function(i) {
    candidates <- unique(lookup[[ids[i]]])
    candidates <- candidates[!is.na(candidates)]
    status <- if (!length(candidates)) {
      if (id_type == "SYMBOL") "unknown_symbol" else "unrecognized_entrez"
    } else if (length(candidates) > 1L) "ambiguous_symbol" else {
      if (candidates %in% annotation$recognized) "recognized" else "unrecognized_entrez"
    }
    data.frame(set = set, index = i, input = ids[i],
      entrez_id = if (length(candidates)) candidates else NA_character_, status = status,
      duplicate_input = duplicates[i], duplicate_entrez = FALSE)
  })
  result <- do.call(rbind, rows)
  accepted <- which(result$status == "recognized")
  result$duplicate_entrez[accepted] <- duplicated(result$entrez_id[accepted])
  result
}

.enrichment_run <- function(database, genes, universe, hallmark) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install optional enrichment dependency: clusterProfiler")
  arguments <- list(gene = genes, universe = universe, pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
  if (database == "GO_BP") {
    arguments$OrgDb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    arguments$keyType <- "ENTREZID"
    arguments$ont <- "BP"
    arguments$readable <- FALSE
    return(do.call(clusterProfiler::enrichGO, arguments))
  }
  if (database == "Reactome") {
    if (!requireNamespace("ReactomePA", quietly = TRUE))
      stop("Install optional enrichment dependency: ReactomePA")
    arguments$organism <- "human"
    arguments$readable <- FALSE
    return(do.call(ReactomePA::enrichPathway, arguments))
  }
  arguments$TERM2GENE <- hallmark
  do.call(clusterProfiler::enricher, arguments)
}
