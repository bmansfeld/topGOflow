################################################################################
## topGOflow — Streamlined GO Enrichment for Non-Model Organisms
##
## Main entry points:
##   read_gaf()                – parse a GAF file into a gene2GO list
##   run_go_enrichment()       – full matched or unmatched analysis in one call
##   plot_matching()           – QC plot of expression matching
##   go_table()                – tidy result table from a topGOdata object
##   genes_in_term()           – which query genes drive a given GO term
##
## Internal helpers (not exported):
##   .match_background()       – MatchIt expression matching
##   .build_topgo_object()     – construct topGOdata for one ontology
##   .run_tests()              – run weight01 Fisher test and call GenTable
################################################################################

#' @importFrom dplyr select mutate group_by summarise filter left_join pull bind_rows
#' @importFrom tidyr %>%
#' @importFrom readr read_tsv
#' @importFrom ggplot2 ggplot aes geom_density scale_x_log10 labs theme_bw
#' @importFrom methods new
#' @importFrom MatchIt matchit
#' @importFrom DESeq2 counts
#' @importFrom topGO runTest GenTable genesInTerm annFUN.gene2GO
NULL


# ── GAF PARSING ───────────────────────────────────────────────────────────────

#' Read a GAF file and return a gene-to-GO list
#'
#' Parses a Gene Association File (GAF v2.x), extracts the gene symbol and
#' GO accession columns, strips isoform suffixes (everything after the first
#' "."), and returns a named list suitable for topGO's `gene2GO` slot.
#'
#' @param gaf_file Path to the GAF file (plain text or gzipped).
#' @param skip Number of header lines to skip (default 2 for standard GAFs).
#' @param strip_suffix Logical. If TRUE (default), remove everything after the
#'   first "." in gene IDs (e.g. "Solyc01g005000.3" → "Solyc01g005000").
#' @param evidence_exclude Character vector of evidence codes to exclude
#'   (default excludes IEA — Inferred from Electronic Annotation).
#'   Set to `NULL` to keep all evidence codes.
#'
#' @return A named list where each element is a character vector of GO term
#'   accessions for that gene. Ready to pass as `GOdb` to [run_go_enrichment()].
#'
#' @examples
#' \dontrun{
#' GOdb <- read_gaf("Solanum_lycopersicum.gaf")
#' head(GOdb, 3)
#' }
#'
#' @export
read_gaf <- function(gaf_file,
                     skip = 2,
                     strip_suffix = TRUE,
                     evidence_exclude = "IEA") {

  gaf_cols <- c(
    "db", "db_object_id", "db_object_symbol", "qualifier",
    "term_accession", "db_reference", "evidence_code", "with",
    "aspect", "db_object_name", "db_object_synonym", "db_object_type",
    "taxon", "date", "assigned_by", "annotation_extension",
    "gene_product_form_id"
  )

  gaf <- readr::read_tsv(
    gaf_file,
    skip        = skip,
    col_names   = gaf_cols,
    col_types   = readr::cols(.default = "c"),
    comment     = "!"
  )

  if (!is.null(evidence_exclude) && length(evidence_exclude) > 0) {
    gaf <- gaf[!gaf$evidence_code %in% evidence_exclude, ]
  }

  if (strip_suffix) {
    gaf$gene_id <- sub("\\..*", "", gaf$db_object_symbol)
  } else {
    gaf$gene_id <- gaf$db_object_symbol
  }

  gene_go <- gaf |>
    dplyr::select(gene_id, term_accession) |>
    dplyr::group_by(gene_id) |>
    dplyr::summarise(GO_terms = paste(unique(term_accession), collapse = "; "),
                     .groups = "drop")

  go_list <- setNames(
    strsplit(gene_go$GO_terms, "; "),
    gene_go$gene_id
  )

  message(sprintf("[topGOflow] GAF loaded: %d genes with GO annotations.", length(go_list)))
  go_list
}


# ── MAIN USER-FACING FUNCTION ─────────────────────────────────────────────────

#' Run GO enrichment analysis with optional expression-matched background
#'
#' The primary function for GO enrichment. Given a query gene set and a
#' gene-to-GO mapping, it will:
#' \enumerate{
#'   \item Optionally select an expression-matched background from a DESeq2
#'         `dds` object (recommended for RNA-seq gene modules or DE gene sets).
#'   \item Build a `topGOdata` object for each requested ontology.
#'   \item Run the `weight01` Fisher test (topology-aware).
#'   \item Return a tidy named list of result data frames, one per ontology.
#' }
#'
#' @param query_genes Character vector of gene IDs to test (e.g. a module,
#'   DE gene list, or any set of interest).
#' @param GOdb Named list mapping gene IDs to character vectors of GO terms.
#'   Produced by [read_gaf()] or built manually.
#' @param ontologies Character vector of ontologies to test. Any combination of
#'   `"BP"` (Biological Process), `"MF"` (Molecular Function), `"CC"`
#'   (Cellular Component). Default is `"BP"`.
#' @param background One of:
#'   \describe{
#'     \item{`"matched"` (default)}{Expression-matched background selected from
#'       all genes in `dds` using [MatchIt]. Requires `dds` and `n_matched`.}
#'     \item{`"full"`}{All genes in `dds` are used as background (no matching).
#'       Requires `dds`.}
#'     \item{A character vector}{Use this explicit set of gene IDs as background.
#'       Does not require `dds`.}
#'   }
#' @param dds A `DESeqDataSet` object. Required when `background` is `"matched"`
#'   or `"full"`. Used for normalised expression values (via `counts(dds, TRUE)`).
#' @param n_matched Integer. Number of background genes to match per query gene
#'   when `background = "matched"`. Default 5.
#' @param node_size Minimum number of genes required in a GO term for it to be
#'   tested. Default 10. Larger values (e.g. 20–50) reduce noise.
#' @param n_top Maximum number of GO terms to return per ontology. Default 200.
#' @param plot Logical. If TRUE, print a [plot_matching()] density plot when
#'   expression matching is used. Default FALSE.
#'
#' @return A named list with one element per ontology (e.g. `$BP`). Each element
#'   is a data frame with columns:
#'   \describe{
#'     \item{GO.ID}{GO term accession.}
#'     \item{Term}{Human-readable GO term name.}
#'     \item{Annotated}{Number of genes annotated to this term in the universe.}
#'     \item{Significant}{Number of query genes annotated to this term.}
#'     \item{Expected}{Expected count under the null.}
#'     \item{Fisher.weight01}{p-value from the topology-aware weight01 test.}
#'   }
#'
#' @examples
#' \dontrun{
#' GOdb   <- read_gaf("tomato.gaf")
#' result <- run_go_enrichment(
#'   query_genes = my_module_genes,
#'   GOdb        = GOdb,
#'   ontologies  = "BP",
#'   background  = "matched",
#'   dds         = dds,
#'   n_matched   = 5,
#'   node_size   = 20
#' )
#' result$BP |> dplyr::filter(Fisher.weight01 < 0.05)
#' }
#'
#' @export
run_go_enrichment <- function(query_genes,
                              GOdb,
                              ontologies  = "BP",
                              background  = "matched",
                              dds         = NULL,
                              n_matched   = 5,
                              node_size   = 10,
                              n_top       = 200,
                              plot        = FALSE) {

  # ── validate ----------------------------------------------------------------
  ontologies <- match.arg(ontologies, c("BP", "MF", "CC"), several.ok = TRUE)

  query_genes <- unique(as.character(query_genes))
  n_query <- length(query_genes)
  message(sprintf("[topGOflow] Query set: %d genes.", n_query))

  # ── determine background ----------------------------------------------------
  if (identical(background, "matched")) {
    if (is.null(dds)) stop("`dds` required when background = 'matched'.")
    bg_genes <- .match_background(query_genes, dds, nR = n_matched)
    message(sprintf("[topGOflow] Expression-matched background: %d genes (ratio ~%dx).",
                    length(bg_genes), n_matched))
    if (plot) print(plot_matching(query_genes, bg_genes, dds))

  } else if (identical(background, "full")) {
    if (is.null(dds)) stop("`dds` required when background = 'full'.")
    bg_genes <- rownames(dds)
    message(sprintf("[topGOflow] Full background: %d genes.", length(bg_genes)))

  } else if (is.character(background)) {
    bg_genes <- unique(background)
    message(sprintf("[topGOflow] User-supplied background: %d genes.", length(bg_genes)))

  } else {
    stop("`background` must be 'matched', 'full', or a character vector of gene IDs.")
  }

  # ── determine universe gene IDs (for building topGO factor) ----------------
  # Must use rownames(dds) when available — the logical subsetting in
  # .build_topgo_object depends on this. Without it, gene IDs can fail to
  # match the GO database and topGO throws "GOBPTerm not found".
  # When dds is absent, fall back to names(GOdb) so IDs are guaranteed to match.
  if (!is.null(dds)) {
    universe_ids <- rownames(dds)
  } else {
    universe_ids <- names(GOdb)
  }

  # ── run per ontology --------------------------------------------------------
  results <- lapply(ontologies, function(ont) {
    message(sprintf("[topGOflow] Building topGOdata for %s ...", ont))
    tgd <- .build_topgo_object(
      query_genes  = query_genes,
      bg_genes     = bg_genes,
      universe_ids = universe_ids,
      GOdb         = GOdb,
      ontology     = ont,
      node_size    = node_size
    )
    message(sprintf("[topGOflow] Running weight01 Fisher test for %s ...", ont))
    .run_tests(tgd, n_top = n_top)
  })

  names(results) <- ontologies
  results
}


# ── HELPER: EXPRESSION MATCHING ───────────────────────────────────────────────

#' @keywords internal
.match_background <- function(query_genes, dds, nR = 5) {
  mean_exp <- rowMeans(DESeq2::counts(dds, normalized = TRUE))
  df <- data.frame(
    sign    = as.integer(rownames(dds) %in% query_genes),
    geneExp = mean_exp
  )
  match_res <- MatchIt::matchit(
    sign ~ geneExp,
    data     = df,
    method   = "nearest",
    distance = "mahalanobis",
    ratio    = nR
  )
  bg <- as.vector(match_res$match.matrix)
  bg <- unique(stats::na.omit(bg))
  bg
}


# ── HELPER: BUILD topGOdata ───────────────────────────────────────────────────

#' @keywords internal
.build_topgo_object <- function(query_genes,
                                bg_genes,
                                universe_ids,
                                GOdb,
                                ontology  = "BP",
                                node_size = 10) {
  # Mirror the logic from the original makeTopGOobject:
  # both logical vectors are defined over the full universe_ids,
  # then in_query is subset BY in_universe before naming.
  in_universe <- universe_ids %in% c(query_genes, bg_genes)
  in_query    <- universe_ids %in% query_genes

  gene_factor <- factor(as.integer(in_query[in_universe]))
  names(gene_factor) <- universe_ids[in_universe]

  methods::new(
    "topGOdata",
    ontology  = ontology,
    allGenes  = gene_factor,
    nodeSize  = node_size,
    annot     = topGO::annFUN.gene2GO,
    gene2GO   = GOdb
  )
}


# ── HELPER: RUN TESTS & FORMAT ────────────────────────────────────────────────

#' @keywords internal
.run_tests <- function(tgd, n_top = 200) {
  res_w01 <- topGO::runTest(tgd, algorithm = "weight01", statistic = "Fisher")

  n_terms <- min(res_w01@geneData[4], n_top)

  tbl <- topGO::GenTable(
    tgd,
    Fisher.weight01 = res_w01,
    orderBy         = "Fisher.weight01",
    topNodes        = n_terms
  )

  # Convert the truncated "< 1e-30" strings to numeric
  tbl$Fisher.weight01 <- ifelse(
    tbl$Fisher.weight01 == "< 1e-30",
    1e-30,
    as.numeric(tbl$Fisher.weight01)
  )

  tbl
}


# ── QC PLOT ───────────────────────────────────────────────────────────────────

#' Plot expression distributions for QC of background matching
#'
#' Produces a density plot of log10(mean normalised expression) for three gene
#' sets: all genes in `dds`, the query set, and the matched background. Use
#' this to verify that the background distribution is similar to the query.
#'
#' @param query_genes Character vector of query gene IDs.
#' @param bg_genes Character vector of background gene IDs (from matching or
#'   user-supplied).
#' @param dds A `DESeqDataSet` object.
#'
#' @return A `ggplot` object.
#'
#' @export
plot_matching <- function(query_genes, bg_genes, dds) {
  mean_exp <- rowMeans(DESeq2::counts(dds, normalized = TRUE))

  df <- dplyr::bind_rows(
    "All genes"  = data.frame(exp = mean_exp[rownames(dds)]),
    "Query"      = data.frame(exp = mean_exp[query_genes]),
    "Background" = data.frame(exp = mean_exp[bg_genes]),
    .id = "Group"
  )
  df$Group <- factor(df$Group, levels = c("All genes", "Background", "Query"))

  ggplot2::ggplot(df, ggplot2::aes(x = exp, colour = Group, fill = Group)) +
    ggplot2::geom_density(alpha = 0.15, linewidth = 0.8) +
    ggplot2::scale_x_log10(name = "Mean normalised expression (log10)") +
    ggplot2::scale_colour_manual(values = c("All genes" = "grey50",
                                            "Background" = "#2196F3",
                                            "Query"      = "#E91E63")) +
    ggplot2::scale_fill_manual(values = c("All genes" = "grey50",
                                          "Background" = "#2196F3",
                                          "Query"      = "#E91E63")) +
    ggplot2::labs(
      title    = "Expression matching QC",
      subtitle = sprintf("Query: %d genes | Background: %d genes",
                         length(query_genes), length(bg_genes)),
      y        = "Density"
    ) +
    ggplot2::theme_bw(base_size = 12)
}


# ── GENE LOOKUP ───────────────────────────────────────────────────────────────

#' Retrieve query genes annotated to a specific GO term
#'
#' After running [run_go_enrichment()], use this to find which of your query
#' genes are driving enrichment of a particular GO term.
#'
#' @param go_term Character. A GO term accession (e.g. `"GO:0006355"`).
#' @param query_genes Character vector of query gene IDs.
#' @param tgd A `topGOdata` object. Obtained internally by
#'   [run_go_enrichment()] — pass the object returned by `.build_topgo_object()`
#'   or use [get_topgo_object()] to retrieve it.
#'
#' @return A character vector of gene IDs from `query_genes` that are
#'   annotated to `go_term`.
#'
#' @export
genes_in_term <- function(go_term, query_genes, tgd) {
  all_annotated <- unlist(topGO::genesInTerm(tgd, go_term))
  intersect(query_genes, all_annotated)
}


# ── CONVENIENCE: FILTER RESULTS ───────────────────────────────────────────────

#' Filter GO enrichment results to significant terms
#'
#' Convenience wrapper to quickly pull significant GO terms from the list
#' returned by [run_go_enrichment()].
#'
#' @param results Named list returned by [run_go_enrichment()].
#' @param pval Numeric. P-value threshold. Default 0.05.
#' @param ontologies Character vector. Which ontologies to include. Defaults to
#'   all ontologies present in `results`.
#' @param min_sig Integer. Minimum number of query genes in the term. Default 2.
#'
#' @return A single combined data frame with an added `Ontology` column, sorted
#'   by p-value.
#'
#' @export
filter_significant <- function(results,
                               pval       = 0.05,
                               ontologies = NULL,
                               min_sig    = 2) {
  if (is.null(ontologies)) ontologies <- names(results)

  combined <- dplyr::bind_rows(
    lapply(ontologies, function(ont) {
      df <- results[[ont]]
      if (is.null(df)) return(NULL)
      df$Ontology <- ont
      df
    })
  )

  combined |>
    dplyr::filter(Fisher.weight01 < pval, Significant >= min_sig) |>
    dplyr::arrange(Fisher.weight01)
}
