################################################################################
## topGOflow — OrgDb building and rrvgo semantic similarity reduction
##
## build_orgdb()       – build a species OrgDb package from a GAF file
## reduce_go_terms()   – semantic similarity reduction via rrvgo
## plot_reduced_go()   – heatmap of reduced GO terms
################################################################################

#' @importFrom dplyr transmute filter distinct left_join mutate
#' @importFrom tidyr replace_na
#' @importFrom rlang check_installed
NULL


# ── OrgDb BUILDER ─────────────────────────────────────────────────────────────

#' Build a species-specific OrgDb package from a GAF file
#'
#' Parses a GAF file and calls [AnnotationForge::makeOrgPackage()] to produce
#' an installable OrgDb R package for use with packages like `rrvgo`, `clusterProfiler`,
#' and `enrichplot`. The resulting package uses `"GID"` as the key type, matching
#' the gene IDs produced by [read_gaf()].
#'
#' Once built, install the package with:
#' ```r
#' install.packages("org.Slycopersicum.eg.db", repos = NULL, type = "source")
#' ```
#'
#' @param gaf_file Path to the GAF file.
#' @param genus Genus name, e.g. `"Solanum"`.
#' @param species Species epithet, e.g. `"lycopersicum"`.
#' @param tax_id NCBI taxonomy ID as a string, e.g. `"4081"` for tomato.
#' @param maintainer Maintainer string in `"Name <email>"` format.
#' @param author Author string in `"Name <email>"` format.
#' @param output_dir Directory where the package will be written. Default `"."`.
#' @param version Package version string. Default `"0.1"`.
#' @param skip Number of GAF header lines to skip. Default 2.
#' @param strip_suffix Logical. Strip isoform suffixes from gene IDs. Default TRUE.
#' @param evidence_exclude Character vector of evidence codes to exclude. Default `"IEA"`.
#' @param gene_desc Optional data frame with columns `geneid` and `description`
#'   to populate the `GENENAME` field. If NULL, all genes get `"No_Annotation"`.
#' @param overwrite Logical. If TRUE, delete any existing package directory
#'   before building. Default FALSE (AnnotationForge will error if it exists).
#'
#' @return Invisibly returns the path to the built package directory.
#'
#' @examples
#' \dontrun{
#' build_orgdb(
#'   gaf_file   = "tomato.gaf",
#'   genus      = "Solanum",
#'   species    = "lycopersicum",
#'   tax_id     = "4081",
#'   maintainer = "Ben Mansfeld <bmansfeld@wustl.edu>",
#'   author     = "Ben Mansfeld <bmansfeld@wustl.edu>",
#'   output_dir = ".",
#'   gene_desc  = my_gene_descriptions   # optional
#' )
#' install.packages("org.Slycopersicum.eg.db", repos = NULL, type = "source")
#' library(org.Slycopersicum.eg.db)
#' }
#'
#' @export
build_orgdb <- function(gaf_file,
                        genus,
                        species,
                        tax_id,
                        maintainer,
                        author,
                        output_dir      = ".",
                        version         = "0.1",
                        skip            = 2,
                        strip_suffix    = TRUE,
                        evidence_exclude = "IEA",
                        gene_desc       = NULL,
                        overwrite       = FALSE) {

  if (!requireNamespace("AnnotationForge", quietly = TRUE)) {
    stop("Package 'AnnotationForge' is required. Install with:\n",
         "  BiocManager::install('AnnotationForge')")
  }

  aspect_to_ont <- c(P = "BP", F = "MF", C = "CC")

  gaf_cols <- c(
    "db", "db_object_id", "db_object_symbol", "qualifier",
    "term_accession", "db_reference", "evidence_code", "with",
    "aspect", "db_object_name", "db_object_synonym", "db_object_type",
    "taxon", "date", "assigned_by", "annotation_extension",
    "gene_product_form_id"
  )

  message("[topGOflow] Reading GAF file ...")
  gaf <- readr::read_tsv(
    gaf_file,
    skip      = skip,
    col_names = gaf_cols,
    col_types = readr::cols(.default = "c"),
    comment   = "!"
  )

  # Filter NOT qualifiers
  gaf_clean <- gaf[is.na(gaf$qualifier) | gaf$qualifier == "" |
                     !grepl("(^|\\|)NOT(\\||$)", gaf$qualifier), ]

  if (!is.null(evidence_exclude) && length(evidence_exclude) > 0) {
    gaf_clean <- gaf_clean[!gaf_clean$evidence_code %in% evidence_exclude, ]
  }

  if (strip_suffix) {
    gaf_clean$gene_id <- sub("\\..*", "", gaf_clean$db_object_symbol)
  } else {
    gaf_clean$gene_id <- gaf_clean$db_object_symbol
  }

  # ── GO table ----------------------------------------------------------------
  message("[topGOflow] Building GO annotation table ...")
  go <- gaf_clean |>
    dplyr::transmute(
      GID      = gene_id,
      GO       = term_accession,
      EVIDENCE = evidence_code,
      ONT      = unname(aspect_to_ont[aspect])
    ) |>
    dplyr::filter(
      !is.na(GID), GID != "",
      !is.na(GO), grepl("^GO:\\d{7}$", GO),
      !is.na(ONT)
    ) |>
    dplyr::distinct() |>
    dplyr::select(GID, GO, EVIDENCE)

  # ── gene info table ---------------------------------------------------------
  message("[topGOflow] Building gene info table ...")
  gene_info <- gaf_clean |>
    dplyr::transmute(GID = gene_id, SYMBOL = gene_id) |>
    dplyr::distinct()

  if (!is.null(gene_desc)) {
    gene_info <- gene_info |>
      dplyr::left_join(gene_desc, by = c("GID" = "geneid")) |>
      dplyr::mutate(GENENAME = tidyr::replace_na(description, "No_Annotation")) |>
      dplyr::select(GID, SYMBOL, GENENAME)
  } else {
    gene_info$GENENAME <- "No_Annotation"
    gene_info <- gene_info[, c("GID", "SYMBOL", "GENENAME")]
  }

  gene_info <- as.data.frame(gene_info)
  go        <- as.data.frame(go)

  # ── optionally clear existing package dir ----------------------------------
  pkg_name <- sprintf("org.%s%s.eg.db",
                      substring(genus, 1, 1),
                      species)
  pkg_path <- file.path(output_dir, pkg_name)

  if (overwrite && dir.exists(pkg_path)) {
    message(sprintf("[topGOflow] Removing existing package at %s ...", pkg_path))
    unlink(pkg_path, recursive = TRUE)
  }

  # ── build ------------------------------------------------------------------
  message(sprintf("[topGOflow] Building OrgDb package '%s' in %s ...", pkg_name, output_dir))
  AnnotationForge::makeOrgPackage(
    gene_info  = gene_info,
    go         = go,
    version    = version,
    maintainer = maintainer,
    author     = author,
    outputDir  = output_dir,
    tax_id     = tax_id,
    genus      = genus,
    species    = species,
    goTable    = "go"
  )

  message(sprintf(
    "[topGOflow] Done. Install with:\n  install.packages('%s', repos = NULL, type = 'source')",
    pkg_path
  ))

  invisible(pkg_path)
}


# ── RRVGO SEMANTIC REDUCTION ──────────────────────────────────────────────────

#' Reduce GO enrichment results using semantic similarity (rrvgo)
#'
#' Wraps [rrvgo::calculateSimMatrix()] and [rrvgo::reduceSimMatrix()] to
#' collapse redundant GO terms in the output of [run_go_enrichment()]. Requires
#' an installed OrgDb package (built by [build_orgdb()] or from Bioconductor).
#'
#' The score used for reduction is `-log10(Fisher.weight01)`, so more significant
#' terms are preferentially kept as cluster representatives.
#'
#' @param go_result A data frame — one element of the list returned by
#'   [run_go_enrichment()], e.g. `results$BP`. Must contain columns `GO.ID`
#'   and `Fisher.weight01`.
#' @param orgdb The OrgDb object or package name string, e.g.
#'   `"org.Slycopersicum.eg.db"`. Must be installed and loadable.
#' @param ontology One of `"BP"`, `"MF"`, or `"CC"`. Must match the ontology
#'   used to generate `go_result`. Default `"BP"`.
#' @param keytype Key type for the OrgDb. Default `"GID"` (matches [build_orgdb()]
#'   output). Use `"ENTREZID"` for standard Bioconductor OrgDbs.
#' @param method Similarity measure passed to [rrvgo::calculateSimMatrix()].
#'   One of `"Rel"`, `"Lin"`, `"Resnik"`, `"Jiang"`. Default `"Rel"`.
#' @param threshold Similarity threshold for [rrvgo::reduceSimMatrix()].
#'   Terms with similarity above this are collapsed. Default `0.9` (aggressive
#'   reduction); lower values (e.g. `0.7`) keep more terms.
#' @param pval_filter Pre-filter `go_result` to terms with p-value below this
#'   before reduction. Default `0.05`.
#'
#' @return The input `go_result` data frame with two additional columns:
#'   \describe{
#'     \item{parentTerm}{The representative term for each cluster.}
#'     \item{score}{The `-log10(Fisher.weight01)` score used for reduction.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(org.Slycopersicum.eg.db)
#'
#' results <- run_go_enrichment(my_genes, GOdb, ontologies = "BP", dds = dds)
#'
#' results_reduced <- reduce_go_terms(
#'   go_result = results$BP,
#'   orgdb     = "org.Slycopersicum.eg.db",
#'   ontology  = "BP",
#'   threshold = 0.9
#' )
#'
#' # Plot
#' plot_reduced_go(results_reduced)
#' }
#'
#' @export
reduce_go_terms <- function(go_result,
                            orgdb,
                            ontology    = "BP",
                            keytype     = "GID",
                            method      = "Rel",
                            threshold   = 0.9,
                            pval_filter = 0.05) {

  if (!requireNamespace("rrvgo", quietly = TRUE)) {
    stop("Package 'rrvgo' is required. Install with:\n",
         "  BiocManager::install('rrvgo')")
  }

  ontology <- match.arg(ontology, c("BP", "MF", "CC"))
  method   <- match.arg(method, c("Rel", "Lin", "Resnik", "Jiang"))

  # Filter to significant terms
  go_sig <- go_result[go_result$Fisher.weight01 < pval_filter, ]

  if (nrow(go_sig) == 0) {
    warning("[topGOflow] No terms pass the p-value filter. Returning input unchanged.")
    return(go_result)
  }

  message(sprintf("[topGOflow] Calculating similarity matrix for %d GO terms ...", nrow(go_sig)))

  sim_matrix <- rrvgo::calculateSimMatrix(
    go_sig$GO.ID,
    orgdb   = orgdb,
    keytype = keytype,
    ont     = ontology,
    method  = method
  )

  scores <- setNames(-log10(go_sig$Fisher.weight01), go_sig$GO.ID)

  message("[topGOflow] Reducing similar terms ...")
  reduced <- rrvgo::reduceSimMatrix(
    sim_matrix,
    scores    = scores,
    threshold = threshold,
    orgdb     = orgdb,
    keytype   = keytype
  )

  # Join parentTerm and score back onto the original filtered table
  go_sig <- merge(
    go_sig,
    reduced[, c("go", "parentTerm", "score")],
    by.x = "GO.ID",
    by.y = "go",
    all.x = TRUE
  )

  go_sig[order(go_sig$Fisher.weight01), ]
}


# ── RRVGO PLOTTING ────────────────────────────────────────────────────────────

#' Heatmap of semantically reduced GO terms
#'
#' Convenience wrapper around [rrvgo::heatmapPlot()]. Pass the output of
#' [reduce_go_terms()] along with the original similarity matrix.
#'
#' Because [rrvgo::calculateSimMatrix()] needs to be called to produce both
#' the reduction and the plot, the recommended pattern is to call
#' [reduce_go_terms()] first, then recompute the sim matrix for plotting —
#' or save it from an intermediate step. See the example below.
#'
#' @param go_result_reduced Data frame returned by [reduce_go_terms()].
#' @param sim_matrix Similarity matrix from [rrvgo::calculateSimMatrix()].
#' @param fontsize Font size for term labels. Default 6.
#' @param ... Additional arguments passed to [rrvgo::heatmapPlot()].
#'
#' @return Called for its side effect (prints a plot). Invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#' # Compute sim matrix once, use for both reduction and plot
#' go_sig <- results$BP[results$BP$Fisher.weight01 < 0.05, ]
#'
#' sim_matrix <- rrvgo::calculateSimMatrix(
#'   go_sig$GO.ID,
#'   orgdb   = "org.Slycopersicum.eg.db",
#'   keytype = "GID",
#'   ont     = "BP",
#'   method  = "Rel"
#' )
#' scores  <- setNames(-log10(go_sig$Fisher.weight01), go_sig$GO.ID)
#' reduced <- rrvgo::reduceSimMatrix(sim_matrix, scores,
#'                                   threshold = 0.9,
#'                                   orgdb = "org.Slycopersicum.eg.db",
#'                                   keytype = "GID")
#'
#' rrvgo::heatmapPlot(sim_matrix, reduced,
#'                    annotateParent  = TRUE,
#'                    annotationLabel = "parentTerm",
#'                    fontsize        = 6)
#' }
#'
#' @export
plot_reduced_go <- function(go_result_reduced,
                            sim_matrix,
                            fontsize = 6,
                            ...) {

  if (!requireNamespace("rrvgo", quietly = TRUE)) {
    stop("Package 'rrvgo' is required. Install with:\n",
         "  BiocManager::install('rrvgo')")
  }

  if (!all(c("parentTerm") %in% colnames(go_result_reduced))) {
    stop("`go_result_reduced` must contain a 'parentTerm' column. ",
         "Run reduce_go_terms() first.")
  }

  reduced_df <- data.frame(
    go         = go_result_reduced$GO.ID,
    parentTerm = go_result_reduced$parentTerm,
    score      = go_result_reduced$score
  )

  rrvgo::heatmapPlot(
    sim_matrix,
    reduced_df,
    annotateParent  = TRUE,
    annotationLabel = "parentTerm",
    fontsize        = fontsize,
    ...
  )

  invisible(NULL)
}
