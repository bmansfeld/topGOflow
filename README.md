# topGOflow

**Streamlined GO term enrichment for non-model organisms.**

`topGOflow` wraps the [`topGO`](https://bioconductor.org/packages/topGO/) package into a clean, minimal interface designed for labs working on non-model species. It handles GAF file parsing, expression-matched background selection (via a DESeq2 `dds` object), and returns tidy, filtered result tables — without requiring a Bioconductor `OrgDb`. Optionally builds an OrgDb from a GAF file for downstream semantic similarity analysis with `rrvgo`.

---

## Why this package?

- **topGO is topology-aware.** The `weight01` algorithm accounts for GO graph structure so that specific terms aren't inflated just because a parent term is also enriched. This is the right test to use.
- **Background matching matters.** For gene module or DE gene sets from RNA-seq, highly expressed genes are more likely to be co-expressed *and* better-annotated. If your background is all genes, you're testing expression level as much as function. `topGOflow` defaults to matching the background set to your query by expression level.
- **No OrgDb required.** You just need a GAF file (standard output from most genome annotation pipelines) or any two-column gene → GO mapping.

---

## Installation

```r
# Install dependencies first (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("topGO", "DESeq2", "GO.db"))
install.packages(c("MatchIt", "dplyr", "tidyr", "readr", "ggplot2"))

# Optional: for OrgDb building and semantic similarity reduction
BiocManager::install(c("AnnotationForge", "rrvgo"))

# Install topGOflow from source
# (once on GitHub, replace with remotes::install_github("yourlab/topGOflow"))
install.packages("path/to/topGOflow_0.1.3.tar.gz", repos = NULL, type = "source")
```

---

## Quick start

```r
library(topGOflow)

# ── Step 1: Load your GO annotations from a GAF file ────────────────────────
GOdb <- read_gaf("Solanum_lycopersicum.SGN.CLEANED.gaf")
# → Named list: each element is a character vector of GO terms for that gene

# ── Step 2: Define your query genes ─────────────────────────────────────────
# Any gene list — a WGCNA module, DE genes, a candidate set, etc.
my_genes <- c("Solyc01g005000", "Solyc01g006550", ...)

# ── Step 3: Run enrichment ───────────────────────────────────────────────────
results <- run_go_enrichment(
  query_genes = my_genes,
  GOdb        = GOdb,
  ontologies  = "BP",       # "BP", "MF", "CC", or a vector of all three
  background  = "matched",  # recommended for RNA-seq; see below
  dds         = dds,
  n_matched   = 5,
  node_size   = 10,
  pval        = 0.05,       # only terms below this are returned
  add_genes   = TRUE,       # include query_genes_in_term column
  plot        = TRUE        # show expression matching QC plot
)

# ── Step 4: Results are already filtered and annotated ───────────────────────
results$BP
#   GO.ID      Term               Annotated Significant Expected Fisher.weight01 query_genes_in_term
#   GO:0006355 transcription reg. 312       18          4.2      0.0003          Solyc01g005000; Solyc03g...
```

---

## Background options

| `background =` | What it does | When to use |
|---|---|---|
| `"matched"` *(default)* | Selects background genes from `dds` matched to query by mean normalised expression | **RNA-seq gene sets** (modules, DE genes). Recommended. |
| `"full"` | All genes in `dds` serve as background | Quick check, or when you have no expression bias concern |
| A character vector | You provide the background gene IDs explicitly | Any situation without a `dds`, e.g. proteomics, custom subsets |

---

## Detailed usage

### Reading a GAF file

```r
GOdb <- read_gaf(
  gaf_file         = "myspecies.gaf",
  skip             = 2,      # standard for GAF v2
  strip_suffix     = TRUE,   # "Gene.1.1" → "Gene.1" (set FALSE to keep as-is)
  evidence_exclude = "IEA"   # exclude electronic annotations (default)
                             # set NULL to keep all evidence codes
)
```

The result is a named list — the same format topGO expects.

### Running enrichment without a dds object

```r
# User-supplied background (e.g. all annotated genes, or a proteome set)
results <- run_go_enrichment(
  query_genes = my_genes,
  GOdb        = GOdb,
  ontologies  = c("BP", "MF", "CC"),
  background  = names(GOdb)   # explicit character vector; no dds needed
)
```

### The results table

Each element of the returned list is a data frame filtered to `Fisher.weight01 < pval`, with columns:

| Column | Description |
|---|---|
| `GO.ID` | GO term accession |
| `Term` | Human-readable GO term name |
| `Annotated` | Genes annotated to this term in the universe |
| `Significant` | Query genes annotated to this term |
| `Expected` | Expected count under the null |
| `Fisher.weight01` | Topology-aware p-value |
| `query_genes_in_term` | Semicolon-separated query genes driving this term (`add_genes = TRUE`) |

To re-filter at a different threshold after the fact, use `filter_significant()`:

```r
# Combine ontologies and apply a stricter threshold
filter_significant(results, pval = 0.01, min_sig = 3)
```

### Expression matching QC

Always check the matching worked:

```r
plot_matching(query_genes = my_genes, bg_genes = my_background, dds = dds)
```

The query (pink) and background (blue) density curves should overlap well. If they don't, try increasing `n_matched`.

---

## Building an OrgDb (optional)

An OrgDb enables downstream semantic similarity analysis with `rrvgo`. Build one from the same GAF file:

```r
build_orgdb(
  gaf_file   = "tomato.gaf",
  genus      = "Solanum",
  species    = "lycopersicum",
  tax_id     = "4081",
  maintainer = "Ben Mansfeld <bmansfeld@wustl.edu>",
  author     = "Ben Mansfeld <bmansfeld@wustl.edu>",
  output_dir = ".",
  gene_desc  = my_gene_descriptions,  # optional: data frame with geneid + description columns
  overwrite  = TRUE
)

install.packages("org.Slycopersicum.eg.db", repos = NULL, type = "source")
library(org.Slycopersicum.eg.db)
```

---

## Semantic similarity reduction with rrvgo

Once you have an OrgDb, reduce redundant GO terms and visualise clusters:

```r
# Reduce: collapses similar terms, keeps most significant as representative
results_reduced <- reduce_go_terms(
  go_result = results$BP,
  orgdb     = "org.Slycopersicum.eg.db",
  ontology  = "BP",
  threshold = 0.9    # lower = keep more terms
)

# Results now have parentTerm and score columns
results_reduced |> dplyr::select(GO.ID, Term, Fisher.weight01, parentTerm)

# Heatmap (requires the sim matrix — see ?plot_reduced_go for the full pattern)
plot_reduced_go(results_reduced, sim_matrix)
```

---

## Full WGCNA module workflow example

```r
library(topGOflow)
library(dplyr)

# Load GO database once per session
GOdb <- read_gaf("tomato.gaf")

# Get high-confidence hub genes from module 1
mod1_genes <- moduleExpSpen |>
  left_join(tidyGeneModuleMembership, by = c("gene", "moduleLabel")) |>
  filter(moduleLabel == 1, kME > 0.90) |>
  pull(gene) |>
  unique()

# Run enrichment — results come back filtered and annotated
mod1_results <- run_go_enrichment(
  query_genes = mod1_genes,
  GOdb        = GOdb,
  ontologies  = "BP",
  background  = "matched",
  dds         = dds,
  n_matched   = 5,
  node_size   = 20,
  pval        = 0.05,
  plot        = TRUE
)

mod1_results$BP

# Optionally reduce with rrvgo
mod1_reduced <- reduce_go_terms(mod1_results$BP, orgdb = "org.Slycopersicum.eg.db")
```

---

## Function reference

| Function | Description |
|---|---|
| `read_gaf()` | Parse a GAF file into a gene-to-GO list |
| `run_go_enrichment()` | Main enrichment function — matched or unmatched, returns filtered table |
| `plot_matching()` | QC density plot for expression matching |
| `filter_significant()` | Re-filter result list at a different threshold or combine ontologies |
| `build_orgdb()` | Build a species OrgDb package from a GAF file |
| `reduce_go_terms()` | Semantic similarity reduction via rrvgo |
| `plot_reduced_go()` | Heatmap of reduced GO terms |

---

## Notes on topGO's `weight01` algorithm

This package uses `weight01` exclusively. This is deliberate:

- **Classic Fisher** treats all GO terms independently and inflates p-values for parent terms that are only enriched because a child term is.
- **`weight01`** upweights specific terms and downweights general ones, reflecting the actual structure of the GO DAG.
- The trade-off is that `weight01` p-values are not correctable with standard FDR methods (they're not independent tests). The conventional threshold is **p < 0.05** used as a filter, not as a corrected value.

---

## Citation

If you use `topGOflow`, please also cite the underlying tools:

- **topGO**: Alexa A, Rahnenführer J (2009). *topGO: Enrichment Analysis for Gene Ontology*. R package, Bioconductor.
- **MatchIt**: Ho DE et al. (2011). *MatchIt: Nonparametric Preprocessing for Parametric Causal Inference*. J. Statistical Software.
- **DESeq2**: Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology.
- **rrvgo**: Sayols S (2023). *rrvgo: a Bioconductor package for interpreting lists of Gene Ontology terms*. microPublication Biology.
- **AnnotationForge**: Carlson M, Pages H (2023). *AnnotationForge: Code for building annotation database packages*. Bioconductor.