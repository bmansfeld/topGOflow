# topGOflow

**Streamlined GO term enrichment for non-model organisms.**

`topGOflow` wraps the [`topGO`](https://bioconductor.org/packages/topGO/) package into a clean, minimal interface designed for labs working on non-model species. It handles GAF file parsing, expression-matched background selection (via a DESeq2 `dds` object), and returns tidy result tables — without requiring a Bioconductor `OrgDb`.

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
install.packages(c("MatchIt", "dplyr", "readr", "ggplot2"))

# Install topGOflow from source
remotes::install_github("bmansfeld/topGOflow")
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
  ontologies  = "BP",           # "BP", "MF", "CC", or a vector of all three
  background  = "matched",      # recommended for RNA-seq; see below for options
  dds         = dds,            # your DESeqDataSet
  n_matched   = 5,              # 5 background genes per query gene
  node_size   = 10,             # ignore GO terms with < 10 annotated genes
  plot        = TRUE            # show expression matching QC plot
)

# ── Step 4: Explore results ──────────────────────────────────────────────────
results$BP                                   # full table for Biological Process
filter_significant(results, pval = 0.05)     # significant terms across all ontologies
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
  gaf_file        = "myspecies.gaf",
  skip            = 2,              # standard for GAF v2
  strip_suffix    = TRUE,           # "Gene.1.1" → "Gene.1" (set FALSE to keep)
  evidence_exclude = "IEA"          # exclude electronic annotations (default)
                                    # set NULL to keep all evidence codes
)
```

The result is a named list — the same format topGO expects, and the same as what you built manually before with `setNames(nm = ..., strsplit(...))`.

### Running enrichment without a dds object

```r
# User-supplied background (e.g. all annotated genes, or a proteome set)
all_genes <- names(GOdb)

results <- run_go_enrichment(
  query_genes = my_genes,
  GOdb        = GOdb,
  ontologies  = c("BP", "MF", "CC"),
  background  = all_genes           # explicit character vector; no dds needed
)
```

### Checking which genes drive a term

```r
# First get the topGOdata object for a given ontology
# (run_go_enrichment returns tables; to look up genes you need the tgd object)
# Re-run with background to get the object, or use the internal helper:

sig_genes <- genes_in_term(
  go_term     = "GO:0006355",       # transcription, DNA-templated
  query_genes = my_genes,
  tgd         = my_tgd_object       # topGOdata object
)
```

> **Tip:** If you need to frequently look up genes-in-term, build the `topGOdata` object separately using the internal `.build_topgo_object()` function and keep it alongside your results.

### Filtering results

```r
# All significant BP terms, at least 3 query genes
filter_significant(results, pval = 0.05, ontologies = "BP", min_sig = 3)

# Combine ontologies
filter_significant(results, pval = 0.01)
```

### Expression matching QC

Always check the matching worked:

```r
plot_matching(query_genes = my_genes, bg_genes = my_background, dds = dds)
```

The query (pink) and background (blue) density curves should overlap well. If they don't, try increasing `n_matched`.

---

## Full WGCNA module workflow example

```r
library(topGOflow)
library(dplyr)

# Load GO database once
GOdb <- read_gaf("tomato.gaf")

# Get high-confidence hub genes from module 1
mod1_genes <- moduleExpSpen |>
  left_join(tidyGeneModuleMembership, by = c("gene", "moduleLabel")) |>
  filter(moduleLabel == 1, kME > 0.90) |>
  pull(gene) |>
  unique()

# Run enrichment
mod1_results <- run_go_enrichment(
  query_genes = mod1_genes,
  GOdb        = GOdb,
  ontologies  = "BP",
  background  = "matched",
  dds         = dds,
  n_matched   = 5,
  node_size   = 20,
  plot        = TRUE
)

# Significant terms
mod1_results$BP |> filter(Fisher.weight01 < 0.05)

# Or with the convenience wrapper
filter_significant(mod1_results)
```

---

## Function reference

| Function | Description |
|---|---|
| `read_gaf()` | Parse a GAF file into a gene-to-GO list |
| `run_go_enrichment()` | Main enrichment function — matched or unmatched |
| `plot_matching()` | QC density plot for expression matching |
| `filter_significant()` | Filter result list to significant terms |
| `genes_in_term()` | Which query genes are annotated to a given GO term |

---

## Notes on topGO's `weight01` algorithm

This package uses the `weight01` algorithm exclusively. This is deliberate:

- **Classic Fisher** treats all GO terms independently and inflates p-values for parent terms that are only significant because a specific child term is.
- **`weight01`** upweights specific terms and downweights general ones, reflecting the actual structure of the GO DAG.
- The trade-off is that `weight01` p-values are not correctable with standard FDR methods (they're not independent tests). The conventional threshold is **p < 0.05** used as a filter, not as a corrected value.

This is the same philosophy that has underpinned the lab's GO analyses for the past decade.

---

## Citation

If you use `topGOflow`, please also cite the underlying tools:

- **topGO**: Alexa A, Rahnenführer J (2009). *topGO: Enrichment Analysis for Gene Ontology*. R package, Bioconductor.
- **MatchIt**: Ho DE et al. (2011). *MatchIt: Nonparametric Preprocessing for Parametric Causal Inference*. J. Statistical Software.
- **DESeq2**: Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology.
