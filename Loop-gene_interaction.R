##############################################################
# Loop filtering and enhancer fragment extraction
#
# Input:
#   --loops     : CSV with loop interactions
#                 cols: chr1,start1,end1,chr2,start2,end2,Name,Genes_F1,Genes_F2
#   --lrt       : LRT_significant_genes_clustered.csv
#   --expr      : geneset_full_classification.csv
#   --outdir    : output directory
#
# Output per cluster/geneset:
#   1) all_loops_<group>.csv     - all loops containing a gene from that group
#   2) enhancer_fragments_<group>.bed - putative enhancer fragments (opposite to gene)
#
# Usage:
#   Rscript loop_filter.R \
#     --loops your_loops.csv \
#     --lrt LRT_significant_genes_clustered.csv \
#     --expr geneset_full_classification.csv \
#     --outdir loop_results
##############################################################

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

option_list <- list(
  make_option("--loops",
              type    = "character",
              default = NULL,
              help    = "Loop interaction CSV file"),
  make_option("--lrt",
              type    = "character",
              default = "LRT_significant_genes_clustered.csv"),
  make_option("--expr",
              type    = "character",
              default = "geneset_full_classification.csv"),
  make_option("--outdir",
              type    = "character",
              default = "loop_results")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$loops)) stop("--loops is required", call. = FALSE)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

##############################################################
# 1. LOAD LOOPS
##############################################################

cat("Loading loop file:", opt$loops, "\n")
loops <- read.csv(opt$loops, stringsAsFactors = FALSE, check.names = FALSE)

# Standardise column names - handle variations
colnames(loops) <- trimws(colnames(loops))

# Expected: chr1,start1,end1,chr2,start2,end2,Name,Genes_F1,Genes_F2
required <- c("chr1","start1","end1","chr2","start2","end2","Name","Genes_F1","Genes_F2")
missing  <- setdiff(required, colnames(loops))
if (length(missing) > 0) {
  cat("Available columns:", paste(colnames(loops), collapse=", "), "\n")
  stop("Missing columns: ", paste(missing, collapse=", "), call. = FALSE)
}

cat("Loops loaded:", nrow(loops), "\n")

# Clean gene columns - replace NA/"." with empty string
loops <- loops %>%
  mutate(
    Genes_F1 = ifelse(is.na(Genes_F1) | Genes_F1 == ".", "", Genes_F1),
    Genes_F2 = ifelse(is.na(Genes_F2) | Genes_F2 == ".", "", Genes_F2)
  )

##############################################################
# 2. HELPER: PARSE MULTI-GENE FIELDS
# Genes_F1/F2 can contain comma-separated gene lists
##############################################################

# Expand loops so each row = one gene per fragment
# Returns a long dataframe with gene, which_fragment, and original loop index
expand_loop_genes <- function(loops_df) {

  # Fragment 1 genes
  f1 <- loops_df %>%
    mutate(loop_idx = row_number()) %>%
    filter(Genes_F1 != "") %>%
    mutate(genes_split = strsplit(Genes_F1, ",\\s*")) %>%
    unnest(genes_split) %>%
    mutate(gene           = trimws(genes_split),
           gene_fragment  = "F1",
           other_chr      = chr2,
           other_start    = start2,
           other_end      = end2,
           gene_chr       = chr1,
           gene_start     = start1,
           gene_end       = end1) %>%
    select(loop_idx, Name, gene, gene_fragment,
           gene_chr, gene_start, gene_end,
           other_chr, other_start, other_end,
           Genes_F1, Genes_F2)

  # Fragment 2 genes
  f2 <- loops_df %>%
    mutate(loop_idx = row_number()) %>%
    filter(Genes_F2 != "") %>%
    mutate(genes_split = strsplit(Genes_F2, ",\\s*")) %>%
    unnest(genes_split) %>%
    mutate(gene           = trimws(genes_split),
           gene_fragment  = "F2",
           other_chr      = chr1,
           other_start    = start1,
           other_end      = end1,
           gene_chr       = chr2,
           gene_start     = start2,
           gene_end       = end2) %>%
    select(loop_idx, Name, gene, gene_fragment,
           gene_chr, gene_start, gene_end,
           other_chr, other_start, other_end,
           Genes_F1, Genes_F2)

  bind_rows(f1, f2) %>%
    filter(gene != "") %>%
    distinct()
}

cat("Expanding gene-loop associations...\n")
loops_expanded <- expand_loop_genes(loops)
cat("Gene-loop pairs:", nrow(loops_expanded), "\n")
cat("Unique genes in loops:", length(unique(loops_expanded$gene)), "\n")

##############################################################
# 3. LOAD GENE SETS
##############################################################

fix_gene_col <- function(df) {
  if (!"gene" %in% colnames(df) && "X" %in% colnames(df)) {
    df <- df %>% rename(gene = X)
  }
  return(df)
}

# LRT clusters
lrt_df <- read.csv(opt$lrt, stringsAsFactors = FALSE) %>%
  fix_gene_col() %>%
  select(gene, cluster) %>%
  mutate(group      = paste0("LRT_cluster_", cluster),
         group_type = "LRT")

# Expression classification
expr_df <- read.csv(opt$expr, stringsAsFactors = FALSE) %>%
  fix_gene_col() %>%
  filter(gene_set != "Unclassified") %>%
  select(gene, gene_set) %>%
  rename(group = gene_set) %>%
  mutate(group_type = "expr_class")

all_gene_sets <- bind_rows(lrt_df %>% select(gene, group, group_type),
                           expr_df %>% select(gene, group, group_type))

cat("Total gene-group assignments:", nrow(all_gene_sets), "\n")
cat("Groups:", paste(unique(all_gene_sets$group), collapse=", "), "\n")

##############################################################
# 4. FILTER LOOPS AND EXTRACT ENHANCER FRAGMENTS PER GROUP
##############################################################

groups <- unique(all_gene_sets$group)

summary_rows <- list()

for (grp in groups) {

  grp_genes <- all_gene_sets$gene[all_gene_sets$group == grp]
  grp_type  <- unique(all_gene_sets$group_type[all_gene_sets$group == grp])[1]

  # Find all loop-gene pairs where gene is in this group
  matched <- loops_expanded %>%
    filter(gene %in% grp_genes)

  if (nrow(matched) == 0) {
    cat("No loops found for:", grp, "\n")
    next
  }

  # Get unique loops (a loop may match multiple genes)
  unique_loop_indices <- unique(matched$loop_idx)
  n_loops <- length(unique_loop_indices)

  cat(grp, "- genes matched:", length(unique(matched$gene)),
      "| loops:", n_loops, "\n")

  # Safe filename
  safe_grp <- gsub("[^A-Za-z0-9_]", "_", grp)
  subdir   <- file.path(opt$outdir, grp_type)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)

  ##############################################################
  # OUTPUT 1: Full loop records for this group
  ##############################################################

  full_loops <- loops[unique_loop_indices, ] %>%
    mutate(matched_genes = sapply(unique_loop_indices, function(idx) {
      paste(unique(matched$gene[matched$loop_idx == idx]), collapse=";")
    }))

  write.csv(full_loops,
            file.path(subdir, paste0("all_loops_", safe_grp, ".csv")),
            row.names = FALSE)

  ##############################################################
  # OUTPUT 2: Enhancer fragments
  # The fragment NOT containing the gene = putative enhancer
  # If gene in F1 -> enhancer = F2 coords
  # If gene in F2 -> enhancer = F1 coords
  # If gene in both fragments of same loop -> skip that loop
  ##############################################################

  # Find loops where gene is only in F1 or only in F2 (not both)
  loop_gene_frags <- matched %>%
    group_by(loop_idx) %>%
    summarise(
      frags_with_gene = paste(sort(unique(gene_fragment)), collapse=","),
      genes_matched   = paste(unique(gene), collapse=";"),
      .groups = "drop"
    )

  # Loops where gene maps to F1 only -> enhancer is F2
  f1_only <- loop_gene_frags %>%
    filter(frags_with_gene == "F1") %>%
    inner_join(matched %>% filter(gene_fragment == "F1") %>%
                 select(loop_idx, other_chr, other_start, other_end,
                        gene_chr, gene_start, gene_end, gene) %>%
                 distinct(loop_idx, .keep_all = TRUE),
               by = "loop_idx") %>%
    mutate(enhancer_chr   = other_chr,
           enhancer_start = other_start,
           enhancer_end   = other_end,
           promoter_chr   = gene_chr,
           promoter_start = gene_start,
           promoter_end   = gene_end)

  # Loops where gene maps to F2 only -> enhancer is F1
  f2_only <- loop_gene_frags %>%
    filter(frags_with_gene == "F2") %>%
    inner_join(matched %>% filter(gene_fragment == "F2") %>%
                 select(loop_idx, other_chr, other_start, other_end,
                        gene_chr, gene_start, gene_end, gene) %>%
                 distinct(loop_idx, .keep_all = TRUE),
               by = "loop_idx") %>%
    mutate(enhancer_chr   = other_chr,
           enhancer_start = other_start,
           enhancer_end   = other_end,
           promoter_chr   = gene_chr,
           promoter_start = gene_start,
           promoter_end   = gene_end)

  n_both_skipped <- sum(loop_gene_frags$frags_with_gene == "F1,F2")
  cat("  Loops skipped (gene in both fragments):", n_both_skipped, "\n")

  enhancer_df <- bind_rows(f1_only, f2_only)

  if (nrow(enhancer_df) == 0) {
    cat("  No enhancer fragments extracted for:", grp, "\n")
    next
  }

  # BED format: chr, start, end, name (loop_Name:gene), score, strand
  enhancer_bed <- enhancer_df %>%
    mutate(
      name  = paste0(loops$Name[loop_idx], ":", genes_matched),
      score = 0,
      strand = "."
    ) %>%
    select(enhancer_chr, enhancer_start, enhancer_end, name, score, strand) %>%
    distinct() %>%
    arrange(enhancer_chr, enhancer_start)

  write.table(enhancer_bed,
              file.path(subdir, paste0("enhancer_fragments_", safe_grp, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Also save a detailed version with both fragments and gene info
  enhancer_detail <- enhancer_df %>%
    mutate(loop_name = loops$Name[loop_idx]) %>%
    select(loop_name, genes_matched,
           enhancer_chr, enhancer_start, enhancer_end,
           promoter_chr, promoter_start, promoter_end) %>%
    distinct()

  write.csv(enhancer_detail,
            file.path(subdir, paste0("enhancer_fragments_detail_", safe_grp, ".csv")),
            row.names = FALSE)

  cat("  Enhancer fragments:", nrow(enhancer_bed), "\n")

  summary_rows[[grp]] <- data.frame(
    group           = grp,
    group_type      = grp_type,
    genes_in_group  = length(grp_genes),
    genes_in_loops  = length(unique(matched$gene)),
    total_loops     = n_loops,
    loops_skipped   = n_both_skipped,
    enhancer_frags  = nrow(enhancer_bed)
  )
}

##############################################################
# 5. SUMMARY TABLE
##############################################################

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  write.csv(summary_df,
            file.path(opt$outdir, "loop_filter_summary.csv"),
            row.names = FALSE)
  cat("\n--- Summary ---\n")
  print(summary_df)
}

cat("\n--- Loop filtering complete ---\n")
cat("Output directory:", opt$outdir, "\n")
cat("Per group:\n")
cat("  all_loops_<group>.csv             - full loop records\n")
cat("  enhancer_fragments_<group>.bed    - putative enhancer coordinates\n")
cat("  enhancer_fragments_detail_<group>.csv - enhancer + promoter + gene info\n")
