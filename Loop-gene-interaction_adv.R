##############################################################
# Loop filtering and enhancer fragment extraction
#
# Input:
#   --loops     : CSV or XLSX with loop interactions
#                 cols: chr1,start1,end1,chr2,start2,end2,Name,Genes_F1,Genes_F2
#   --bed       : BED file with gene promoter regions (for unlooped output)
#   --lrt       : LRT_significant_genes_clustered.csv
#   --expr      : geneset_full_classification.csv
#   --outdir    : output directory
#
# Output per cluster/geneset:
#   1) all_loops_<group>.csv              - all loops containing a gene from that group
#   2) enhancer_fragments_<group>.bed       - putative enhancer fragments (opposite to gene)
#   3) unlooped_promoters_<group>.bed       - promoters with no loop evidence (requires --bed)
#
# Usage:
#   Rscript Loop-gene_interaction.R \
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
  library(readxl)
})

option_list <- list(
  make_option("--loops",
              type    = "character",
              default = NULL,
              help    = "Loop interaction file (.csv or .xlsx)"),
  make_option("--lrt",
              type    = "character",
              default = "LRT_significant_genes_clustered.csv"),
  make_option("--expr",
              type    = "character",
              default = "geneset_full_classification.csv"),
  make_option("--outdir",
              type    = "character",
              default = "loop_results"),
  make_option("--bed",
              type    = "character",
              default = NULL,
              help    = "BED file with gene promoter regions (chr,start,end,gene,score,strand). Used for unlooped promoter output.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$loops)) stop("--loops is required", call. = FALSE)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

##############################################################
# 1. LOAD LOOPS
##############################################################

cat("Loading loop file:", opt$loops, "\n")

# Determine file type by extension
loops_ext <- tolower(tools::file_ext(opt$loops))

if (loops_ext %in% c("xlsx", "xls")) {
  # Read from Excel
  loops <- read_xlsx(opt$loops)
  loops <- as.data.frame(loops, stringsAsFactors = FALSE)
} else {
  # Default: robust CSV loading to handle environments where read.csv
  # treats the entire header as a single column name
  loops <- tryCatch(
    {
      read.csv(opt$loops, stringsAsFactors = FALSE, check.names = FALSE)
    },
    error = function(e) {
      cat("Standard read.csv failed, trying manual CSV parsing...\n")
      NULL
    }
  )

  # Fallback: manually parse comma-separated file if we ended up
  # with a single column whose name contains commas
  if (is.null(loops) || (ncol(loops) == 1 && grepl(",", colnames(loops)[1]))) {
    cat("Detected single-column table with comma-separated header; reparsing manually.\n")
    raw_lines <- readLines(opt$loops)
    if (length(raw_lines) < 2) {
      stop("Loop file appears to be empty or missing data rows.", call. = FALSE)
    }

    header_fields <- trimws(strsplit(raw_lines[1], ",", fixed = TRUE)[[1]])
    body_lines    <- raw_lines[-1]
    split_body    <- strsplit(body_lines, ",", fixed = TRUE)
    max_len       <- max(lengths(split_body))

    split_body <- lapply(split_body, function(x) {
      length(x) <- max_len
      x
    })

    mat <- do.call(rbind, split_body)
    loops <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(loops) <- header_fields
  }
}

# Standardise column names - handle variations
colnames(loops) <- trimws(colnames(loops))

# Drop columns with empty or NA names (e.g. trailing commas in CSV header)
valid_cols <- !is.na(colnames(loops)) & colnames(loops) != ""
if (!all(valid_cols)) {
  cat("Dropping unnamed columns at positions:",
      paste(which(!valid_cols), collapse = ", "), "\n")
  loops <- loops[, valid_cols, drop = FALSE]
}

# Expected: chr1,start1,end1,chr2,start2,end2,Name,Genes_F1,Genes_F2
required <- c("chr1","start1","end1","chr2","start2","end2","Name","Genes_F1","Genes_F2")
missing  <- setdiff(required, colnames(loops))
if (length(missing) > 0) {
  cat("Available columns:", paste(colnames(loops), collapse=", "), "\n")
  stop("Missing columns: ", paste(missing, collapse=", "), call. = FALSE)
}

cat("Loops loaded:", nrow(loops), "\n")

# Ensure gene columns are character and clean values
if ("Genes_F1" %in% colnames(loops)) {
  loops$Genes_F1 <- as.character(loops$Genes_F1)
}
if ("Genes_F2" %in% colnames(loops)) {
  loops$Genes_F2 <- as.character(loops$Genes_F2)
}

# Clean gene columns - replace NA, ".", "NA"/"na", or empty/whitespace with ""
loops <- loops %>%
  mutate(
    Genes_F1 = ifelse(
      is.na(Genes_F1) |
        trimws(Genes_F1) == "" |
        trimws(Genes_F1) %in% c(".", "NA", "na"),
      "",
      trimws(Genes_F1)
    ),
    Genes_F2 = ifelse(
      is.na(Genes_F2) |
        trimws(Genes_F2) == "" |
        trimws(Genes_F2) %in% c(".", "NA", "na"),
      "",
      trimws(Genes_F2)
    )
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

##############################################################
# 4. LOAD BED FILE (GENE PROMOTERS) FOR UNLOOPED REGIONS
##############################################################

bed_df <- NULL
if (!is.null(opt$bed)) {
  cat("Loading BED file:", opt$bed, "\n")
  bed_df <- read.table(opt$bed, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE)
  colnames(bed_df)[1:4] <- c("chr", "start", "end", "gene")
  if (ncol(bed_df) >= 6) colnames(bed_df)[5:6] <- c("score", "strand")
  bed_df$gene <- trimws(bed_df$gene)
  cat("BED regions loaded:", nrow(bed_df), "\n")
  cat("Unique genes in BED:", length(unique(bed_df$gene)), "\n")
} else {
  cat("No --bed file provided - unlooped promoter output will be skipped\n")
}

# All genes with any loop across the whole dataset
genes_in_any_loop <- unique(loops_expanded$gene)

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
  cat("  Loops with genes in both fragments:", n_both_skipped, "\n")

  both_loop_idxs <- loop_gene_frags %>%
    filter(frags_with_gene == "F1,F2") %>%
    pull(loop_idx)

  n_promoter_promoter <- 0
  n_intragenic        <- 0

  if (length(both_loop_idxs) > 0) {

    # Get gene annotations per fragment per loop
    f1_genes_both <- matched %>%
      filter(loop_idx %in% both_loop_idxs, gene_fragment == "F1") %>%
      group_by(loop_idx) %>%
      summarise(genes_F1 = paste(sort(unique(gene)), collapse = ";"),
                .groups = "drop")

    f2_genes_both <- matched %>%
      filter(loop_idx %in% both_loop_idxs, gene_fragment == "F2") %>%
      group_by(loop_idx) %>%
      summarise(genes_F2 = paste(sort(unique(gene)), collapse = ";"),
                .groups = "drop")

    both_loops_coords <- loops[both_loop_idxs, ] %>%
      mutate(loop_idx = both_loop_idxs) %>%
      left_join(f1_genes_both, by = "loop_idx") %>%
      left_join(f2_genes_both, by = "loop_idx") %>%
      select(loop_idx, Name,
             chr1, start1, end1, genes_F1,
             chr2, start2, end2, genes_F2)

    # Classify: intragenic = same gene(s) on both sides
    #           promoter-promoter = different genes on each side
    both_loops_coords <- both_loops_coords %>%
      mutate(
        is_intragenic = mapply(function(g1, g2) {
          s1 <- unlist(strsplit(g1, ";"))
          s2 <- unlist(strsplit(g2, ";"))
          length(intersect(s1, s2)) > 0 && setequal(s1, s2)
        }, genes_F1, genes_F2)
      )

    intragenic_coords  <- both_loops_coords %>% filter( is_intragenic)
    prom_prom_coords   <- both_loops_coords %>% filter(!is_intragenic)

    cat("    Promoter-promoter loops:", nrow(prom_prom_coords), "\n")
    cat("    Intragenic loops:       ", nrow(intragenic_coords), "\n")

    ##################################################################
    # PROMOTER-PROMOTER OUTPUT
    # Same as enhancer output format:
    # Gene of interest fragment = "anchor", opposite fragment = output BED
    # For each group gene in F1 -> output F2 HiChIP fragment coords
    # For each group gene in F2 -> output F1 HiChIP fragment coords
    ##################################################################

    if (nrow(prom_prom_coords) > 0) {

      pp_f1 <- matched %>%
        filter(loop_idx %in% prom_prom_coords$loop_idx,
               gene_fragment == "F1") %>%
        distinct(loop_idx, gene, .keep_all = TRUE) %>%
        left_join(prom_prom_coords %>%
                    select(loop_idx, genes_F1, genes_F2),
                  by = "loop_idx") %>%
        mutate(paired_promoter_genes = genes_F2)

      pp_f2 <- matched %>%
        filter(loop_idx %in% prom_prom_coords$loop_idx,
               gene_fragment == "F2") %>%
        distinct(loop_idx, gene, .keep_all = TRUE) %>%
        left_join(prom_prom_coords %>%
                    select(loop_idx, genes_F1, genes_F2),
                  by = "loop_idx") %>%
        mutate(paired_promoter_genes = genes_F1)

      pp_all <- bind_rows(pp_f1, pp_f2)

      # BED: opposite HiChIP fragment coordinates (same format as enhancer BED)
      pp_bed <- pp_all %>%
        mutate(
          name   = paste0(loops$Name[loop_idx], ":", gene,
                          "<->", paired_promoter_genes),
          score  = 0,
          strand = ".",
          other_start = format(as.integer(other_start),
                               scientific = FALSE, trim = TRUE),
          other_end   = format(as.integer(other_end),
                               scientific = FALSE, trim = TRUE)
        ) %>%
        select(other_chr, other_start, other_end, name, score, strand) %>%
        distinct() %>%
        arrange(other_chr, other_start)

      write.table(pp_bed,
                  file.path(subdir,
                            paste0("promoter_promoter_fragments_", safe_grp, ".bed")),
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

      # Detail CSV
      pp_detail <- pp_all %>%
        mutate(loop_name = loops$Name[loop_idx]) %>%
        select(loop_name, gene, paired_promoter_genes,
               other_chr, other_start, other_end,
               gene_chr, gene_start, gene_end) %>%
        rename(paired_fragment_chr   = other_chr,
               paired_fragment_start = other_start,
               paired_fragment_end   = other_end,
               gene_fragment_chr     = gene_chr,
               gene_fragment_start   = gene_start,
               gene_fragment_end     = gene_end) %>%
        distinct()

      write.csv(pp_detail,
                file.path(subdir,
                          paste0("promoter_promoter_detail_", safe_grp, ".csv")),
                row.names = FALSE)

      n_promoter_promoter <- nrow(pp_bed)
      cat("  Promoter-promoter fragments (BED):", n_promoter_promoter, "\n")
    }

    ##################################################################
    # INTRAGENIC OUTPUT
    # Same gene on both sides -> output BED coordinates from input
    # BED file for that gene (not the HiChIP fragment)
    ##################################################################

    if (nrow(intragenic_coords) > 0 && !is.null(bed_df)) {

      intragenic_genes <- unique(unlist(
        strsplit(intragenic_coords$genes_F1, ";")
      ))

      intragenic_bed <- bed_df %>%
        filter(gene %in% intragenic_genes) %>%
        distinct() %>%
        arrange(chr, start)

      write.table(intragenic_bed,
                  file.path(subdir,
                            paste0("intragenic_loops_", safe_grp, ".bed")),
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

      n_intragenic <- nrow(intragenic_bed)
      cat("  Intragenic loop gene regions (BED):", n_intragenic, "\n")

    } else if (nrow(intragenic_coords) > 0) {
      # No BED file - save gene list only
      intragenic_genes_df <- data.frame(
        gene = unique(unlist(strsplit(intragenic_coords$genes_F1, ";")))
      )
      write.csv(intragenic_genes_df,
                file.path(subdir,
                          paste0("intragenic_loops_", safe_grp, ".csv")),
                row.names = FALSE)
      cat("  Intragenic genes saved (provide --bed for coordinates)\n")
      n_intragenic <- nrow(intragenic_genes_df)
    }
  }

  enhancer_df <- bind_rows(f1_only, f2_only)

  if (nrow(enhancer_df) == 0) {
    cat("  No enhancer fragments extracted for:", grp, "\n")
    next
  }

  # BED format: chr, start, end, name (loop_Name:gene), score, strand
  # Use integer coordinates without scientific notation (required by deepTools etc.)
  enhancer_bed <- enhancer_df %>%
    mutate(
      name  = paste0(loops$Name[loop_idx], ":", genes_matched),
      score = 0,
      strand = ".",
      enhancer_start = format(as.integer(enhancer_start), scientific = FALSE, trim = TRUE),
      enhancer_end   = format(as.integer(enhancer_end),   scientific = FALSE, trim = TRUE)
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


  ##############################################################
  # OUTPUT 3: UNLOOPED PROMOTERS
  # Genes in this group that are in the BED file but have
  # NO loop link anywhere in the full loop dataset
  ##############################################################

  n_unlooped <- 0
  if (!is.null(bed_df)) {
    genes_without_any_loop <- setdiff(grp_genes, genes_in_any_loop)
    unlooped_bed <- bed_df %>%
      filter(gene %in% genes_without_any_loop) %>%
      distinct() %>%
      arrange(chr, start)

    n_unlooped <- nrow(unlooped_bed)

    if (n_unlooped > 0) {
      write.table(unlooped_bed,
                  file.path(subdir, paste0("unlooped_promoters_", safe_grp, ".bed")),
                  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      cat("  Unlooped promoters:", n_unlooped,
          "(genes with no loop in dataset)\n")
    } else {
      cat("  Unlooped promoters: 0\n")
    }
  }

  summary_rows[[grp]] <- data.frame(
    group              = grp,
    group_type         = grp_type,
    genes_in_group     = length(grp_genes),
    genes_in_loops     = length(unique(matched$gene)),
    genes_unlooped     = n_unlooped,
    total_loops        = n_loops,
    promoter_promoter  = n_promoter_promoter,
    intragenic_loops   = n_intragenic,
    enhancer_frags     = nrow(enhancer_bed)
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
cat("  all_loops_<group>.csv                  - full loop records\n")
cat("  enhancer_fragments_<group>.bed         - putative enhancer coordinates\n")
cat("  enhancer_fragments_detail_<group>.csv  - enhancer + promoter + gene info\n")
cat("  unlooped_promoters_<group>.bed         - promoters with no loop evidence\n")
cat("  promoter_promoter_fragments_<group>.bed - opposite HiChIP fragment for promoter-promoter loops\n")
cat("  promoter_promoter_detail_<group>.csv   - paired promoter gene info\n")
cat("  intragenic_loops_<group>.bed           - BED coords of genes with intragenic loops\n")