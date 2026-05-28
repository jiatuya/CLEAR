suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(optparse)
  library(Rcpp)
})

option_list <- list(
  make_option(c("-g", "--gene_file"), type = "character"),
  make_option(c("-p", "--preprocess_data"), type = "character", default = NULL),
  make_option(c("-a", "--annotation"), type = "character"),
  make_option(c("-o", "--output_filename"), type = "character"),
  
  make_option(c("-c", "--lower_cutoff"), type = "integer", default = 20),
  make_option(c("-C", "--upper_cutoff"), type = "integer", default = 500),
  
  make_option(c("-i", "--n_iterations"), type = "integer", default = 1000000),
  make_option(c("-b", "--burn_in"), type = "integer", default = 500000),
  make_option(c("-l", "--record_likelihoods"), type = "integer", default = 1000),
  make_option(c("-s", "--seed"), type = "integer", default = 1),
  
  make_option(c("-K", "--n_components"), type = "integer", default = 3),
  make_option(c("-u", "--mixture_update_every"), type = "integer", default = 5),
  make_option(c("-r", "--return_params"), action = "store_true", default = FALSE)
)

args <- parse_args(OptionParser(option_list = option_list))

set.seed(args$seed)

# -----------------------------
# 1. Read annotation
# -----------------------------
annt <- read.csv(args$annotation, header = TRUE, sep = "\t")
colnames(annt) <- c("term", "T")

T_list <- apply(annt["T"], 1, function(x) {
  x <- gsub("\\[|\\]|'", "", x)
  strsplit(x, ", ")[[1]]
})

annt2 <- data.frame(
  term = rep(annt$term, sapply(T_list, length)),
  gene = unlist(T_list)
)

annt2 <- annt2 %>%
  group_by(term) %>%
  filter(n() >= args$lower_cutoff, n() <= args$upper_cutoff) %>%
  ungroup()

GO <- split(annt2$gene, annt2$term)

# -----------------------------
# 2. Read gene-level data
# -----------------------------
if (is.null(args$preprocess_data)) {
  gene_file <- read.table(args$gene_file, header = FALSE, sep = "\t")
  if (ncol(gene_file) > 2) {
    stop("Gene file has more than 2 columns. Please provide -p preprocessing script.")
  }
} else {
  source(args$preprocess_data)
  gene_file <- preprocess_data(args$gene_file)
}

gene_file$V2 <- as.numeric(gene_file$V2)
gene_file$V3 <- as.numeric(gene_file$V3)
gene_file$V4 <- as.numeric(gene_file$V4)

genes <- as.character(gene_file$V1)
p_values <- gene_file$V3

# -----------------------------
# 3. Run NEAT p-value mode
# -----------------------------
source("R/NEATv4.R")

set.seed(args$seed)

res <- NEAT(
  genes = genes,
  stats = p_values,
  GO = GO,
  data_type = "p-values",
  
  n_iterations = args$n_iterations,
  burn_in = args$burn_in,
  record_likelihoods = args$record_likelihoods,
  
  n_components = args$n_components,
  dirichlet_alpha = 1,
  kappa0 = 0.01,
  sigma2_a0 = 2,
  mixture_update_every = args$mixture_update_every,
  initialize_from_all_genes = TRUE,
  
  verbose = FALSE
)

# -----------------------------
# 4. Write pipeline output
# -----------------------------
final_df <- data.frame(
  ID = names(res$on_frequency),
  on_frequency = as.numeric(res$on_frequency)
)

write.table(
  final_df,
  file = args$output_filename,
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

# Optional: save diagnostics for convergence checking
if (isTRUE(args$return_params)) {
  param_file <- gsub("\\.tsv$", "_params.rds", args$output_filename)
  saveRDS(res, file = param_file)
}