#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript DESEQ2_Analysis_Ld.R <counts.tsv> <colData.tsv> [--cond X --ctrl Y --cond_title A --ctrl_title B]", call. = FALSE)
}

counts_file <- args[1]
colData_file <- args[2]

# ---- parse optional flags ----
get_opt <- function(flag, default=NULL) {
  if (!(flag %in% args)) return(default)
  i <- match(flag, args)
  if (i == length(args)) return(default)
  return(args[i+1])
}

cond <- get_opt("--cond", NULL)
ctrl <- get_opt("--ctrl", NULL)
cond_title <- get_opt("--cond_title", cond)
ctrl_title <- get_opt("--ctrl_title", ctrl)

if (is.null(cond) || is.null(ctrl)) {
  stop("You must provide --cond and --ctrl (matching values in colData$Condition).", call. = FALSE)
}
if (is.null(cond_title)) cond_title <- cond
if (is.null(ctrl_title)) ctrl_title <- ctrl

# ---- install/load packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("dplyr", "ggplot2", "ggrepel")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

bioc_pkgs <- c("DESeq2")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
options(ggrepel.max.overlaps = Inf)

# ---- repo paths for DB ----
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_args[grep("^--file=", cmd_args)]
if (length(file_arg) == 0) {
  # fallback: assume script is run from repo root
  script_dir <- getwd()
} else {
  script_path <- sub("^--file=", "", file_arg[1])
  script_dir <- dirname(normalizePath(script_path))
}
repo_dir <- normalizePath(file.path(script_dir, ".."))

desc_file <- file.path(repo_dir, "DB", "LD_Id_Description.txt")

# ---- read inputs ----
counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# read colData with comment.char default '#' so commented samples are ignored
colData <- read.table(colData_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

stopifnot("SampleID" %in% colnames(colData))
stopifnot("Condition" %in% colnames(colData))

# Filter colData to only samples present in counts
colData <- colData[colData$SampleID %in% colnames(counts), , drop = FALSE]
counts <- counts[, colData$SampleID, drop = FALSE]
rownames(colData) <- colData$SampleID

stopifnot(all(colnames(counts) == rownames(colData)))

colData$Condition <- as.factor(colData$Condition)
if ("Batch" %in% colnames(colData)) colData$Batch <- as.factor(colData$Batch)

# Ensure cond/ctrl exist after filtering
if (!(cond %in% levels(colData$Condition))) stop(paste("Condition not found in colData after filtering:", cond))
if (!(ctrl %in% levels(colData$Condition))) stop(paste("Condition not found in colData after filtering:", ctrl))

# ---- DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Condition)
dds <- DESeq(dds)

rld <- varianceStabilizingTransformation(dds, blind = TRUE)
p <- plotPCA(rld, intgroup = "Condition") +
  geom_text_repel(aes(label = colnames(rld)), vjust = 2, size = 2.0, show.legend = FALSE)

ggsave("PCA_plot.png", plot = p, width = 8, height = 6, units = "in", dpi = 600, bg = "white")
cat("PCA plot saved as: PCA_plot.png\n")

# ---- results ----
outputPrefix <- paste0(cond_title, "_vs_", ctrl_title)
res <- results(dds, contrast = c("Condition", cond, ctrl))

resdata <- merge(as.data.frame(res),
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names", sort = FALSE)
names(resdata)[1] <- "gene"

if ("log2FoldChange" %in% names(resdata)) names(resdata)[names(resdata) == "log2FoldChange"] <- "log2FC"

resdata$FC <- sign(resdata$log2FC) * 2^abs(resdata$log2FC)
resdata$FoldChange <- ifelse(resdata$log2FC >= 0, 2^resdata$log2FC, 1 / (2^abs(resdata$log2FC)))

core_cols <- c("gene","baseMean","log2FC","FoldChange","FC","lfcSE","stat","pvalue","padj")
core_cols <- core_cols[core_cols %in% names(resdata)]
resdata <- resdata[, c(core_cols, setdiff(names(resdata), core_cols))]

# Add descriptions if available
if (file.exists(desc_file)) {
  desc_data <- read.table(desc_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  colnames(desc_data) <- c("gene", "Description")
  resdata <- merge(resdata, desc_data, by = "gene", all.x = TRUE)
} else {
  cat("Warning: description file not found (skipping): ", desc_file, "\n")
}

out_csv <- paste0(outputPrefix, "-results-with-normalized.csv")
write.csv(resdata, file = out_csv, row.names = FALSE)
cat("DESeq2 results written to: ", out_csv, "\n")
