#! /usr/bin/env Rscript

###############################################################################
#
# DESCRIPTION
#
# Agentic variant of deseq.R from github.com/akshayparopkari/RNAseq.
# Same DESeq2 logic as the original script. The only substantive change:
# the reference level is a parameter instead of a hardcoded
# `factor(dds$Condition, levels = c("WT", "Mutant"))`, so the agent's
# planning step can choose it deliberately based on the actual Condition
# values in the metadata, rather than the script assuming a specific
# two-level WT/Mutant experiment.
#
# USAGE
#
# Rscript --vanilla deseq_agentic.R gene_counts.txt metadata.xlsx \
#         deseq2_lfc.txt MA_plot.pdf [reference_level]
#
# If reference_level is omitted, it defaults to whichever Condition value
# sorts first alphabetically -- NOT necessarily "WT". This is a real
# behavior change from the original script and should be called out
# explicitly wherever this script is used, since it can flip the sign of
# log2FoldChange relative to the original published results.
#
###############################################################################

library("DESeq2", warn.conflicts = FALSE, quietly = T, verbose = F)
library("readxl", warn.conflicts = FALSE, quietly = T, verbose = F)
library("IHW", warn.conflicts = FALSE, quietly = T, verbose = F)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript deseq.R gene_raw_counts.txt metadata.xlsx deseq2_lfc.txt MA_plot.pdf [reference_level]")
}

coldata <- as.data.frame(read_excel(args[2]))
raw.counts <- as.matrix(read.delim(file = args[1], row.names = 1))
raw.counts <- raw.counts[, coldata$Sample_ID]

# MUST BE TRUE, same sanity checks as the original script
stopifnot(all(coldata$Sample_ID %in% colnames(raw.counts)))
stopifnot(all(coldata$Sample_ID == colnames(raw.counts)))

dds <- DESeqDataSetFromMatrix(countData = raw.counts,
                              colData = coldata,
                              design = ~ Condition)

# --- this block replaces the original hardcoded factor levels ---
condition_levels <- unique(as.character(coldata$Condition))
reference_level <- if (length(args) >= 5) args[5] else sort(condition_levels)[1]

if (!(reference_level %in% condition_levels)) {
  stop(paste0("Reference level '", reference_level,
              "' not found in Condition column. Found: ",
              paste(condition_levels, collapse = ", ")))
}

ordered_levels <- c(reference_level, setdiff(condition_levels, reference_level))
dds$Condition <- factor(dds$Condition, levels = ordered_levels)
message(paste0("Using '", reference_level, "' as reference level. Full order: ",
               paste(ordered_levels, collapse = " -> ")))
# --- end replacement ---

print(resultsNames(dds))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
resIHW <- results(dds, filterFun = ihw, alpha = 0.05)
print(summary(resIHW))
print(sum(resIHW$padj < 0.05, na.rm = TRUE))
print(metadata(resIHW)$ihwResult)

resLFC <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")
write.table(resLFC[order(resLFC$padj), ],
            file = args[3],
            quote = FALSE,
            sep = "\t")
pdf(file = args[4])
plotMA(resLFC, ylim = c(-2, 2))
abline(h = c(-1, 1), col = "dodgerblue", lwd = 2)
dev.off()
