# TEtranscripts post-processing plots
# Run this script from the practical4 directory, which contains the results/ link

# First, let's load the complete output table for the DE analysis
DEresults <- read.table("results/piwiKD_vs_control_gene_TE_analysis.txt")
# we can look at the top of the table
head(DEresults)
# and the end of the table
tail(DEresults)
# Let's create a column that indicates if the feature is a gene or a TE
DEresults$feature <- ifelse(grepl(pattern = "FBgn", rownames(DEresults)), # test if the name contains "FBgn"
                            "gene", # if yes it is a gene
                            "TE")   # else, it is a TE!
DEresults$feature <- factor(DEresults$feature, levels = c("gene", "TE")) # genes in black, TEs in red
DEresults <- DEresults[order(DEresults$feature), ] # genes first, so the TEs are drawn on top
# DESeq2 gives no adjusted P-value (NA) to the features it filtered out, mostly those with low counts; we leave them out
DEresults <- subset(DEresults, !is.na(padj))
# We make a subset of this table that only includes the TEs (easier to use with R-base plot function)
DETEs <- subset(DEresults, feature == "TE")
# and a subset of the TEs with Padj < 0.05 and a fold change of more than 4 (|log2FC| > 2)
DETEssig <- subset(DETEs, abs(log2FoldChange) > 2 & padj < 0.05)

# Now let's plot the very anticipated volcano plot!
# First we set the window for two plots
par(mfrow = c(1, 2))
# Then we plot a general volcano plot representing all genes and TEs
plot(x = DEresults$log2FoldChange, y = -log10(DEresults$padj),
     pch = 16,                              # dots instead of circles
     col = DEresults$feature,               # color according to gene or TE
     main = "Piwi KD vs Control",           # Title and axis names
     xlab = "log2FC",
     ylab = "-log10(Padj)")
abline(h = -log10(0.05), col = "red", lty = 2) # plot the Padj < 0.05 threshold
abline(v = -2, col = "blue", lty = 2)          # plot the |log2FC| > 2 thresholds
abline(v = 2, col = "blue", lty = 2)
legend("topright", levels(DEresults$feature), col = 1:2, pch = 16) # add a legend

# And now, another plot that focuses on TEs
plot(x = DETEs$log2FoldChange, y = -log10(DETEs$padj), # first plot all TEs
     pch = 16,                                         # dots instead of circles
     col = "grey",                                     # color all in grey as background
     main = "Piwi KD vs Control (TE only)",            # Title and axis names
     xlab = "log2FC",
     ylab = "-log10(Padj)",
     xlim = c(-5, 7))                                  # make space for TE labels
abline(h = -log10(0.05), col = "red", lty = 2)         # plot the Padj < 0.05 threshold
abline(v = -2, col = "blue", lty = 2)                  # plot the |log2FC| > 2 thresholds
abline(v = 2, col = "blue", lty = 2)
points(x = DETEssig$log2FoldChange, y = -log10(DETEssig$padj), # add the significant data points
       pch = 16,                                               # dots instead of circles
       col = "red")                                            # color them in red
text(x = DETEssig$log2FoldChange, y = -log10(DETEssig$padj),   # add labels for significant TEs
     labels = rownames(DETEssig), pos = 4, cex = .7)
legend("topleft", c("other TE families", "Padj < 0.05 and |log2FC| > 2"),
       col = c("grey", "red"), pch = 16, cex = .8) # add a legend
