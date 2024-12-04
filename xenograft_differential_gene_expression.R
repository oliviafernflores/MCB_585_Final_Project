# Load libraries
library(limma)
library(tidyverse)

# Read raw data
df <- read.csv("/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_data_normalized.csv")

# Get a matrix of the expression data that only has numerica columns
expr_data <- df %>% column_to_rownames(var = "ID_REF") %>%
  select(-c(Sample_Type, Time_Status)) %>%
  mutate(across(everything(), ~as.numeric(.))) %>%
  t()
# expr_data <- expr_data[, -ncol(expr_data)]

# Calculate differential gene expression
group <- factor(c(rep("long_TTL", 7), rep("short_TTL", 5)))
design <- model.matrix(~group)
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, number=Inf)

# Filter the genes that are up or down regulated and statistically significant
upregulated_genes <- results %>%
  filter(logFC > 0 & P.Value < 0.05)

downregulated_genes <- results %>%
  filter(logFC < 0 & P.Value < 0.05)

# Turn the gene filters into tables that can be saved, formatting is nicer too
upregulated_table <- upregulated_genes %>%
  mutate(Probe = rownames(upregulated_genes)) %>%
  mutate(Probe = sub("^X", "", Probe)) %>% 
  select(Probe, Gene = Gene_ID, logFC)
rownames(upregulated_table) <- seq_along(rownames(upregulated_table))

downregulated_table <- downregulated_genes %>%
  mutate(Probe = rownames(downregulated_genes)) %>%
  mutate(Probe = sub("^X", "", Probe)) %>%
  select(Probe, Gene = Gene_ID, logFC)
rownames(downregulated_table) <- seq_along(rownames(downregulated_table))

# Making the log fold change volcano plot

# Assign a default color for genes
results$color <- "grey"

# Color points green if they meet the up regulation criteria, red for down regulation
results$color[results$logFC > 0 & results$P.Value < 0.05] <- "green"
results$color[results$logFC < 0 & results$P.Value < 0.05] <- "red"

# Make a factor for the color labels that will be used for the legend
results$color_status <- factor(results$color, levels = c("green", "red", "grey"), 
                               labels = c("Upregulated", "Downregulated", "Not Significant"))

# Plot the volcano plot
ggplot(results, aes(x = logFC, y = -log10(P.Value))) + 
  geom_point(aes(color = color_status), 
             position = position_jitter(width = 0.2, height = 0.1), alpha = 0.6, size = 1) + 
  scale_color_manual(values = c("Upregulated" = "green", "Downregulated" = "red", "Not Significant" = "grey")) +
  theme_minimal() + 
  labs(
    title = "Differential Gene (Probe) Expression",
    x = "Log2 Fold Change",
    y = "-log10(P-value)",
    color = "Gene Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  )

# Save the differential gene expression results as csv files.
write_csv(downregulated_table, "/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_downregulated.csv")
write_csv(upregulated_table, "/Users/olivia/Documents/Fall 2024/MCB 585/Final Project/xenograft_upregulated.csv")
