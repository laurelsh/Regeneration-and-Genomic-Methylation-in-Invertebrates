###############################################################################
# Phylogenetic Comparative Analysis
# Runs PGLS on methylation vs regeneration for one tree.
###############################################################################

library(ape)
library(phytools)
library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)
library(rcompanion)

# Load data and tree ----------------------------------------------------------
# Set this to your local folder containing the data and trees:
base_dir <- "WRITE_YOUR_PATH_HERE"

# Choose which tree to analyze (uncomment one)
tree_file <- file.path(base_dir, "data", "Tree_calibrated.nwk")
# tree_file <- file.path(base_dir, "data", "Tree_grafted.nwk")

data_file <- file.path(base_dir, "data", "Methylation_data.csv")
out_dir   <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Match species between tree and data -----------------------------------------
tree$tip.label <- tolower(trimws(tree$tip.label))
df$Species     <- tolower(trimws(df$Species))
df <- df[df$Regen_Ability %in% c(1,2,3), ]

common <- intersect(tree$tip.label, df$Species)
tree <- drop.tip(tree, setdiff(tree$tip.label, common))
df   <- df[match(tree$tip.label, df$Species), ]
df$Species <- tree$tip.label

# Transform methylation -------------------------------------------------------
lambda_tuk <- transformTukey(df$Methylation, returnLambda = TRUE, plotit = FALSE)
df$Methylation_Tukey <- transformTukey(df$Methylation, returnLambda = FALSE, plotit = FALSE)

# Recode regeneration categories ----------------------------------------------
df$Regen <- factor(df$Regen_Ability, levels = c(1,2,3),
                   labels = c("None","Partial","Whole-body"))

# Fit PGLS model --------------------------------------------------------------
tree$edge.length[tree$edge.length < 1e-6] <- 1e-4
tree <- phytools::force.ultrametric(tree, method = "nnls")

corPag <- corPagel(0.5, phy = tree, form = ~Species, fixed = FALSE)
mod_full <- gls(Methylation_Tukey ~ Regen, data = df,
                correlation = corPag, method = "ML")
mod_null <- gls(Methylation_Tukey ~ 1, data = df,
                correlation = corPag, method = "ML")

anova(mod_null, mod_full)

# Estimated marginal means ----------------------------------------------------
emm_ci <- as.data.frame(confint(emmeans(mod_full, ~Regen, data = df)))
contr  <- as.data.frame(pairs(emmeans(mod_full, ~Regen, data = df)))

# Ensure consistent naming
names(emm_ci)[1] <- "Regen"

# Plot ------------------------------------------------------------------------
p <- ggplot(df, aes(x = Regen, y = Methylation_Tukey)) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 1.5, color = "gray40") +
  geom_errorbar(
    data = emm_ci,
    aes(x = Regen, ymin = lower.CL, ymax = upper.CL),
    width = 0.1, color = "black",
    inherit.aes = FALSE   # <--- key line
  ) +
  geom_point(
    data = emm_ci,
    aes(x = Regen, y = emmean),
    size = 3, color = "black",
    inherit.aes = FALSE   # <--- key line
  ) +
  theme_classic(base_size = 11) +
  labs(
    x = "Regeneration capacity",
    y = "CpG methylation (Tukey-transformed)"
  )

ggsave(file.path(base_dir, "results", "PGLS_grafted_simple.pdf"),
       plot = p, width = 95/25.4, height = 85/25.4, units = "in")

# Save contrast table ---------------------------------------------------------
write.csv(contr, file.path(base_dir, "results", "Contrasts_grafted_simple.csv"),
          row.names = FALSE)

cat("\nDone. Figure and contrasts saved.\n")
