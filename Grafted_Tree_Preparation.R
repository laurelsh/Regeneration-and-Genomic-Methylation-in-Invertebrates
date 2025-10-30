# Grafting Species to Make Full Calibrated Tree

###############################################################################
## STEP 1: Load Libraries & Prepare Calibrated Tree
###############################################################################

library(ape)
library(phytools)
library(diversitree)
library(treeio)
library(ggtree)
library(ggplot2)
library(rgbif)
library(dplyr)
library(stringr)
library(pbapply)

# 1.1 Read the existing time-calibrated tree
existing_tree <- read.tree("Calibrated_tree.nwk")

# 1.2 Replace extremely small or zero branch lengths
small_branches <- which(existing_tree$edge.length < 1e-6)
if (length(small_branches) > 0) {
  existing_tree$edge.length[small_branches] <- 1e-4
}
zero_branches <- which(existing_tree$edge.length == 0)
if (length(zero_branches) > 0) {
  existing_tree$edge.length[zero_branches] <- 1e-4
}

# 1.3 Force the tree to be ultrametric
tree_ultrametric <- phytools::force.ultrametric(existing_tree, method = "nnls")

# Rename for clarity: we keep working with 'tree_ultrametric'
existing_tree <- tree_ultrametric

# 1.4 Function to safely remove/initialize node labels
if (is.null(existing_tree$node.label)) {
  existing_tree$node.label <- rep("", existing_tree$Nnode)
} else {
  # Optional: clear them if you want a fresh start
  existing_tree$node.label <- rep("", existing_tree$Nnode)
}

###############################################################################
## STEP 2: Retrieve Taxonomy for the Existing Tips
###############################################################################

# 2.1 Grab tip labels
existing_tips <- existing_tree$tip.label

# 2.2 Create a function to retrieve taxonomy from GBIF
retrieve_taxonomy_rgbif <- function(species_name) {
  species_query <- str_replace_all(species_name, "_", " ")
  out <- tryCatch(
    {
      name_backbone(name = species_query, kingdom = "Animalia")
    },
    error = function(e) {
      message("Error for ", species_name, ": ", e$message)
      return(NULL)
    }
  )
  
  # If the query returned NULL or 0 rows, return all NA
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame(
      Species = species_name,
      Genus   = NA_character_,
      Family  = NA_character_,
      Order   = NA_character_,
      Class   = NA_character_,
      Phylum  = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  # A small helper to safely extract columns if they exist
  safe_extract <- function(df, col) {
    if (!col %in% names(df)) return(NA_character_)
    df[1, col, drop=TRUE]  # always row 1
  }
  
  data.frame(
    Species = species_name,
    Genus   = safe_extract(out, "genus"),
    Family  = safe_extract(out, "family"),
    Order   = safe_extract(out, "order"),
    Class   = safe_extract(out, "class"),
    Phylum  = safe_extract(out, "phylum"),
    stringsAsFactors = FALSE
  )
}

# 2.3 Retrieve taxonomy for the existing tips
taxonomy_list_existing <- pblapply(existing_tips, retrieve_taxonomy_rgbif)
taxonomy_existing_df   <- do.call(rbind, taxonomy_list_existing)

###############################################################################
## STEP 3: Label Internal Nodes at Multiple Ranks (Family, Order, etc.)
###############################################################################
#
# We'll define a function that, given a rank ("Family", "Order", "Class", "Phylum"),
# labels all nodes in 'existing_tree' for which at least 2 tips share that rank.

label_nodes_by_rank <- function(tree, taxonomy_df, rank_name) {
  # 1) Keep only tips that have a known value for rank_name
  sub_df <- taxonomy_df %>% filter(!is.na(.data[[rank_name]]))
  
  # 2) Group by that rank
  rank_to_tips <- sub_df %>%
    group_by(.data[[rank_name]]) %>%
    summarise(Tips = list(Species), .groups = "drop")
  
  # 3) For each group, find the MRCA and label it with the rank value
  num_tips <- length(tree$tip.label)
  
  # Initialize an empty vector to store node labels for this rank
  # We'll store them *in parallel* to tree$node.label, so we keep them separate
  node_labels_this_rank <- rep("", tree$Nnode)
  
  for (i in seq_len(nrow(rank_to_tips))) {
    group_label <- rank_to_tips[[rank_name]][i]  # e.g., the Family name
    these_tips  <- unlist(rank_to_tips$Tips[i])
    
    if (length(these_tips) < 2) {
      next  # can't label a single tip's node
    }
    
    # Find the MRCA
    mrca_node <- findMRCA(tree, tips = these_tips, type = "node")
    if (!is.null(mrca_node) && !is.na(mrca_node)) {
      # Overwrite or fill in the label
      node_labels_this_rank[mrca_node - num_tips] <- group_label
    }
  }
  return(node_labels_this_rank)
}

# 3.1 Label for each rank
genus_labels <- label_nodes_by_rank(existing_tree, taxonomy_existing_df, "Genus")
family_labels <- label_nodes_by_rank(existing_tree, taxonomy_existing_df, "Family")
order_labels  <- label_nodes_by_rank(existing_tree, taxonomy_existing_df, "Order")
class_labels  <- label_nodes_by_rank(existing_tree, taxonomy_existing_df, "Class")
phylum_labels <- label_nodes_by_rank(existing_tree, taxonomy_existing_df, "Phylum")

# Now you have 4 vectors, each length = tree$Nnode, that contain 
# Family, Order, Class, Phylum labels (or ""), e.g.:
# family_labels[i] might be "Asteriidae" if node i is that family's MRCA.

# If you want to store them in the tree object, you could do so in a list:
node_label_data <- data.frame(
  Genus  = genus_labels,
  Family = family_labels,
  Order  = order_labels,
  Class  = class_labels,
  Phylum = phylum_labels
)
# This 'node_label_data' has as many rows as internal nodes, in the same order.

###############################################################################
## STEP 4: Estimate Speciation Rates with Diversitree for All Ranks
###############################################################################

# 4.1 A helper function to fit a Yule (pure-birth) model to a node
#     in the tree, identified by internal-node index (1..Nnode).
library(diversitree)

estimate_speciation_rate <- function(tree, node_idx) {
  # Convert node_idx (1..Nnode) to "ape" node numbering
  node_num <- node_idx + length(tree$tip.label)
  
  # Extract the clade
  sub_clade <- extract.clade(tree, node_num)
  
  # If fewer than 2 tips, we cannot estimate a Yule rate
  if (length(sub_clade$tip.label) < 2) {
    return(NA_real_)
  }
  
  # Fit Yule model
  lik <- make.yule(sub_clade)
  fit <- tryCatch(
    find.mle(lik, x.init = 0.1, lower = 0),
    error = function(e) {
      message("Error fitting Yule at node ", node_num, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(fit)) return(NA_real_)
  
  # Return the speciation rate (lambda)
  fit$par["lambda"]
}

# 4.3 Create a data frame that will store all ranks’ speciation rates
all_clade_rates <- data.frame(
  Rank      = character(),
  NodeIndex = integer(),    # 1..Nnode
  Label     = character(),  # e.g. "Asteriidae", "Drosophila", etc.
  SpecRate  = numeric(),
  stringsAsFactors = FALSE
)

# 4.4 We define which ranks we want to handle
all_ranks <- c("Genus", "Family", "Order", "Class", "Phylum")

# 4.5 Loop over each rank and each internal node
num_nodes <- nrow(node_label_data)  # should be tree$Nnode
for (rank_name in all_ranks) {
  
  cat("\nEstimating Yule rates for rank:", rank_name, "\n")
  
  # The vector of labels for this rank
  label_vector <- node_label_data[[rank_name]]
  
  # For each node (1..Nnode) in the tree:
  for (i in seq_len(num_nodes)) {
    lab <- label_vector[i]
    if (lab == "" || is.na(lab)) {
      # No label => skip
      next
    }
    
    # Estimate the speciation rate
    rate_val <- estimate_speciation_rate(existing_tree, i)
    
    # Store the result
    row_to_add <- data.frame(
      Rank      = rank_name,
      NodeIndex = i,
      Label     = lab,
      SpecRate  = rate_val,
      stringsAsFactors = FALSE
    )
    all_clade_rates <- rbind(all_clade_rates, row_to_add)
  }
}

cat("Finished estimating speciation rates across all ranks.\n")
cat("Number of labeled clades (with possible rates):", nrow(all_clade_rates), "\n")

# 4.6 (Optional) Handle missing rates
# You can fill them with the mean across all ranks, or the mean by rank:
grand_mean <- mean(all_clade_rates$SpecRate, na.rm = TRUE)
all_clade_rates$SpecRate[is.na(all_clade_rates$SpecRate)] <- grand_mean

# Now you have a single data frame 'all_clade_rates' with columns:
#   Rank, NodeIndex, Label, SpecRate

###############################################################################
## STEP 5: Identify & Retrieve Taxonomy for New Species
###############################################################################

full_topology_tree <- read.tree("Full_topology.nwk")
full_tips <- full_topology_tree$tip.label

# 5.1 Figure out which species are new
new_species <- setdiff(full_tips, existing_tree$tip.label)
cat("Number of new species to add:", length(new_species), "\n")
print(new_species)

# 5.2 Retrieve taxonomy for those new species
taxonomy_list_new <- pblapply(new_species, retrieve_taxonomy_rgbif)
new_tax_df <- do.call(rbind, taxonomy_list_new)
head(new_tax_df)

###############################################################################
## STEP 6: Fallback Approach to Graft Each New Species
###############################################################################
#
# We'll define a function that attempts to find a node label in the order:
# 1) Genus, 2) Family, 3) Order, 4) Class, 5) Phylum
# If *none* is found, we skip or attach at root. 
# Once a node label is found, we can look up the speciation rate from that rank’s 
# rates data frame, or default to a global mean.

# Let's define a helper that tries to find a node with a given rank_label in 
# the relevant vector (family_labels, order_labels, etc.). 
# Then we can combine them in a single fallback function.

find_node_for_label <- function(rank_vector, rank_label) {
  # rank_vector is something like family_labels
  idx <- which(rank_vector == rank_label)
  if (length(idx) > 0) {
    return(idx[1])  # just return the first match
  } 
  return(NA_integer_)
}

# We'll define a function to get the speciation rate for that node or a default
# from the relevant rates data frame. Let's do it for families as an example.
get_family_rate <- function(fam, family_rates) {
  row_idx <- which(family_rates$Family == fam)
  if (length(row_idx) > 0) {
    return(family_rates$SpecRate[row_idx[1]])
  }
  return(NA_real_)
}

# Now the main fallback function:
graft_species_fallback <- function(tree, species_name, tax_row, 
                                   node_label_data,
                                   family_rates, default_rate = 0.01) {
  # tax_row: a single row with (Genus, Family, Order, Class, Phylum)
  # node_label_data: data frame with columns: Family, Order, Class, Phylum
  # family_rates: data frame with columns: (NodeIndex, Family, SpecRate)
  
  # Fallback order:
  fallback_ranks <- c("Genus", "Family", "Order", "Class", "Phylum")
  
  # Try each rank in turn
  rank_found <- NA
  node_idx   <- NA
  
  for (rank in fallback_ranks) {
    rank_val <- tax_row[[rank]]
    if (is.na(rank_val)) next  # skip if unknown
    
    # Figure out which node_label_vector to use
    if (rank == "Genus") {
      next 
    } else if (rank == "Family") {
      node_idx <- find_node_for_label(node_label_data$Family, rank_val)
    } else if (rank == "Order") {
      node_idx <- find_node_for_label(node_label_data$Order, rank_val)
    } else if (rank == "Class") {
      node_idx <- find_node_for_label(node_label_data$Class, rank_val)
    } else if (rank == "Phylum") {
      node_idx <- find_node_for_label(node_label_data$Phylum, rank_val)
    }
    
    if (!is.na(node_idx) && node_idx > 0) {
      rank_found <- rank
      break
    }
  }
  
  if (is.na(node_idx) || node_idx < 1) {
    # No rank was found, so attach at root or skip
    # Let's attach at root for demonstration
    root_node <- length(tree$tip.label) + 1
    # We'll just use a default small branch length
    updated_tree <- bind.tip(tree,
                             tip.label = species_name,
                             where = root_node,
                             position = 0.001)
    cat("Warning: No rank found for", species_name, 
        " attaching at root.\n")
    return(updated_tree)
  }
  
  # If we *did* find node_idx, let's find a speciation rate
  node_rate <- default_rate
  
  # If rank_found == "Family", we can look up in family_rates
  if (rank_found == "Family") {
    fam_rate <- get_family_rate(tax_row$Family, family_rates)
    if (!is.na(fam_rate)) {
      node_rate <- fam_rate
    }
  }
  
  # If we wanted to store or compute separate rates for Order, Class, etc.
  # you'd define similar data frames for them and do a similar lookup.
  
  # Convert node_idx to actual node number
  node_number <- node_idx + length(tree$tip.label)
  
  # Branch length ~ 1 / spec_rate is one approach
  # Or a small fraction if you don't want huge branches
  new_branch_length <- 1 / node_rate
  
  updated_tree <- bind.tip(tree,
                           tip.label = species_name,
                           where = node_number,
                           position = new_branch_length)
  
  cat("Grafted species", species_name, "onto", rank_found, "=", tax_row[[rank_found]],
      "using rate =", node_rate, "\n")
  return(updated_tree)
}

# This function iterates over all existing tips in the tree 
# (for which you already have taxonomy in taxonomy_existing_df)
# and finds the tip whose taxonomy matches the new species on 
# the largest number of ranks.

find_closest_sibling_in_tree <- function(tax_row, taxonomy_existing_df,
                                         rank_order = c("Genus","Family","Order","Class","Phylum")) {
  # tax_row is something like data.frame(Genus="xxx", Family="yyy", ...)
  # taxonomy_existing_df is the dataframe of your existing tips’ taxonomy 
  #   with columns c("Species", "Genus","Family","Order","Class","Phylum")
  # rank_order = the ranks you want to compare in priority order

  # We'll store (tip_name, match_count)
  best_tip       <- NA_character_
  best_num_match <- -1
  
  for (i in seq_len(nrow(taxonomy_existing_df))) {
    candidate_tip   <- taxonomy_existing_df$Species[i]
    candidate_taxon <- taxonomy_existing_df[i, rank_order, drop=FALSE]
    
    # Count how many ranks match
    num_match <- sum(
      !is.na(tax_row[rank_order]) & 
      !is.na(candidate_taxon[rank_order]) & 
      (tax_row[rank_order] == candidate_taxon[rank_order])
    )
    
    # Keep track of the best so far
    if (num_match > best_num_match) {
      best_num_match <- num_match
      best_tip       <- candidate_tip
    }
  }
  
  return(best_tip)  
}

###############################################################################
## STEP 7: Graft the New Species (Using all_clade_rates for Speciation Rates)
###############################################################################
library(phytools)  # for getParent()

graft_species_fallback <- function(tree, species_name, tax_row, 
                                   node_label_data, all_clade_rates, default_rate,
                                   taxonomy_existing_df) {
  # We'll check these ranks in order (you can change the order if you prefer)
  fallback_ranks <- c("Genus","Family","Order","Class","Phylum")
  
  node_idx   <- NA_integer_
  rank_found <- NA_character_
  
  #--------------------------------------------------------------------
  # 1) Try each rank -> find node_idx in node_label_data
  #--------------------------------------------------------------------
  for (rank in fallback_ranks) {
    rank_val <- tax_row[[rank]]
    if (is.na(rank_val) || rank_val == "") {
      next
    }
    
    # Identify which label vector to use
    label_vector <- switch(
      rank,
      "Genus"  = node_label_data$Genus,
      "Family" = node_label_data$Family,
      "Order"  = node_label_data$Order,
      "Class"  = node_label_data$Class,
      "Phylum" = node_label_data$Phylum,
      character()  # default empty if none matched
    )
    
    if (length(label_vector) == 0) {
      next
    }
    
    idx <- which(label_vector == rank_val)
    if (length(idx) > 0) {
      node_idx   <- idx[1]  # pick the first matching node
      rank_found <- rank
      break
    }
  }
  
  #--------------------------------------------------------------------
  # 2) If no node found, attach to the best "closest sibling" tip 
  #    instead of the root.
  #--------------------------------------------------------------------
  if (is.na(node_idx) || node_idx < 1) {
    # (A) Find the best sibling tip among existing tips
    best_tip_label <- find_closest_sibling_in_tree(
      tax_row              = tax_row,
      taxonomy_existing_df = taxonomy_existing_df,
      rank_order           = fallback_ranks
    )
    
    if (is.na(best_tip_label)) {
      # if for some reason that fails, revert to root
      cat("No node found and couldn't identify a sibling for", species_name,
          "; attaching at root.\n")
      root_node <- length(tree$tip.label) + 1
      updated_tree <- bind.tip(
        tree, 
        tip.label = species_name, 
        where     = root_node, 
        position  = 0
      )
      return(updated_tree)
    } else {
      # (B) Attach as a sister to that tip’s parent node
      tip_index   <- which(tree$tip.label == best_tip_label)
      parent_node <- getParent(tree, tip_index)

      # We'll choose a small branch length or use your default_rate logic
      node_rate          = default_rate
      new_branch_length  = 1 / node_rate

      # Just to be safe, clamp if longer than the parent's edge
      parent_edge_length <- tree$edge.length[ which(tree$edge[,2] == parent_node) ]
      if (length(parent_edge_length) == 0) {
        # Sometimes the parent is root & doesn't have an edge above
        new_branch_length <- 0
      } else if (new_branch_length > parent_edge_length) {
        cat(sprintf(
          "Clamping branch for %s from %.2f to 0 (polytomy)\n",
          species_name, new_branch_length
        ))
        new_branch_length <- 0
      }

      cat(sprintf(
        "No clade match for %s; attaching as sibling to tip '%s'.\n",
        species_name, best_tip_label
      ))

      updated_tree <- bind.tip(
        tree,
        tip.label = species_name,
        where     = parent_node,
        position  = new_branch_length
      )
      return(updated_tree)
    }
  }
  
  #--------------------------------------------------------------------
  # 3) Otherwise, we DID find a node_idx in some rank => proceed
  #--------------------------------------------------------------------
  
  # Convert node_idx -> 'ape' node_number
  node_number <- node_idx + length(tree$tip.label)
  
  # Determine speciation rate (lambda) from all_clade_rates
  node_rate <- default_rate   # fallback
  
  if (!is.na(rank_found)) {
    # e.g. subset by Rank=="Family" & Label=="Felidae"
    sub_df <- subset(all_clade_rates, 
                     Rank == rank_found & Label == tax_row[[rank_found]])
    if (nrow(sub_df) > 0) {
      rate_val <- sub_df$SpecRate[1]  # or mean(sub_df$SpecRate)
      if (!is.na(rate_val)) {
        node_rate <- rate_val
      }
    }
  }
  
  new_branch_length <- 1 / node_rate
  
  # If node_number is the root, clamp to polytomy
  root_node <- length(tree$tip.label) + 1
  if (node_number == root_node) {
    new_branch_length <- 0
  } else {
    # clamp if bigger than parent's edge
    parent_node <- getParent(tree, node_number)
    parent_edge_length <- tree$edge.length[ which(tree$edge[,2] == node_number) ]
    if (length(parent_edge_length) == 0) {
      new_branch_length <- 0
    } else if (new_branch_length > parent_edge_length) {
      cat(sprintf(
        "Clamping branch for %s from %.2f to 0 (polytomy)\n",
        species_name, new_branch_length
      ))
      new_branch_length <- 0
    }
  }
  
  updated_tree <- bind.tip(
    tree,
    tip.label = species_name,
    where     = node_number,
    position  = new_branch_length
  )
  
  cat("Grafted", species_name, "onto", rank_found, "=", tax_row[[rank_found]], 
      "; used rate =", node_rate, "-> BL =", new_branch_length, "\n")
  
  return(updated_tree)
}

## Now the loop that adds the new species
updated_tree <- existing_tree  # Start from your labeled tree

for (i in seq_len(nrow(new_tax_df))) {
  sp      <- new_tax_df$Species[i]
  tax_row <- new_tax_df[i, c("Genus","Family","Order","Class","Phylum")]
  
  if (is.null(updated_tree$root.edge)) {
    updated_tree$root.edge <- 0
  }
  
  updated_tree <- graft_species_fallback(
    tree               = updated_tree,
    species_name       = sp,
    tax_row            = tax_row,
    node_label_data    = node_label_data,
    all_clade_rates    = all_clade_rates,
    default_rate       = 0.01,
    taxonomy_existing_df = taxonomy_existing_df   # pass the existing tips' taxonomy
  )
}

cat("Now the updated_tree has", length(updated_tree$tip.label), "tips.\n")

###############################################################################
## STEP 8: (Optional) Re-ultrametrize
###############################################################################

# If desired, you can re-ultrametrize the final tree to ensure all root-to-tip
# distances are equal
updated_tree <- force.ultrametric(updated_tree, method = "nnls")

###############################################################################
## STEP 9: Save & Visualize
###############################################################################

write.tree(updated_tree, file = "final_grafted_tree.newick")
