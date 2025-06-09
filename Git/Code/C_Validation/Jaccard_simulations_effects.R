## Random Jaccard simualtions based on absolute numbers (not proportions) ##
##
## number and identity of sites occupied in 1 (number: from 1 to max_vector_length, identitiy: ID of the occupied sites that was assigned to vector elements from 1:max_vector_length)
## number of sites occupied in 2 (number and identity, similar to above)
## number of flips (from 1 to max_vector_length) (should be linear i.e., not exponential)
## type of flip: either: "random" or "from_1_to_0"
## I want to simulate a set of inital vector based on different parameter combinations,
## modify them based on num_of_flips and type of flip,
## extract IDs for elements that have now a 1, where IDs (i.e, index of vector elements) = set vectors 1 and 2


set.seed(123) # For reproducibility



simulate_jaccard <- function(lengths, occu1, flip_types, n_runs = 1) {
  # Expand grid for parameter combinations
  params <- expand.grid(lengths, occu1, flip_types, stringsAsFactors = FALSE)
  colnames(params) <- c("length", "occu1", "flip_type")

  # Initialize storage
  jaccard_results <- list()

  # Loop through parameter combinations
  for (row_i in seq_len(nrow(params))) {
    length_i <- params$length[row_i]
    occu_1_i <- params$occu1[row_i]
    flip_type <- params$flip_type[row_i]

    # Ensure the number of occupied sites does not exceed the vector length
    occu_1_i <- min(occu_1_i, length_i)

    # Repeat the simulation `n_runs` times
    for (run in seq_len(n_runs)) {
      # Generate set1 (initially occupied sites)
      set1_i <- sample(1:length_i, occu_1_i, replace = FALSE)

      # Sequence of number of changes (from 1 to max_vector_length)
      n_changes <- seq(from = 1, to = length_i, by = 1)

      # Dataframe to store results for this run
      jaccard_temp <- data.frame(num_changes = n_changes, jaccard = NA, run = run)

      set2_temp_list <- vector("list", length(n_changes))

      # Modify set1 stepwise and store set2 states
      for (n_i in seq_along(n_changes)) {
        set2_temp <- set1_i
        changes_i <- n_changes[n_i]

        if (flip_type == "random") {
          flip_sites <- sample(1:length_i, size = changes_i, replace = FALSE)

          set2_i <- union(                    # Combine left overs after 1 to 0 flips with new cells from 0 to 1 flips
            setdiff(set2_temp, flip_sites),   # Remove indices in `flip_sites` already in `set2_temp`
            intersect(flip_sites, setdiff(1:length_i, set2_temp)) # Add indices not in `set2_temp` (0 to 1; colonizations)
          )

        } else if (flip_type == "from_1_to_0") {
          sample_mat <- if (length(set1_i) == 1) rep(set1_i, 2) else set1_i
          flip_sites <- sample(x = sample_mat, size = min(changes_i, length(set1_i)), replace = FALSE)
          set2_i <- setdiff(set2_temp, flip_sites) # Remove flipped sites from set1
          set2_i <- if (length(set2_i) == 0) NA else set2_i
        }

        # Store modified set2
        set2_temp_list[[n_i]] <- set2_i

        # Compute Jaccard similarity
        jaccard_res <- my_jaccard_abc(set1_i, set2_i)

        jaccard_temp$jaccard[n_i] <- jaccard_res$j
        jaccard_temp$intersect[n_i] <- jaccard_res$intersect
        jaccard_temp$unique_1[n_i] <- jaccard_res$unique_1
        jaccard_temp$unique_2[n_i] <- jaccard_res$unique_2
        jaccard_temp$N_occu_1[n_i] <- jaccard_res$N_occu_1
        jaccard_temp$N_occu_2[n_i] <- jaccard_res$N_occu_2
        jaccard_temp$length[n_i] <- length_i
        jaccard_temp$occu1[n_i] <- occu_1_i
        jaccard_temp$flip_type[n_i] <- flip_type
      }

      # Append results for this run
      jaccard_results[[length(jaccard_results) + 1]] <- jaccard_temp
    }
  }

  # Combine Jaccard results into a single dataframe for visualization
  jaccard_df <- do.call(rbind, jaccard_results)
  return(jaccard_df)
}

# Example usage
lengths <- c(10, 100, 1000)
occu1 <- c(2^(0:floor(log2(1000))), 999)
flip_types <- c("random", "from_1_to_0")
n_runs <- 100  # Repeat each parameter combination 10 times



jaccard_df <- simulate_jaccard(lengths, occu1, flip_types, n_runs)
saveRDS(jaccard_df, "Code/Sensitivity/Jaccard_Prevalence_Simulations/res/jaccards_sim_absolute_100runs.rds")


