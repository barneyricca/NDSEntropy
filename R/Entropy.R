census <- function(data_vec, codes) {
  rep(0, length(codes)) -> census_vec
  codes -> names(census_vec)

  for(i in seq_along(data_vec)) {
    census_vec[data_vec[i]] + 1 -> census_vec[data_vec[i]]
  }
  return(census_vec)
}

win_ent <- function(data_seq,
                    parts = 16,
                    lag_max = 50,
                    alpha = 0.05) {
  # Follows Wiltshire, Fiore, & Butner (2017), sort of
  # Returns a data frame (or NULL, if failed) of two time series:
  #  1. A (numeric) windowed entropy
  #  2. A marker of maxima (TRUE for a maxima, FALSE for not)
  #
  require(tseriesChaos)

  if(mode(data_seq) == "integer") {
    as.character(data_seq) -> data_seq
  }

  if(mode(data_seq) == "numeric") {
    max(data_seq) -> max_data
    min(data_seq) -> min_data
    as.character(1:parts) -> codes
    ceiling( parts *
               (data_seq - min_data) /
               (max_data - min_data)) -> partition
    1 -> partition[which(partition == 0)]  # partition goes 0:parts, not 1:parts
    codes[partition] -> data_seq
  }

  if(mode(data_seq) != "character") {
    return(NULL)
  }

  require(data.table)
  unique(data_seq) -> codes
  length(codes) -> n_codes

  # Because some sequences may be short
  n_codes -> parts
  min(lag_max, length(data_seq)-n_codes) -> lag_max

  table(data_seq) -> code_counts
  1:n_codes -> nums
  codes -> names(nums)

  code_counts / sum(code_counts) -> code_probs

  # The appropriate lag, according to Fraser & Swinney, 1986, is at
  #  the first minimum of mutual information.
  if(n_codes > 1) {
    which.min(mutual(unname(nums[data_seq]),
                     partitions = n_codes,       # Arbitrary partitions needed                                                     #  for continuous data; any
                     #  partition >= n_codes works
                     #  for categorical data.
                     lag.max = lag_max,
                     plot = FALSE)) -> wind_sz
  } else {
    return(NULL)
  }

  # Calculate the windowed entropy, ent[]:
  # This is different than what Wiltshire, Fiore, & Butner (2017) do. However,
  #  I think that they overestimate the issue (and they probably would think
  #  that I underestimate the issue.) They use the natural log, and only use
  #  the probabilities within the window; it is the latter that I find problematic.
  # However, there are also problems with this approach, so look at this more carefully.
  #
  # They also do some smoothing; that may be helpful
  length(data_seq) -> seq_len
  rep(0.0, length = (seq_len - wind_sz + 1)) -> ent

  # Old way
  #  for(i in 1:(length(ent))) {
  #    for(j in 0:(wind_sz-1)) {
  #      ent[i] - (code_probs[data_seq[i+j]] *
  #                  log2(code_probs[data_seq[i+j]])) -> ent[i]
  #    }
  #  }

  # A faster way, it appears, than the old way
  for(i in 1:(length(ent))) {
    ent[i] - sum(code_probs[data_seq[i:(i+wind_sz-1)]] *
                   log2(code_probs[data_seq[i:(i+wind_sz-1)]])) -> ent[i]
  }

  # Find the appropriate maxima:
  # Pass 1 - find all maxima:
  vector("logical", length(ent)) -> maxes  # All entries are FALSE by default

  if(length(ent) > 2) {
    for (i in 2:(length(ent)-1)) {
      if(ent[i] >= ent[i-1] &&
         ent[i] >= ent[i+1]) {
        TRUE -> maxes[i]
      }
    }
  } else {
    return(NULL)
  }

  # Pass 2 - which maxima indicate code probability changes from before to after a
  #  maximum located in Pass 1
  1 -> i_old
  length(data_seq) -> n_max

  for(i in 1:length(maxes)) {
    if(maxes[i] == TRUE) {
      0 -> cont_tab
      if((i+1) < (n_max - n_codes - wind_sz)) {   # Skip if too close to the end
        census(data_seq[i_old:i], codes) -> before
        census(data_seq[(i+1):n_max], codes) -> after
        rbind(before,after) -> cont_tab

        # The fisher.test() is appropriate for the next one, but it sometimes
        #  runs out of memory. The chi-squared sometimes fails for really
        #  small cell values (and zeroes); hence the is.na().
        #  The chisq.test() is also very biased, and it will generate warnings if any
        #  cell is < 5. Hence, use the fisher.test() is any cell size < 10.
        #  Probably should just use the chisq.test() if all(cont_tab > 100); will do
        #  that on the 3rd pass.
        if (sum(cont_tab) > 0) {
          if (any(cont_tab < 10)) {
            fisher.test(cont_tab) -> p_temp
          } else {
            chisq.test(cont_tab) -> p_temp
          }
          if (is.na(p_temp$p.value)) {
            FALSE -> maxes[i]
          } else {
            # The suppressWarnings seems to be necessary for about
            #  a quarter of the vector entries. I don't know why yet:
            if(p_temp$p.value > 0.10) { # 0.1 is arbitrary
              FALSE -> maxes[i]  # No change
            } else {
              i -> i_old         # Only include info from this sub-segment
            }
          }
        } else {
          FALSE -> maxes[i]
        }
      }
    }
  }
  # In principle, this refinement should be iterated until it converges. However,
  #  as we know this approach is biased, we won't try to refine something that
  #  will still be a bit biased.

  # However, let's do a 3rd pass:
  #  Fisher test for probably all, but subset by subset, with updating.

  which(maxes == TRUE) -> dividing_points
  length(dividing_points) -> num_possibles
  if(num_possibles > 0) {

    1 -> begin_a
    length(data_seq) -> end_b
    dividing_points[1] -> end_a -> begin_b
    census(data_seq[begin_a:end_a], codes) -> before
    census(data_seq[begin_b:end_b], codes) -> after
    rbind(before, after) -> cont_tab

    if(any(cont_tab < 100)) {
      fisher.test(cont_tab)$p.value -> p_value
    } else {
      chisq.test(cont_tab)$p.value -> p_value
    }

    if(p_value > alpha) {
      FALSE -> maxes[dividing_points[1]]
      1 -> dividing_points[1]
    }

    if (num_possibles > 1) {
      TRUE -> continue
      2 -> index
      while (continue) {
        dividing_points[index - 1] -> begin_a
        dividing_points[index] -> begin_b -> end_a
        census(data_seq[begin_a:end_a], codes) -> before
        census(data_seq[begin_b:end_b], codes) -> after
        rbind(before, after) -> cont_tab
        if (any(cont_tab < 100)) {
          if(length(which(colSums(cont_tab)>0))>1) {
            fisher.test(cont_tab,
                        simulate.p.value = TRUE)$p.value -> p_value
          } else {
            0 -> p_value
          }
        } else {
          chisq.test(cont_tab)$p.value -> p_value
        }
        if (p_value > alpha) {
          FALSE -> maxes[dividing_points[index]]
          dividing_points[-index] -> dividing_points
          index - 1 -> index
        }
        index + 1 -> index
        if (index > length(dividing_points)) {
          FALSE -> continue
        }
      }
    }

    #  if (length(which(maxes == TRUE)) > 1) {
    #    length(data_seq) -> end_b
    #    max(which(maxes == TRUE)) -> begin_b -> end_a
    #    if (length(which(maxes == TRUE)) > 2) {
    #      which(maxes == TRUE) -> temp_maxes
    #      temp_maxes[length(temp_maxes) - 1] -> begin_a
    #    } else {
    #      1 -> begin_a
    #    }
    #  } else {
    #    1 -> begin_a
    #    length(data_seq) -> end_b
    #    dividing_points[num_possibilities] -> begin_b -> end_a
    #  }

    #  unname(table(data_seq[begin_a:end_a])) -> before
    #  unname(table(data_seq[begin_b:end_b])) -> after
    #  rbind(before, after) -> cont_tab
    #  if (any(cont_tab < 100)) {
    #    fisher.test(cont_tab)$p.value -> p_value
    #  } else {
    #    chisq.test(cont_tab)$p.value -> p_value
    #  }
    #  if (p_value > alpha) {
    #    FALSE -> maxes[dividing_points[num_possibiles]]
    #  }
  }

  data.frame("Windowed_Entropy" = ent,    # Entropy
             "Maxima" = maxes) -> ent_df  # True if a local maximum of entropy
  return(ent_df)
}
