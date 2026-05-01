# =============================================================================
# helpers_DRL.R
# Additional helper functions for SCOS navigation, two-stage learning (DRL) 
# =============================================================================

# lower bound of path cost from s (current vertex) to g (optimistic)
lower_bound <- function(G_original, obstacle_df, 
                        intersect_mat, edge_df, 
                        s, g){
  # INPUT
  # CAT I 
  # G_original: igraph object
  # obstacle_df: obstacle information (x, y, disambiguation_cost, prob, status)
  # CAT II
  # intersect_mat: edge-obstacle intersection matrix
  # edge_df: pre-computed get.data.frame(G_original, what="edges")
  # CAT III
  # s: current vertex
  # g: goal vertex
  
  local_weights <- edge_df$weight
  
  # confirmed blockage
  true_obs <- which(obstacle_df$prob == 1)
  if(length(true_obs) > 0){
    blocked_edges <- which(rowSums(intersect_mat[, true_obs, drop = F]) > 0)
    local_weights[blocked_edges] <- Inf
  }
  # ambiguous obstacles
  false_obs <- which(obstacle_df$prob != 1)
  if(length(false_obs) > 0){
    cost_additions <- intersect_mat[, false_obs, drop = FALSE] %*%
      (0.5 * obstacle_df$disambiguation_cost[false_obs])
    local_weights  <- local_weights + as.numeric(cost_additions)
  }
  
  # shortest path under optimistic weights
  s_idx <- which(V(G_original)$name == as.character(s))
  g_idx <- which(V(G_original)$name == as.character(g))
  
  path   <- shortest_paths(G_original,
                           from    = s_idx,
                           to      = g_idx,
                           weights = local_weights,
                           output  = "epath")
  
  E_list <- as.numeric(path$epath[[1]])
  
  # OUTPUT
  # lower bound for traversal cost under optimistic assumption
  if(length(E_list) == 0){
    return(Inf)
  } else {
    return(sum(local_weights[E_list]))
  }
}


# upper bound of path cost from s (current vertex) to g (pessimistic)
upper_bound <- function(G_original, obstacle_df, 
                        intersect_mat, edge_df, 
                        s, g){
  # INPUT
  # CAT I 
  # G_original: igraph object
  # obstacle_df: obstacle information (x, y, disambiguation_cost, prob, status)
  # CAT II
  # intersect_mat: edge-obstacle intersection matrix
  # edge_df: pre-computed get.data.frame(G_original, what="edges")
  # CAT III
  # s: current vertex
  # g: goal vertex
  
  local_weights <- edge_df$weight
  
  # pessimistic belief — block everything except confirmed free (prob==0)
  block_obs <- which(obstacle_df$prob != 0)
  if (length(block_obs) > 0) {
    blocked_edges <- which(rowSums(intersect_mat[, block_obs, drop = FALSE]) > 0)
    local_weights[blocked_edges] <- Inf
  }
  
  # shortest path under pessimistic weights
  s_idx <- which(V(G_original)$name == as.character(s))
  g_idx <- which(V(G_original)$name == as.character(g))
  
  path   <- shortest_paths(G_original,
                           from    = s_idx,
                           to      = g_idx,
                           weights = local_weights,
                           output  = "epath")
  
  E_list <- as.numeric(path$epath[[1]])
  
  # OUTPUT
  # Upper bound of traversal cost under pessimistic assumption
  if (length(E_list) == 0){
    return(Inf)
  } else {
    return(sum(local_weights[E_list]))
  }
}

# initialize value distribution for a new state
init_state_value <- function(s, g, obstacle_df, intersect_mat,
                             G_original, edge_df, replace_condi) {
  # INPUT
  # s            : vertex to initialise
  # g            : goal vertex
  # obstacle_df  : obstacle information
  # intersect_mat: edge-obstacle intersection matrix
  # G_original   : igraph object
  # edge_df      : pre-computed edge data frame
  # replace_condi: support grid spacing delta
  
  # goal state — deterministic zero cost
  if (s == g) {
    return(list(V_function  = 0,
                V_prob      = 1,
                V_indicator = TRUE))
  }
  # find bounds
  lower <- lower_bound(G_original, obstacle_df, intersect_mat, edge_df, s, g)
  upper <- upper_bound(G_original, obstacle_df, intersect_mat, edge_df, s, g)
  # grid
  step      <- min(replace_condi, max(upper - lower, replace_condi))
  grid      <- seq(lower, upper, by = step)
  grid_size <- length(grid)
  
  # OUTPUT
  output <- list(V_function  = grid,
                 V_prob      = rep(1, grid_size),
                 V_indicator = c(TRUE, rep(FALSE, max(0, grid_size - 2)), TRUE))
  
  return(output)
}

# find where in the support grid each observed cost falls
find_neighbors_indices <- function(val_current, bell_temp){
  # INPUT
  # val_current: sorted numeric vector — support grid of V(s)
  # bell_temp: numeric vector — shifted observed costs to interpolate
  
  # OUTPUT
  # list of length(bell_temp), each element is a list with:
  #   $indices: integer vector of 1 or 2 neighboring support indices
  #   $weights: numeric vector of interpolation weights (sums to 1)
  
  n     <- length(val_current)
  n_obs <- length(bell_temp)
  
  lapply(seq_len(n_obs), function(k) {
    bv <- bell_temp[k]
    
    # exact match — assign full weight to that point
    exact <- which(val_current == bv)
    if (length(exact) > 0) {
      return(list(indices = exact[1],
                  weights = 1))
    }
    
    # below minimum — assign full weight to lowest support point
    if (bv <= val_current[1]) {
      return(list(indices = 1,
                  weights = 1))
    }
    
    # above maximum — assign full weight to highest support point
    if (bv >= val_current[n]) {
      return(list(indices = n,
                  weights = 1))
    }
    
    # interior — find bracketing interval [x_j, x_{j+1}]
    j       <- findInterval(bv, val_current)   # index of x_j
    x_j     <- val_current[j]
    x_j1    <- val_current[j + 1]
    width   <- x_j1 - x_j
    # interpolation weights from projected Bellman update
    w_j  <- (x_j1 - bv) / width   # weight for lower neighbour
    w_j1 <- (bv - x_j) / width    # weight for upper neighbour
    return(list(indices = c(j, j + 1),
                weights = c(w_j, w_j1)))
  })
}

# support refinement
update_support <- function(V_function_s, V_prob_s, V_indicator_s,
                           val_observe, replace_condi){
  # INPUT
  # V_function_s: support grid for state s (sorted numeric vector)
  # V_prob_s : Dirichlet parameters for state s
  # V_indicator_s : observed/boundary flags for state s
  # val_observe : new observed cumulative cost
  # replace_condi : neighbourhood threshold delta
  # OUTPUT
  # list(V_function_s, V_prob_s, V_indicator_s) — updated support
  
  # exact match - increment count directly
  exact_idx <- which(V_function_s == val_observe)
  if (length(exact_idx) > 0) {
    V_prob_s[exact_idx[1]] <- V_prob_s[exact_idx[1]] + 1
    return(list(V_function  = V_function_s,
                V_prob      = V_prob_s,
                V_indicator = V_indicator_s))
  }
  
  n <- length(V_function_s)
  
  # find immediate left and right neighbours
  j <- findInterval(val_observe, V_function_s)
  
  # left neighbour: index j (0 if val_observe below minimum)
  # right neighbour: index j+1 (n+1 if val_observe above maximum)
  left_idx  <- if (j >= 1) j      else NA_integer_
  right_idx <- if (j <  n)  j + 1 else NA_integer_
  
  # check which immediate neighbours are unobserved and within replace_condi
  replace_idx <- c()
  if (!is.na(left_idx) &&
      !V_indicator_s[left_idx] &&
      abs(V_function_s[left_idx] - val_observe) <= replace_condi) {
    replace_idx <- c(replace_idx, left_idx)
  }
  if (!is.na(right_idx) &&
      !V_indicator_s[right_idx] &&
      abs(V_function_s[right_idx] - val_observe) <= replace_condi) {
    replace_idx <- c(replace_idx, right_idx)
  }
  
  if (length(replace_idx) == 0) {
    # no unobserved immediate neighbours within delta — insert as new point
    V_function_s  <- c(V_function_s,  val_observe)
    V_prob_s      <- c(V_prob_s,      1)
    V_indicator_s <- c(V_indicator_s, TRUE)
  } else {
    # replace immediate unobserved neighbours — absorb their counts
    alpha_base    <- sum(V_prob_s[replace_idx]) + 1
    V_function_s  <- V_function_s[-replace_idx]
    V_prob_s      <- V_prob_s[-replace_idx]
    V_indicator_s <- V_indicator_s[-replace_idx]
    
    V_function_s  <- c(V_function_s,  val_observe)
    V_prob_s      <- c(V_prob_s,      alpha_base)
    V_indicator_s <- c(V_indicator_s, TRUE)
  }
  
  # maintain sorted order
  ord           <- order(V_function_s)
  V_function_s  <- V_function_s[ord]
  V_prob_s      <- V_prob_s[ord]
  V_indicator_s <- V_indicator_s[ord]
  
  return(list(V_function  = V_function_s,
              V_prob      = V_prob_s,
              V_indicator = V_indicator_s))
  
}

# simulate - DRL approach (offline)
simulate_trajectory_DRL <- function(s_start, g, obs_start,
                                    obstacle_df, intersect_mat,
                                    G_original, edge_df, vertice_list,
                                    V_states, V_function, V_prob,
                                    V_indicator, N_states,
                                    replace_condi, G_values){
  # INPUT
  # CAT I: Start and Goal
  # s_start  : starting vertex name
  # g        : goal vertex name
  # obs_start: first obstacle index to process (NA if none)
  #
  # CAT II: Graph Related
  # G_original   : igraph object
  # obstacle_df  : local copy with pseudo_status already sampled
  # intersect_mat: edge-obstacle intersection matrix (local copy)
  # edge_df      : pre-computed get.data.frame(G_original, what="edges")
  # vertice_list : lattice vertex coordinates
  #
  # CAT III: Value Function (distributional)
  # V_states   : vertex indices with value estimates
  # V_function : list of support grids per state
  # V_prob     : list of Dirichlet parameters per state
  # V_indicator: list of observed/boundary flags per state
  # N_states   : visit counts per state
  #
  # CAT IV: Agent
  # replace_condi: support refinement threshold delta
  # G_values     : info gain bonus per obstacle 
  
  # --- STEP 1 Process starting obstacle ---------------------------------------
  if (!is.na(obs_start)){
    obstacle_df$prob[obs_start]                <- as.numeric(obstacle_df$pseudo_status[obs_start])
    obstacle_df$disambiguation_cost[obs_start] <- 0
    if (obstacle_df$pseudo_status[obs_start] == 0) {
      intersect_mat[, obs_start] <- 0
    }
  }
  
  # --- STEP 2 Initialization tracking -----------------------------------------
  max_steps      <- 500
  traj_states    <- integer(max_steps)
  traj_costs     <- numeric(max_steps)
  traj_states[1] <- s_start
  step           <- 1
  s_current      <- s_start
  
  # --- STEP 3 Main loop -------------------------------------------------------
  while (s_current != g && step < max_steps){
    # --- 3(a) candidate set ---------------------------------------------------
    cs <- candidate_set(G_original    = G_original,
                        obstacle_df   = obstacle_df,
                        intersect_mat = intersect_mat,
                        edge_df       = edge_df,
                        vertice_list  = vertice_list,
                        g             = g,
                        s_current     = s_current,
                        cutoff_opt    = 1)
    
    Set_state    <- cs$state_set
    Set_obstacle <- cs$obstacle_set
    Cost_imm     <- cs$cost_immediate
    n_cand       <- length(Set_state)
    
    # --- 3(b) early stop if only goal -----------------------------------------
    if (n_cand == 1 && Set_state[1] == g) {
      step                  <- step + 1
      traj_states[step]     <- g
      traj_costs[step - 1] <- Cost_imm[1]
      break
    }
    
    # --- 3(c) add disambiguation cost and info gain ---------------------------
    non_goal     <- seq_len(n_cand - 1)
    obs_non_goal <- Set_obstacle[non_goal]
    valid        <- !is.na(obs_non_goal)
    
    if (any(valid)) {
      vi <- non_goal[valid]
      oi <- obs_non_goal[valid]
      Cost_imm[vi] <- Cost_imm[vi] +
        obstacle_df$disambiguation_cost[oi] -
        G_values[oi]
    }
    
    # --- 3(d) sampling — greedy under sampled values -----------------
    # known states   : sample from posterior
    # unknown states : lower_bound
    # goal state     : always 0
    future_val <- vapply(Set_state, function(sv) {
      idx <- match(sv, V_states)
      if (!is.na(idx)) {
        # known — Thompson sampling
        as.numeric(rdirichlet(1, V_prob[[idx]])) %*% V_function[[idx]]
      } else {
        # unknown — midpoint estimate
        lower_bound(G_original, obstacle_df, intersect_mat, edge_df, sv, g)
      }
    }, numeric(1))
    future_val[n_cand] <- 0 # g is always zero
    
    # greedy selection under sampled value estimates
    idx_chosen <- which.min(Cost_imm + future_val)
    s_next     <- Set_state[idx_chosen]
    obs_next   <- Set_obstacle[idx_chosen]
    step_cost  <- Cost_imm[idx_chosen]
    
    # --- 3(e) initialise s_next if first visit --------------------------------
    idx_next <- match(s_next, V_states)
    if (is.na(idx_next)) {
      init        <- init_state_value(s_next, g, obstacle_df,
                                      intersect_mat, G_original,
                                      edge_df, replace_condi)
      V_states    <- c(V_states,   s_next)
      V_function  <- append(V_function,  list(init$V_function))
      V_prob      <- append(V_prob,      list(init$V_prob))
      V_indicator <- append(V_indicator, list(init$V_indicator))
      N_states    <- c(N_states,   0)
      idx_next    <- length(V_states)
    }
    
    # --- 3(f) distributional Bellman update -----------------------------------
    idx_current  <- match(s_current, V_states)
    prob_next    <- V_prob[[idx_next]]
    prob_next    <- prob_next / sum(prob_next) 
    
    bell_temp    <- V_function[[idx_next]] + step_cost
    prob_current <- V_prob[[idx_current]]
    neighbors    <- find_neighbors_indices(V_function[[idx_current]],bell_temp)
    
    for (k in seq_along(neighbors)) {
      nb                       <- neighbors[[k]]
      prob_current[nb$indices] <- prob_current[nb$indices] + prob_next[k] * nb$weights
    }
    
    V_prob[[idx_current]] <- prob_current
    N_states[idx_current] <- N_states[idx_current] + 1
    
    # --- 3(g) record step -----------------------------------------------------
    step                 <- step + 1
    traj_states[step]    <- s_next
    traj_costs[step - 1] <- step_cost
    
    # --- 3(h) update local obstacle_df and intersect_mat ----------------------
    if (!is.na(obs_next)) {
      obstacle_df$prob[obs_next]                <- as.numeric(obstacle_df$pseudo_status[obs_next])
      obstacle_df$disambiguation_cost[obs_next] <- 0
      if (obstacle_df$pseudo_status[obs_next] == 0) {
        intersect_mat[, obs_next] <- 0
      }
    }
    
    s_current <- s_next
  }
  
  # --- STEP 4: Support refinement (given full trajectory) ---------------------
  n_steps   <- step - 1
  traj_s    <- traj_states[seq_len(step)]
  traj_c    <- traj_costs[seq_len(n_steps)]
  costs_cum <- rev(cumsum(rev(traj_c)))
  
  for (i in seq_len(n_steps)) {
    idx               <- match(traj_s[i], V_states)
    upd               <- update_support(V_function[[idx]], V_prob[[idx]],
                                        V_indicator[[idx]], costs_cum[i],
                                        replace_condi)
    V_function[[idx]]  <- upd$V_function
    V_prob[[idx]]      <- upd$V_prob
    V_indicator[[idx]] <- upd$V_indicator
  }
  
  # OUTPUT
  output <- list(trajectory_states = traj_s,
                 trajectory_costs  = traj_c,
                 V_states          = V_states,
                 V_function        = V_function,
                 V_prob            = V_prob,
                 V_indicator       = V_indicator,
                 N_states          = N_states)
  return(output)
}

# simulate trajectory with DRL approach (online: no support refinement, 
# which will be done all together after online iterations)
simulate_trajectory_DRL_rollout <- function(s_start, g, obs_start,
                                            obstacle_df, intersect_mat,
                                            G_original, edge_df, vertice_list,
                                            V_states, V_function, V_prob,
                                            V_indicator, N_states,
                                            replace_condi, G_values){
  # INPUT
  # CAT I: Start and Goal
  # s_start  : starting vertex name
  # g        : goal vertex name
  # obs_start: first obstacle index to process (NA if none)
  #
  # CAT II: Graph Related
  # G_original   : igraph object
  # obstacle_df  : local copy with pseudo_status already sampled
  # intersect_mat: edge-obstacle intersection matrix (local copy)
  # edge_df      : pre-computed get.data.frame(G_original, what="edges")
  # vertice_list : lattice vertex coordinates
  #
  # CAT III: Value Function (distributional)
  # V_states   : vertex indices with value estimates
  # V_function : list of support grids per state
  # V_prob     : list of Dirichlet parameters per state
  # V_indicator: list of observed/boundary flags per state
  # N_states   : visit counts per state
  #
  # CAT IV: Agent
  # replace_condi: support refinement threshold delta
  # G_values     : info gain bonus per obstacle
  
  # --- STEP 1 Process starting obstacle ---------------------------------------
  if (!is.na(obs_start)){
    obstacle_df$prob[obs_start]                <- as.numeric(obstacle_df$pseudo_status[obs_start])
    obstacle_df$disambiguation_cost[obs_start] <- 0
    if (obstacle_df$pseudo_status[obs_start] == 0) {
      intersect_mat[, obs_start] <- 0
    }
  }
  
  # --- STEP 2 Initialization tracking -----------------------------------------
  max_steps      <- 500
  traj_states    <- integer(max_steps)
  traj_costs     <- numeric(max_steps)
  traj_states[1] <- s_start
  step           <- 1
  s_current      <- s_start
  
  # --- STEP 3 Main loop -------------------------------------------------------
  while (s_current != g && step < max_steps){
    # --- 3(a) candidate set ---------------------------------------------------
    cs <- candidate_set(G_original    = G_original,
                        obstacle_df   = obstacle_df,
                        intersect_mat = intersect_mat,
                        edge_df       = edge_df,
                        vertice_list  = vertice_list,
                        g             = g,
                        s_current     = s_current,
                        cutoff_opt    = 1)
    
    Set_state    <- cs$state_set
    Set_obstacle <- cs$obstacle_set
    Cost_imm     <- cs$cost_immediate
    n_cand       <- length(Set_state)
    
    # --- 3(b) early stop if only goal -----------------------------------------
    if (n_cand == 1 && Set_state[1] == g) {
      step                  <- step + 1
      traj_states[step]     <- g
      traj_costs[step - 1] <- Cost_imm[1]
      break
    }
    
    # --- 3(c) add disambiguation cost and info gain ---------------------------
    non_goal     <- seq_len(n_cand - 1)
    obs_non_goal <- Set_obstacle[non_goal]
    valid        <- !is.na(obs_non_goal)
    
    if (any(valid)) {
      vi <- non_goal[valid]
      oi <- obs_non_goal[valid]
      Cost_imm[vi] <- Cost_imm[vi] +
        obstacle_df$disambiguation_cost[oi] -
        G_values[oi]
    }
    
    # --- 3(d) Thompson sampling — greedy under sampled values -----------------
    # known states   : sample from posterior
    # unknown states : lower_bound
    # goal state     : always 0
    future_val <- vapply(Set_state, function(sv) {
      idx <- match(sv, V_states)
      if (!is.na(idx)) {
        # known — Thompson sampling
        as.numeric(rdirichlet(1, V_prob[[idx]])) %*% V_function[[idx]]
      } else {
        # unknown — midpoint estimate
        lower_bound(G_original, obstacle_df, intersect_mat, edge_df, sv, g)
      }
    }, numeric(1))
    future_val[n_cand] <- 0 # g is always zero
    
    # greedy selection under sampled value estimates
    idx_chosen <- which.min(Cost_imm + future_val)
    s_next     <- Set_state[idx_chosen]
    obs_next   <- Set_obstacle[idx_chosen]
    step_cost  <- Cost_imm[idx_chosen]
    
    # --- 3(e) initialise s_next if first visit --------------------------------
    idx_next <- match(s_next, V_states)
    if (is.na(idx_next)) {
      init        <- init_state_value(s_next, g, obstacle_df,
                                      intersect_mat, G_original,
                                      edge_df, replace_condi)
      V_states    <- c(V_states,   s_next)
      V_function  <- append(V_function,  list(init$V_function))
      V_prob      <- append(V_prob,      list(init$V_prob))
      V_indicator <- append(V_indicator, list(init$V_indicator))
      N_states    <- c(N_states,   0)
      idx_next    <- length(V_states)
    }
    
    # --- 3(f) distributional Bellman update -----------------------------------
    idx_current  <- match(s_current, V_states)
    prob_next    <- V_prob[[idx_next]]
    prob_next    <- prob_next / sum(prob_next) 
    
    bell_temp    <- V_function[[idx_next]] + step_cost
    prob_current <- V_prob[[idx_current]]
    neighbors    <- find_neighbors_indices(V_function[[idx_current]],bell_temp)
    
    for (k in seq_along(neighbors)) {
      nb                       <- neighbors[[k]]
      prob_current[nb$indices] <- prob_current[nb$indices] + prob_next[k] * nb$weights
    }
    
    V_prob[[idx_current]] <- prob_current
    N_states[idx_current] <- N_states[idx_current] + 1
    
    # --- 3(g) record step -----------------------------------------------------
    step                 <- step + 1
    traj_states[step]    <- s_next
    traj_costs[step - 1] <- step_cost
    
    # --- 3(h) update local obstacle_df and intersect_mat ----------------------
    if (!is.na(obs_next)) {
      obstacle_df$prob[obs_next]                <- as.numeric(obstacle_df$pseudo_status[obs_next])
      obstacle_df$disambiguation_cost[obs_next] <- 0
      if (obstacle_df$pseudo_status[obs_next] == 0) {
        intersect_mat[, obs_next] <- 0
      }
    }
    
    s_current <- s_next
  }
  
  # OUTPUT
  n_steps   <- step - 1
  output <- list(trajectory_states = traj_states[seq_len(step)],
                 trajectory_costs  = traj_costs[seq_len(n_steps)],
                 V_states          = V_states,
                 V_function        = V_function,
                 V_prob            = V_prob,
                 V_indicator       = V_indicator,
                 N_states          = N_states)
  return(output)
}

