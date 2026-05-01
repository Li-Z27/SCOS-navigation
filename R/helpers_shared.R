# =============================================================================
# helpers_shared.R
# Shared helper functions for SCOS navigation, two-stage learning 
# =============================================================================

########################################################################
#####                 GRAPH CONSTRUCTION                           #####
########################################################################
# create lattice vertex grid
Lattice_Vertices <- function(x, y) {
  # INPUT
  # x: grid width  (vertices from 0 to x)
  # y: grid height (vertices from 0 to y)
  #
  # OUTPUT
  # data.frame with columns Var1 (x coords), Var2 (y coords)
  # row i corresponds to vertex i in the graph
  
  expand.grid(0:x, 0:y)
}

# create discretized graph
Graph_Discretized <- function(x, y){
  # INPUT
  # x: grid width  (vertices from 0 to x)
  # y: grid height (vertices from 0 to y)
  
  n_v <- (x + 1) * (y + 1)
  
  # vertex index helper — converts (i,j) grid coords to vertex index
  # i in 0:x, j in 0:y
  vid <- function(i, j) 1 + i + j * (x + 1)
  
  # pre-allocate edge lists
  # max edges: horizontal + vertical + diagonal
  # horizontal: x*(y+1), vertical: (x+1)*y, diagonal: x*y*2
  max_edges <- x * (y + 1) + (x + 1) * y + 2 * x * y
  from_v  <- integer(max_edges)
  to_v    <- integer(max_edges)
  weights <- numeric(max_edges)
  ptr     <- 1
  
  for (j in 0:y) {
    for (i in 0:x) {
      v <- vid(i, j)
      
      # horizontal edge (i, j) -> (i+1, j)
      if (i < x) {
        from_v[ptr]  <- v
        to_v[ptr]    <- vid(i + 1, j)
        weights[ptr] <- 1
        ptr <- ptr + 1
      }
      
      # vertical edge (i, j) -> (i, j+1)
      if (j < y) {
        from_v[ptr]  <- v
        to_v[ptr]    <- vid(i, j + 1)
        weights[ptr] <- 1
        ptr <- ptr + 1
      }
      
      # diagonal edge (i, j) -> (i+1, j+1)
      if (i < x && j < y) {
        from_v[ptr]  <- v
        to_v[ptr]    <- vid(i + 1, j + 1)
        weights[ptr] <- sqrt(2)
        ptr <- ptr + 1
        
        # diagonal edge (i+1, j) -> (i, j+1)
        from_v[ptr]  <- vid(i + 1, j)
        to_v[ptr]    <- vid(i, j + 1)
        weights[ptr] <- sqrt(2)
        ptr <- ptr + 1
      }
    }
  }
  
  # trim to actual size
  actual <- ptr - 1
  from_v  <- from_v[seq_len(actual)]
  to_v    <- to_v[seq_len(actual)]
  weights <- weights[seq_len(actual)]
  
  # build graph directly from edge list
  G <- make_empty_graph(n = n_v, directed = FALSE)
  G <- set_vertex_attr(G, "name", value = as.character(seq_len(n_v)))
  G <- add_edges(G, rbind(from_v, to_v))
  G <- set_edge_attr(G, "weight", value = weights)
  
  return(G)
}

# compute edge-obstacle intersection matrix
Intersect_Info <- function(G_original, x, y, obstacle_df, r) {
  # INPUT
  # G_original  : igraph object
  # x, y        : grid dimensions
  # obstacle_df : obstacle information (x, y, disambiguation_cost, prob, status)
  # r           : obstacle radius
  
  edge_df    <- get.data.frame(G_original, what = "edges")
  n_edges    <- nrow(edge_df)
  n_obs      <- nrow(obstacle_df)
  vertice_list <- Lattice_Vertices(x, y)
  
  intersect_mat <- matrix(0, nrow = n_edges, ncol = n_obs)
  
  for (i in seq_len(n_obs)) {
    
    # skip obstacles with zero disambiguation cost — they never block
    if (obstacle_df$disambiguation_cost[i] == 0) next
    
    cx <- obstacle_df$x[i]
    cy <- obstacle_df$y[i]
    
    # bounding box around obstacle — limits vertex search
    x_min <- floor(cx - r - sqrt(2))
    x_max <- ceiling(cx + r + sqrt(2))
    y_min <- floor(cy - r - sqrt(2))
    y_max <- ceiling(cy + r + sqrt(2))
    
    # find vertices near obstacle boundary — two rings:
    # ring 1 (outside): r < dist <= r + sqrt(2)
    # ring 2 (inside) : r - sqrt(2) < dist <= r
    bbox_idx <- which(vertice_list[, 1] >= x_min &
                        vertice_list[, 1] <= x_max &
                        vertice_list[, 2] >= y_min &
                        vertice_list[, 2] <= y_max)
    
    if (length(bbox_idx) == 0) next
    
    dx_v  <- vertice_list[bbox_idx, 1] - cx
    dy_v  <- vertice_list[bbox_idx, 2] - cy
    dist2 <- dx_v^2 + dy_v^2
    r2    <- r^2
    r_out <- (r + sqrt(2))^2
    r_in  <- (r - sqrt(2))^2
    
    # outer ring: just outside obstacle boundary
    outer_ring <- bbox_idx[dist2 > r2 & dist2 <= r_out]
    # inner ring: just inside obstacle boundary
    inner_ring <- bbox_idx[dist2 > r_in & dist2 <= r2]
    
    if (length(outer_ring) == 0 || length(inner_ring) == 0) next
    
    # find edges crossing obstacle boundary:
    # edges connecting one outer vertex to one inner vertex
    # with length 1 or sqrt(2)
    for (o in outer_ring) {
      for (inn in inner_ring) {
        dx_e <- vertice_list[o, 1] - vertice_list[inn, 1]
        dy_e <- vertice_list[o, 2] - vertice_list[inn, 2]
        dist_e <- dx_e^2 + dy_e^2
        
        # edge length must be 1 or sqrt(2)
        if (dist_e == 1 || dist_e == 2) {
          edge_idx <- get.edge.ids(G_original, c(o, inn))
          if (edge_idx > 0) {
            intersect_mat[edge_idx, i] <- 1
          }
        }
      }
    }
  }
  
  return(intersect_mat)
}

# fnd interior vertices
find_vertices_in_obstacle <- function(obstacle_center, r, vertice_list){
  # INPUT
  # obstacle_cernter: dataframe row with x, y coordinates
  # r: obstacle radius
  # vertice_list: lattice vertice coordiates (x, y columns)
  
  x_center <- as.numeric(obstacle_center[1])
  y_center <- as.numeric(obstacle_center[2])
  
  # bounding box 
  x_min <- ceiling(x_center - r)
  x_max <- floor(x_center + r)
  y_min <- ceiling(y_center - r)
  y_max <- floor(y_center + r)
  bbox_idx <- which(vertice_list[, 1] >= x_min &
                      vertice_list[, 1] <= x_max &
                      vertice_list[, 2] >= y_min &
                      vertice_list[, 2] <= y_max)
  
  if (length(bbox_idx) == 0){
    return(integer(0))
  }
  
  # find distance
  dx <- vertice_list[bbox_idx, 1] - x_center
  dy <- vertice_list[bbox_idx, 2] - y_center
  inside <- bbox_idx[dx^2 + dy^2 <= r^2]
  
  # OUTPUT
  return(inside)
}

########################################################################
#####                  BELIEF UPDATE                               #####
########################################################################
# prior covariance matrix 
create_prior_cov <- function(obstacle_df, sigma_f, l) {
  # INPUT:
  # obstacle_df : obstacle information dataframeo
  # sigma_f     : marginal standard deviation 
  # l           : length scale 
  #
  # KERNEL: k(xi, xj) = sigma_f^2 * exp(-||xi-xj||^2 / (2*l^2))
  
  coords <- as.matrix(obstacle_df[, c("x", "y")])
  
  diff_x <- outer(coords[, 1], coords[, 1], "-")
  diff_y <- outer(coords[, 2], coords[, 2], "-")
  dist2  <- diff_x^2 + diff_y^2   # squared distance — no sqrt needed
  
  cov_matrix       <- sigma_f^2 * exp(-dist2 / (2 * l^2))
  diag(cov_matrix) <- sigma_f^2 + 0.5   # nugget = 0.5
  
  return(cov_matrix)
}

# get neighbors list
get_neighbors_list <- function(obstacle_df, d_max){
  # INPUT:
  # obstacle_df: obstacle_info dataframe (x, y, disambiguation_cost, prob, status)
  # d_max      : maximum distance threshold
  
  coords <- as.matrix(obstacle_df[, c("x", "y")])
  
  # vectorised pairwise distance matrix — no loop
  diff_x   <- outer(coords[, 1], coords[, 1], "-")
  diff_y   <- outer(coords[, 2], coords[, 2], "-")
  dist_mat <- sqrt(diff_x^2 + diff_y^2)
  
  # for each obstacle, find indices within d_max
  lapply(seq_len(nrow(coords)), function(i) {
    which(dist_mat[i, ] <= d_max)
  })
}

# distance from state to all obstacles 
calculate_distances_to_state <- function(s_current, obstacle_df) {
  # INPUT
  # s_current  : coordinate vector c(x, y) of current position
  # obstacle_df: obstacle information
  
  sqrt((obstacle_df$x - s_current[1])^2 +
         (obstacle_df$y - s_current[2])^2)
}

# transformation
logit <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)  # Ensure in (0, 1) range
  return(log(p / (1 - p)))
}
inv_log <- function(log_p){
  return(exp(log_p)/(1+exp(log_p)))
}

# posterior update
posterior_update <- function(obstacle_df, Record_obstacle, 
                             prior_mean, prior_cov,
                             observation,
                             noise_variance, lambda){
  # INPUT
  # obstacle_df: obstacle information
  # x, y
  # disambiguation_cost
  # prob
  # status
  # Record_obstacle: vector of obstacle status - TRUE once disambiguated
  # prior_mean, prior_cov: n-vector of prior log-odds means and n*n covariance matrix
  # observation: n-vector of beta probability readings, NA - not sensed
  # noise_variance: n-vector of per-obstacle noise variance sigma^2
  # lambda: beta distribution parameter
  
  n <- length(prior_mean)
  eps <- 1e-8
  
  # --- STEP 1 Build Observations ------------------------------
  log_LH <- log((dbeta(observation, 4+lambda, 4-lambda)+eps)/
                  (dbeta(observation, 4-lambda, 4+lambda)+eps))
  obs_logodds <- prior_mean +  log_LH
  
  # --- STEP 2 Dimension Handling ------------------------------
  # NA observation -> infinite variance / 0 precision 
  dim_missing <- which(is.na(observation))
  noise_variance[dim_missing] <- Inf
  obs_logodds[dim_missing] <- prior_mean[dim_missing]
  
  # active dimensions
  dim_dis <- which(Record_obstacle) # disambiguated ones
  if (length(dim_dis) > 0) {
    obs_logodds[dim_dis]    <- ifelse(obstacle_df$status[dim_dis] == 1, 5, -5)
    noise_variance[dim_dis] <- 1e-4  # near-zero noise = hard observation
  }
  
  # all dimensions active 
  dim_active  <- seq_len(n)
  n_active    <- n
  
  mu_active    <- prior_mean
  K_active     <- prior_cov
  y_active     <- obs_logodds
  sigma2_active <- noise_variance
  
  # --- STEP 3 Posterior ---------------------------------------
  M <- K_active + diag(sigma2_active, n_active, n_active)
  diag(M) <- diag(M)+eps
  
  cholM <- tryCatch(
    chol(M),
    error = function(e) {
      M_retry <- M
      diag(M_retry) <- diag(M_retry) + 1e-4
      chol(M_retry)
    }
  )
  
  # back-substitution (avoid inversion)
  Minv_y <- backsolve(cholM, forwardsolve(t(cholM), y_active))
  Minv_K <- backsolve(cholM, forwardsolve(t(cholM), K_active))
  
  posterior_mean_active <- as.numeric(K_active %*% Minv_y)
  posterior_cov_active <- K_active - K_active %*% Minv_K
  
  # --- STEP 4 Restore full length -----------------------------
  posterior_mean_full <- posterior_mean_active
  if (length(dim_dis) > 0) {
    posterior_mean_full[dim_dis] <- ifelse(obstacle_df$status[dim_dis] == 1, 5, -5)
  }
  
  posterior_cov_full <- as.matrix(posterior_cov_active)
  
  # OUTPUT
  output <- list(pos_mean = posterior_mean_full,
                 pos_cov = posterior_cov_full)
  
  return(output)
}

# sense and update
sense_and_update <- function(s_current, vertice_list,
                             obstacle_df, intersect_mat, Record_obstacle, New_obstacle, Reduce_obstacle,
                             mean_prior, cov_prior, noise_var,
                             sensing_range, sense_noise){
  # INPUT:
  # CAT I
  # s_current: current vertex name
  # vertice_list: lattice vertex coordinates
  
  # CAT II
  # obstacle_df: obstacle information (x, y, disambiguation_cost, prob, status)
  # intersect_mat: edge-obstacle intersection matrix
  # Record_obstacle: logical vector, TRUE once disambiguated
  # New_obstacle: index of obstacle just disambiguated
  # Reduce_obstacle: indices of obstacles confirmed as always-blocking
  
  # CAT III
  # mean_prior: current prior mean vector (log-odds)
  # cov_prior: current prior covariance matrix
  # noise_var: per-obstacle noise variance scalar
  # sensing_range: sensor detection radius
  # sense_noise: lambda parameter for Beta observation noise
  
  # --- STEP 1: process disambiguation from previous step ---------------
  if (length(New_obstacle) > 0) {
    obstacle_df$disambiguation_cost[New_obstacle] <- 0
    obstacle_df$prob[New_obstacle]                <- obstacle_df$status[New_obstacle]
    
    # remove intersection for false obstacles
    false_new <- New_obstacle[obstacle_df$status[New_obstacle] == FALSE]
    if (length(false_new) > 0) {
      intersect_mat[, false_new] <- 0
    }
    Record_obstacle[New_obstacle] <- TRUE
  }
  
  # pin confirmed blocking obstacles
  if (length(Reduce_obstacle) > 0) {
    obstacle_df$prob[Reduce_obstacle] <- 1
  }
  
  # --- STEP 2: sense obstacles within range ----------------------------
  s_coor <- as.numeric(vertice_list[s_current, ])
  dx     <- obstacle_df$x - s_coor[1]
  dy     <- obstacle_df$y - s_coor[2]
  sense_index <- which(dx^2 + dy^2 <= sensing_range^2 & !Record_obstacle) 
  
  # --- STEP 3: build observation vector --------------------------------
  observation <- rep(NA_real_, nrow(obstacle_df))
  
  if (length(sense_index) > 0) {
    # vectorised Beta sampling based on true status
    true_idx  <- sense_index[obstacle_df$status[sense_index] == TRUE]
    false_idx <- sense_index[obstacle_df$status[sense_index] == FALSE]
    
    if (length(true_idx) > 0) {
      observation[true_idx]  <- rbeta(length(true_idx),
                                      4 + sense_noise, 4 - sense_noise)
    }
    if (length(false_idx) > 0) {
      observation[false_idx] <- rbeta(length(false_idx),
                                      4 - sense_noise, 4 + sense_noise)
    }
  }
  
  # --- STEP 4: posterior update ----------------------------------------
  # only run if there are new sensor readings
  if (length(sense_index) > 0 || length(New_obstacle) > 0) {
    noise_variance <- rep(noise_var, nrow(obstacle_df))
    post           <- posterior_update(obstacle_df, Record_obstacle,
                                       mean_prior, cov_prior,
                                       observation,
                                       noise_variance, sense_noise)
    mean_prior <- post$pos_mean
    cov_prior  <- post$pos_cov
    
    # --- STEP 5: update probabilities ----------------------------------
    prob_pos   <- as.numeric(inv_log(mean_prior))
    unresolved <- which(!Record_obstacle)
    obstacle_df$prob[unresolved] <- prob_pos[unresolved]
    
    # re-pin confirmed blocking obstacles after prob update
    if (length(Reduce_obstacle) > 0) {
      obstacle_df$prob[Reduce_obstacle] <- 1
    }
  }
  
  # OUTPUT
  output <- list(mean_prior      = mean_prior,
                 cov_prior       = cov_prior,
                 obstacle_df     = obstacle_df,
                 intersect_mat   = intersect_mat,
                 Record_obstacle = Record_obstacle)
  return(output)
  
}

# sampling 
prediction_sample <- function(obstacle_df, 
                              pos_mean_log_odds, pos_cov){
  # INPUT:
  # obstacle_df: obstacle information (x, y, disambiguation_cost, prob, status)
  # pos_mean_log_odds: posterior mean vector (log-odds scale)
  # pos_cov: posterior covariance matrix (for log-odds scale)
  
  # obstacles with fixed status — no sampling needed
  adjust_index <- which(obstacle_df$prob == 0 | obstacle_df$prob == 1)
  sample_index <- setdiff(seq_len(nrow(obstacle_df)), adjust_index)
  
  # initialise pseudo_status
  obstacle_df$pseudo_status <- NA
  
  # fixed obstacles — set directly from prob
  if (length(adjust_index) > 0) {
    obstacle_df$pseudo_status[adjust_index] <- obstacle_df$prob[adjust_index]
  }
  
  # ambiguous obstacles — sample from posterior distribution
  # sampling from full posterior captures uncertainty and spatial correlation
  if (length(sample_index) > 0) {
    sampled_log_odds <- mvrnorm(1,
                                mu    = pos_mean_log_odds[sample_index],
                                Sigma = pos_cov[sample_index,
                                                sample_index,
                                                drop = FALSE])
    sampled_probs <- as.numeric(inv_log(sampled_log_odds))
    obstacle_df$pseudo_status[sample_index] <- rbinom(length(sampled_probs), 1, sampled_probs)
  }
  
  return(obstacle_df)
}

########################################################################
#####                 CANDIDATE SET GENERATION                     #####
########################################################################
# compute candidate set
candidate_set <- function(G_original, obstacle_df, intersect_mat,
                          edge_df, vertice_list, 
                          g, s_current, cutoff_opt){
  # INPUT
  # CAT I - Graph Related
  # G_original: igraph object (edge weights = Euclidean lengths)
  # obstacle_df: obstacle_information (x, y, disambiguation_cost, prob, status)
  # intersect_mat: edge-obstacle intersection matrix
  
  # CAT II - Edge/Vertex Related
  # edge_df: pre-computed get.data.frame(G_original, what="edges")
  # vertice_list: lattice vertex coordinates
  
  # CAT III - Start & Objective
  # g: index of goal vertex
  # s_current: current vertex index
  
  # CAT IV - Blockage Probability
  # cutoff_opt: true (prob > cutoff_opt, default = 0.99) blockage based on posterior belief
  
  # --- STEP 1 optimistic belief & graph construct ----------------------
  pseudo_status <- as.numeric(obstacle_df$prob >= cutoff_opt)
  local_weights <- edge_df$weight
  
  # block confirmed obstacles and adjust edge weights
  true_obs <- which(pseudo_status == 1)
  if(length(true_obs) > 0){
    blocked_edges <- which(rowSums(intersect_mat[,true_obs,drop = FALSE]) > 0)
    local_weights[blocked_edges] <- Inf
  }
  
  # add disambiguation costs for ambiguous obstacles
  false_obs <- which(pseudo_status == 0)
  if(length(false_obs) > 0){
    cost_additions <- intersect_mat[,false_obs, drop=F] %*% (0.5*obstacle_df$disambiguation_cost[false_obs])
    local_weights <- local_weights + as.numeric(cost_additions)
  }
  
  # --- STEP 2 initialization -------------------------------------------
  state_set      <- c()
  obstacle_set   <- c()
  cost_immediate <- c()
  path_record    <- list()
  H              <- c()
  itr            <- 1
  
  s_idx <- which(V(G_original)$name == as.character(s_current))
  g_idx <- which(V(G_original)$name == as.character(g))
  
  # --- STEP 3 main loop ------------------------------------------------
  while (TRUE){
    # --- 3(a) get path and track ---------------------------------------
    path_result <- shortest_paths(G_original,
                                  from     = s_idx,
                                  to       = g_idx,
                                  weights  = local_weights,
                                  output   = "both")
    
    V_list <- as.numeric(attributes(path_result$vpath[[1]])$names)
    E_list <- as.numeric(path_result$epath[[1]])
    
    # record optimistic path costs at this iteration
    H[itr] <- sum(local_weights[E_list])
    
    # no path exists or already at goal
    if(length(V_list) == 0 || s_current == g){
      state_set          <- c(state_set, g)
      obstacle_set       <- c(obstacle_set, NA)
      cost_immediate     <- c(cost_immediate, 0)
      path_record[[itr]] <- c()
      break
    }
    
    # --- 3(b) follow path, find obstacle -------------------------------
    seg_len       <- 0
    free_path     <- c()
    s_next        <- NULL
    obstacle_next <- NULL
    
    for(i in seq(2, length(V_list))){
      v1 <- V_list[i-1]
      v2 <- V_list[i]
      
      # edge lookup using vertex name
      edge_idx <- get.edge.ids(G_original, c(v1, v2))
      
      obs_hit <- which(intersect_mat[edge_idx,] == 1)
      if(length(obs_hit) > 0){
        # multiple obstacles on same edge - pick closest
        if(length(obs_hit) > 1){
          v1_coords <- as.numeric(vertice_list[v1,])
          dists <- rowSums((as.matrix(obstacle_df[obs_hit, c("x","y")]) - 
                              matrix(v1_coords, nrow = length(obs_hit), ncol = 2, byrow = T))^2)
          obs_hit <- obs_hit[which.min(dists)]
        }
        s_next <- v1 
        obstacle_next <- obs_hit
        break
      }
      # cumulate path cost and path vertices
      seg_len <- seg_len + edge_df$weight[edge_idx]
      free_path <- c(free_path, v2)
    }
    
    # --- 3(c): record candidate ----------------------------------------
    # path reaches goal with no obstacle hit
    if(is.null(s_next)){
      state_set          <- c(state_set,g)
      obstacle_set       <- c(obstacle_set, NA)
      cost_immediate     <- c(cost_immediate, seg_len)
      path_record[[itr]] <- free_path
      break
    }
    
    # add candidate and record
    state_set          <- c(state_set, s_next)
    obstacle_set       <- c(obstacle_set, obstacle_next)
    cost_immediate     <- c(cost_immediate, seg_len)
    path_record[[itr]] <- free_path
    
    # temporarily block this obstacle's edges for next iteration
    obs_edges <- which(intersect_mat[,obstacle_next] > 0)
    local_weights[obs_edges] <- Inf
    
    itr <- itr + 1
  }
  
  # OUTPUT
  output <- list(state_set      = state_set,
                 obstacle_set   = obstacle_set,
                 cost_immediate = cost_immediate,
                 path_record    = path_record,
                 H              = H)
  return(output)
}

# search space reduction
reduce_candidate_set <- function(G_original, obstacle_df, intersect_mat,
                                 edge_df, vertice_list,
                                 g, s_current, 
                                 candidate_output, r, cutoff_opt){
  # INPUT
  # CAT I - Graph Related
  # G_original      : igraph object
  # obstacle_df     : obstacle information (x, y, disambiguation_cost, prob, status)
  # intersect_mat   : edge-obstacle intersection matrix
  #
  # CAT II - Edge/Vertex Related
  # edge_df         : pre-computed get.data.frame(G_original, what="edges")
  # vertice_list    : lattice vertex coordinates
  #
  # CAT III - Start & Objective
  # g               : goal vertex name
  # s_current       : current vertex name
  #
  # CAT IV - Candidate Set & Radius
  # candidate_output: output from candidate_set()
  # r               : obstacle radius
  #
  # CAT V - Blockage Probability
  # cutoff_opt: true (prob > cutoff_opt, default = 0.99) blockage based on posterior belief
  
  # --- STEP 1 Extract Candidate Set Info --------------------------------
  state_set <- candidate_output$state_set
  obstacle_set <- candidate_output$obstacle_set
  cost_immediate <- candidate_output$cost_immediate
  
  upper_bound <- cost_immediate[length(cost_immediate)]
  
  # --- STEP 2 Identify Obstacle to Evaluate -----------------------------
  candidate_obs <- unique(obstacle_set[!is.na(obstacle_set)])
  true_obs <- which(obstacle_df$prob == 1)
  check_obs <- setdiff(seq_len(nrow(obstacle_df)),union(true_obs, candidate_obs))
  
  if(length(check_obs) == 0){
    return(integer(0))
  }
  
  # --- STEP 3 Build Optimistic Weight Vector ----------------------------
  pseudo_status <- as.numeric(obstacle_df$prob >= cutoff_opt)
  local_weights <- edge_df$weight
  
  true_obs_idx <- which(pseudo_status == 1)
  if(length(true_obs_idx) > 0){
    blocked_edges <- which(rowSums(intersect_mat[, true_obs_idx, drop = F]) > 0)
    local_weights[blocked_edges] <- Inf
  }
  
  false_obs_idx <- which(pseudo_status == 0)
  if(length(false_obs_idx) > 0){
    cost_additions <- intersect_mat[, false_obs_idx, drop = F] %*% (0.5*obstacle_df$disambiguation_cost[false_obs_idx])
    local_weights <- local_weights +  as.numeric(cost_additions)
  }
  
  # --- STEP 4 Initialization --------------------------------------------
  obstacle_remove <- c()
  
  # --- STEP 5 Main Loop -------------------------------------------------
  for (obs in check_obs){
    # --- 5(a) highest probability obs -----------------------------------
    if (pseudo_status[obs] == 1){
      obstacle_remove <- c(obstacle_remove, obs)
      next
    }
    # --- 5(b) pseudo vertex & graph construction ------------------------
    # find interior vertices of obstacle
    vertex_inside <- find_vertices_in_obstacle(obstacle_df[obs, c("x", "y")],
                                               r, vertice_list)
    if (length(vertex_inside) == 0) {
      obstacle_remove <- c(obstacle_remove, obs)
      next
    }
    
    # pseudo vertex
    new_edges <- data.frame(from   = c(vertex_inside, rep("p", length(vertex_inside))),
                            to     = c(rep("p", length(vertex_inside)), vertex_inside),
                            weight = rep(0, 2*length(vertex_inside)))
    edge_df_p <- rbind(edge_df[, c("from", "to", "weight")], new_edges)
    edge_df_p$weight[seq_len(nrow(edge_df))] <- local_weights
    
    # graph with pseudo vertex
    G_p <- graph.data.frame(edge_df_p, directed= F)
    p_idx   <- which(V(G_p)$name == "p")
    s_idx_p <- which(V(G_p)$name == as.character(s_current))
    g_idx_p <- which(V(G_p)$name == as.character(g))
    
    # --- 5(c) lower bound -----------------------------------------------
    # lower bound part 1
    path_1 <- shortest_paths(G_p,
                             from    = s_idx_p,
                             to      = p_idx,
                             weights = edge_df_p$weight,
                             output  = "both")
    E_list_1 <- as.numeric(path_1$epath[[1]])
    V_names_1 <- attributes(path_1$vpath[[1]])$names
    
    if(length(E_list_1) == 0){
      obstacle_remove <- c(obstacle_remove, obs)
      next
    }
    
    cost_1 <- sum(edge_df_p$weight[E_list_1])
    
    # lower bound part 2
    v_entry <- V_names_1[length(V_names_1) - 1]
    v_entry_idx <- which(V(G_p)$name == v_entry)
    path_2 <- shortest_paths(G_p,
                             from    = v_entry_idx,
                             to      = g_idx_p,
                             weights = edge_df_p$weight,
                             output  = "epath")
    E_list_2 <- as.numeric(path_2$epath[[1]])
    
    if (length(E_list_2) == 0) {
      obstacle_remove <- c(obstacle_remove, obs)
      next
    }
    
    cost_2      <- sum(edge_df_p$weight[E_list_2])
    lower_bound <- cost_1 + cost_2
    
    # check remove condition
    if(lower_bound >= upper_bound){
      obstacle_remove <- c(obstacle_remove, obs)
    }
  }
  # OUTPUT
  output <- obstacle_remove
  return(output)
}

########################################################################
#####                 INFORMATION GAIN                             #####
########################################################################
# information gain - simulation
compute_info_gain_approx <- function(neighbors_list, Record_obstacle,
                                     cov_prior, noise_variance, 
                                     alpha, I_sum){
  # INPUT:
  # neighbors_list: pre-computed list of neighbouring obstacle indices
  #                 (list, for all obstacles)
  # Record_obstacle: vector indicating if disambiguated
  # cov_prior: current posterior covariance matrix
  # noise_variance : per-obstacle observation noise variances sigma^2
  # alpha: scaling factor for info gain
  # I_sum: cumulative info gathered so far
  
  vapply(seq_along(Record_obstacle), function(i) {
    
    # skip already disambiguated obstacles
    if (Record_obstacle[i]) return(0)
    
    # keep only unresolved neighbours
    nb_valid <- neighbors_list[[i]][!Record_obstacle[neighbors_list[[i]]]]
    
    if (length(nb_valid) == 0) return(0)
    
    # I(Y_t) = 0.5 * log|I + Sigma^-1 K_sub|
    n <- length(nb_valid)
    K_sub <- cov_prior[nb_valid, nb_valid, drop = FALSE]
    Sigma_inv  <- diag(1 / noise_variance[nb_valid], n, n)
    I_t   <- 0.5 * as.numeric(determinant(diag(n) + Sigma_inv %*% K_sub,
                                          logarithm = TRUE)$modulus)
    
    alpha * (sqrt(I_t + I_sum) - sqrt(I_sum))
    
  }, numeric(1))
}

# information gain - implementation
compute_info_gain_exact <- function(state, vertice_list,
                                    obstacle_df, Record_obstacle,
                                    cov_prior, noise_variance, 
                                    sensing_range, I_sum, alpha){
  # INPUT:
  # state: stopping vertex
  # vertice_list: lattice vertex coordinates
  # obstacle_df: obstacle information df
  # Record_obstacle: vector indicating if disambiguated
  # cov_prior: current posterior covariance matrix
  # noise_variance : per-obstacle observation noise variances sigma^2
  # sensing_range  : sensor detection radius (to obstacle center)
  # I_sum: cumulative info gathered so far
  # alpha: scaling factor for info gain
  
  # --- STEP 1 obstacles in the range ---------------------------
  s_coor <- as.numeric(vertice_list[state,])
  dx <- obstacle_df$x - s_coor[1]
  dy <- obstacle_df$y - s_coor[2]
  obs_index <- which(dx^2 + dy^2 <= sensing_range^2 & !Record_obstacle)
  # early stop condition
  if(length(obs_index) == 0){
    output <- list(gain_info = 0,
                   I_sum     =  I_sum)
    return(output)
  }
  
  # --- STEP 2 info gain compute --------------------------------
  n         <- length(obs_index)
  K_sub     <- cov_prior[obs_index, obs_index, drop = FALSE]
  Sigma_inv <- diag(1 / noise_variance[obs_index], n, n)
  I_t       <- 0.5 * as.numeric(determinant(diag(n) + Sigma_inv %*% K_sub,
                                            logarithm = TRUE)$modulus)
  gain      <- alpha * (sqrt(I_t + I_sum) - sqrt(I_sum))  
  I_sum     <- I_sum + I_t
  
  output <- list(gain_info = gain,
                 I_sum     = I_sum)
  return(output)
}

########################################################################
#####                 SIMULATION (MC)                              #####
########################################################################
# simulate trajectory
simulate_trajectory <- function(s_start, g, obs_start,
                                G_original, n_vertices, obstacle_df, intersect_mat,
                                edge_df, vertice_list,
                                V_states, V_function,
                                r, epsilon,
                                method_decision, G_values){
  # INPUT:
  # CAT I: Start and Goal Vertex
  # s_start         : starting vertex 
  # g               : goal vertex 
  # obs_start       : first obstacle index to process (NA if none)
  
  # CAT II: Graph Related
  # G_original   : igraph object
  # n_vertices   : total number of vertices in graph
  # obstacle_df  : local copy with pseudo_status already sampled
  # intersect_mat: edge-obstacle intersection matrix (local copy)
  
  # CAT III: Edge/Vertice Related
  # edge_df      : pre-computed get.data.frame(G_original, what="edges")
  # vertice_list : lattice vertex coordinates
  
  # CAT IV: Agent Related
  # V_states        : vector of states with value estimates
  # V_function      : associated value estimates
  # r               : obstacle radius
  # epsilon         : exploration rate for epsilon-greedy
  # method_decision : 'greedy'/'epsilon-greedy'/'softmax'
  # G_values        : information gain bonus per obstacle
  
  # --- STEP 1: Process starting obstacles ------------------------------
  if(!is.na(obs_start)){
    obstacle_df$prob[obs_start] <- as.numeric(obstacle_df$pseudo_status[obs_start])
    obstacle_df$disambiguation_cost[obs_start] <- 0
    if(obstacle_df$pseudo_status[obs_start] == 0){
      intersect_mat[,obs_start] <- 0
    }
  }
  
  # --- STEP 2: Initialization ------------------------------------------
  # value lookup
  V_lookup           <- numeric(n_vertices)
  V_lookup[V_states] <- V_function
  
  # track
  max_steps      <- 500
  traj_states    <- integer(max_steps)
  traj_costs     <- numeric(max_steps)
  traj_states[1] <- s_start
  step           <- 1
  s_current      <- s_start
  
  # --- STEP 3: Main loop -----------------------------------------------
  while (s_current != g && step < max_steps){
    # --- 3(a) candidate set --------------------------------------------
    cs <- candidate_set(G_original    = G_original, 
                        obstacle_df   = obstacle_df[,c("x", "y", "disambiguation_cost", "prob", "status")], 
                        intersect_mat = intersect_mat,
                        edge_df       = edge_df, 
                        vertice_list  = vertice_list, 
                        g             = g, 
                        s_current     = s_current,
                        cutoff_opt    = 1)
    
    set_state_sim    <- cs$state_set
    set_obstacle_sim <- cs$obstacle_set
    cost_imm_sim     <- cs$cost_immediate
    n_cand           <- length(set_state_sim)
    
    if(n_cand == 1 && set_state_sim[1] == g){
      # early stop condition 
      step <- step + 1
      traj_states[step] <- g
      traj_costs[step-1] <- cost_imm_sim[1]
      break
    }
    
    # --- 3(b) add disambiguation cost and info gain ----------------
    for (i in seq_len(n_cand - 1)){
      obs_i <- set_obstacle_sim[i]
      if(!is.na(obs_i)){
        cost_imm_sim[i] <- cost_imm_sim[i] + obstacle_df$disambiguation_cost[obs_i] - G_values[obs_i]
      }
    }
    
    # --- 3(c) add future value function --------------------------------
    future_val <- vapply(set_state_sim, 
                         function(sv){
                           if (sv >= 1 && sv <= n_vertices) V_lookup[sv] else 0
                         }, numeric(1))
    future_val[n_cand] <- 0  # goal state always has 0 value
    
    decision_val <- cost_imm_sim + future_val
    
    # --- 3(d) method_decision ------------------------------------------
    idx_chosen <- switch (method_decision,
                          "greedy" = {
                            which.min(decision_val)
                          },
                          "epsilon-greedy" = {
                            if (runif(1) < epsilon) {
                              # explore: choose uniformly from non-greedy candidates
                              greedy_idx    <- which.min(decision_val)
                              non_greedy    <- setdiff(seq_len(n_cand), greedy_idx)
                              if (length(non_greedy) > 0) sample(non_greedy, 1)
                              else greedy_idx   # only one candidate, no choice
                            } else {
                              # exploit: choose greedily
                              which.min(decision_val)
                            }
                          },
                          "softmax" = {
                            ev   <- -decision_val
                            ev   <- ev - max(ev)   # numerical stability
                            prob <- exp(ev) / sum(exp(ev))
                            sample(n_cand, 1, prob = prob)
                          }
    )
    
    s_next_sim <- set_state_sim[idx_chosen]
    obs_next_sim <- set_obstacle_sim[idx_chosen]
    
    # --- 3(e) record step and update info ------------------------------
    step <- step + 1
    traj_states[step] <- s_next_sim
    traj_costs[step-1] <- cost_imm_sim[idx_chosen]
    
    if(!is.na(obs_next_sim)){
      obstacle_df$prob[obs_next_sim] <- as.numeric(obstacle_df$pseudo_status[obs_next_sim])
      obstacle_df$disambiguation_cost[obs_next_sim] <- 0
      if(obstacle_df$pseudo_status[obs_next_sim] == 0){
        intersect_mat[, obs_next_sim] <- 0
      }
    }
    s_current <- s_next_sim
  }
  
  # OUTPUT - information gain adjusted trajectory history
  used <- seq_len(step)
  output <- list(trajectory_states = traj_states[used],
                 trajectory_costs  = traj_costs[seq_len(step - 1)])
  return(output)
}
