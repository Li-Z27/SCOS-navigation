# =============================================================================
# main_TS.R
# Two-stage policy (MC)
# =============================================================================
Policy_TS <- function(obstacle_df, x, y, r, s, g, status, 
                      method_decision, sensing_range, sense_noise,
                      noise_var_scale, l, sigma_f,
                      alpha_ig, div_IG){
  # INPUT
  # --- CAT I: Graph-related ----
  # obstacle_df: obstacle information
    # x_coor, y_coor
    # disambiguation_cost
    # prob
  # x, y: grid dimension
  # r: obstacle radius
  # s, g: start and goal vertex
  # status: ground truth status vector for all obstacles
  
  # --- CAT II: Agent-related ---
  # method_decision: greedy/epsilon-greedy/softmax
  # sensing_range: the range in which the sensor can get prob reading
  # sense_noise: lambda value (higher lambda, higher precision level)
  
  # --- CAT III: Belief-related ---
  # noise_var_scale: scale used based on variance from delta method
  # l: length scale
  # sigma_f: marginal standard deviation
  
  # --- CAT IV: InformationGain-related ---
  # alpha_ig: initial information gain scaling factor
  # div_IG: scaling factor adjusted based on observed values
  
  # --- STEP 0 Hyperparameters --------------------------------------
  epsilon           <- 0.1   # exploration rate for epsilon-greedy
  n_rollouts        <- 30    # online rollouts per decision step (adjusted based on problem size)
  N_itr_max         <- 2000  # max offline MC iterations
  conv_tol          <- 0.005 # convergence tolerance
  conv_patience     <- 10    # consecutive under-tolerance iters before stopping
  
  # --- STEP 1 Graph Construction -----------------------------------
  obstacle_df$status <- status
  obstacle_df        <- subset(obstacle_df, select = -prob_true)
  
  vertice_list       <- Lattice_Vertices(x, y)
  G_original         <- Graph_Discretized(x, y)
  intersect_mat      <- Intersect_Info(G_original, x, y, obstacle_df, r)
  
  edge_df            <- get.data.frame(G_original, what = "edges")
  n_vertices         <- vcount(G_original)
  
  # --- STEP 2 Belief Initialization --------------------------------
  noise_var      <- (64 / ((4+sense_noise) * (4-sense_noise) * 9))*noise_var_scale
  mean_prior     <- rep(0, nrow(obstacle_df))
  cov_prior      <- create_prior_cov(obstacle_df, sigma_f, l)
  noise_variance <- rep(noise_var, nrow(obstacle_df))
  # neighbors list for approximation info gain
  neighbors_list <- get_neighbors_list(obstacle_df, sensing_range + r)
  
  # --- STEP 3 Initialization (tracking) ----------------------------
  s_current       <- s
  Record_obstacle <- rep(FALSE, nrow(obstacle_df)) # TRUE once disambiguated
  Reduce_obstacle <- integer(0) # obstacles confirmed as always-blocking/unnecessary
  New_obstacle    <- integer(0) # obstacle just disambiguated at recent step
  
  V_states   <- integer(0) # visited states with estimated values
  V_function <- numeric(0) # associated estimated values
  N_states   <- integer(0) # visit counts per state
  
  I_sum        <- 0 # cumulative information
  Total_length <- 0 # path length (Euclidean)
  Total_cost   <- 0 # disambiguation cost accumulated
  Record_path  <- integer(0) # ordered vertex indices of executed path

  idx_step     <- 1
  
  # --- STEP 4 Main Loop --------------------------------------------
  while(s_current != g){
    # stop once reach the goal vertex
    
    # --- 4(a) Sense & update ---------------------------------------
    belief <- sense_and_update(s_current       = s_current,
                               vertice_list    = vertice_list,
                               obstacle_df     = obstacle_df,
                               intersect_mat   = intersect_mat,
                               Record_obstacle = Record_obstacle,
                               New_obstacle    = New_obstacle,
                               Reduce_obstacle = Reduce_obstacle,
                               mean_prior      = mean_prior,
                               cov_prior       = cov_prior,
                               noise_var       = noise_var,
                               sensing_range   = sensing_range,
                               sense_noise     = sense_noise)
    
    mean_prior      <- belief$mean_prior
    cov_prior       <- belief$cov_prior
    obstacle_df     <- belief$obstacle_df
    intersect_mat   <- belief$intersect_mat
    Record_obstacle <- belief$Record_obstacle
    
    # --- 4(b) Decision set & reduction -----------------------------
    cs <- candidate_set(G_original    = G_original,
                        obstacle_df   = obstacle_df,
                        intersect_mat = intersect_mat,
                        edge_df       = edge_df,
                        vertice_list  = vertice_list,
                        g             = g,
                        s_current     = s_current,
                        cutoff_opt    = 0.99)
    
    Set_state    <- cs$state_set
    Set_obstacle <- cs$obstacle_set
    Cost_imm     <- cs$cost_immediate
    Path_record  <- cs$path_record
    H            <- cs$H
    
    # reduction
    rcs <- reduce_candidate_set(G_original       = G_original,
                                obstacle_df      = obstacle_df,
                                intersect_mat    = intersect_mat,
                                edge_df          = edge_df,
                                vertice_list     = vertice_list,
                                g                = g,
                                s_current        = s_current,
                                candidate_output = cs,
                                r                = r,
                                cutoff_opt       = 0.99)
    Reduce_obstacle <- unique(c(Reduce_obstacle, rcs))
    if (length(Reduce_obstacle) > 0) {
      obstacle_df$prob[Reduce_obstacle] <- 1
    }
    
    n_cand <- length(Set_state)
    
    # early stop condition
    if(n_cand == 1 && Set_state[1] == g){
      Record_path  <- c(Record_path, Path_record[[1]])
      Total_length <- Total_length + Cost_imm[1]
      s_current    <- g
      break
    }
    
    # --- 4(c) Initialize & Pre-compute info gain --------------------------------
    # initialise value states for new candidates
    new_states <- Set_state[!(Set_state %in% V_states)]
    if (length(new_states) > 0) {
      V_states   <- c(V_states,   new_states)
      V_function <- c(V_function, rep(0, length(new_states)))
      N_states   <- c(N_states,   rep(1, length(new_states)))
    }
    
    # info gain
    G_values <- compute_info_gain_approx(neighbors_list  = neighbors_list,
                                         Record_obstacle = Record_obstacle,
                                         cov_prior       = cov_prior,
                                         noise_variance  = noise_variance,
                                         alpha           = alpha_ig,
                                         I_sum           = I_sum)
    
    # --- 4(d) Offline loop -----------------------------------------
    if(idx_step == 1){
      ChangeValue_state <- Set_state[-n_cand] # remove g (no need to track)
      ChangeValue       <- rep(Inf, length(ChangeValue_state))
      extra_count       <- 0
      extra_triggered   <- FALSE
      n_offline         <- 0
      
      repeat{
        # sample environment
        sample_obs <- prediction_sample(obstacle_df       = obstacle_df, 
                                        pos_mean_log_odds = mean_prior, 
                                        pos_cov           = cov_prior)
        
        # weighted start state sampling
        nc <- N_states[match(Set_state[-n_cand], V_states)]
        nc <- ifelse(is.na(nc), 1, nc)
        w  <- 1/(nc + 0.1)
        w  <- w/sum(w)
        idx_start <- sample(n_cand-1, 1, prob=w)
        
        # simulate trajectory
        traj <- simulate_trajectory(s_start         = Set_state[idx_start],
                                    g               = g,
                                    obs_start       = Set_obstacle[idx_start],
                                    G_original      = G_original,
                                    n_vertices      = n_vertices,
                                    obstacle_df     = sample_obs,
                                    intersect_mat   = intersect_mat,
                                    edge_df         = edge_df,
                                    vertice_list    = vertice_list,
                                    V_states        = V_states,
                                    V_function      = V_function,
                                    r               = r,
                                    epsilon         = epsilon,
                                    method_decision = method_decision,
                                    G_values        = G_values)
        states_exp <- traj$trajectory_states
        costs_exp  <- rev(cumsum(rev(traj$trajectory_costs)))
        
        # incremental value update
        for (i in seq_len(length(states_exp) - 1)){
          sv   <- states_exp[i]
          cost <- costs_exp[i]
          idx  <- match(sv, V_states)
          
          if(!is.na(idx)){
            N_states[idx] <- N_states[idx] + 1
            a             <- 1 / N_states[idx]
            new_val       <- (1 - a) * V_function[idx] + a * cost
            if(sv %in% ChangeValue_state){
              ci <- match(sv, ChangeValue_state)
              ChangeValue[ci] <- abs(V_function[idx] - new_val) / (abs(V_function[idx]) + 1e-8)
            }
            V_function[idx] <- new_val
          } else{
            V_states   <- c(V_states, sv)
            N_states   <- c(N_states, 1)
            V_function <- c(V_function, cost)
            if (sv %in% ChangeValue_state) {
              ChangeValue[match(sv, ChangeValue_state)] <- cost
            }
          }
        }
        n_offline <- n_offline + 1
        
        # convergence check
        if (max(ChangeValue) <= conv_tol) {
          extra_triggered <- TRUE
          extra_count     <- extra_count + 1
          if (extra_count >= conv_patience) {
            break
          }
        } else if (extra_triggered) {
          extra_triggered <- FALSE
          extra_count     <- 0
        }
        if (n_offline >= N_itr_max) {
          break
        }
      }
    }
    
    # --- 4(e) online loop ------------------------------------------
    if(idx_step > 1){
      
      # pre-allocate rollout storage
      init_cap    <- n_rollouts * (n_cand - 1) * 15
      state_buf   <- integer(init_cap)
      cost_buf    <- numeric(init_cap)
      cap         <- init_cap
      ptr         <- 1
      
      # loop
      for (roll in seq_len(n_rollouts)){
        sample_obs <- prediction_sample(obstacle_df       = obstacle_df, 
                                        pos_mean_log_odds = mean_prior, 
                                        pos_cov           = cov_prior)
        results <- lapply(seq_len((n_cand - 1)), function(si){
          traj <- simulate_trajectory(s_start         = Set_state[si],
                                      g               = g,
                                      obs_start       = Set_obstacle[si],
                                      G_original      = G_original,
                                      n_vertices      = n_vertices,
                                      obstacle_df     = sample_obs,
                                      intersect_mat   = intersect_mat,
                                      edge_df         = edge_df,
                                      vertice_list    = vertice_list,
                                      V_states        = V_states,
                                      V_function      = V_function,
                                      r               = r,
                                      epsilon         = epsilon,
                                      method_decision = "greedy",
                                      G_values        = G_values)
          ns <- length(traj$trajectory_states) - 1
          if(ns < 1){
            return(NULL)
          }
          list(states = traj$trajectory_states[seq_len(ns)],
               costs  = rev(cumsum(rev(traj$trajectory_costs)))[seq_len(ns)],
               n      = ns)
        })
        
        for(res in Filter(Negate(is.null), results)){
          need <- ptr +  res$n - 1
          if (need > cap){
            cap <- max(need, cap*2)
            length(state_buf) <- cap
            length(cost_buf) <- cap
          }
          idx_range <- ptr : (ptr + res$n - 1)
          state_buf[idx_range] <- res$states
          cost_buf[idx_range] <- res$costs
          ptr <- ptr + res$n
        }
      }
      
      # aggregate rollout results
      final_n   <- ptr - 1
      state_buf <- state_buf[seq_len(final_n)]
      cost_buf  <- cost_buf[seq_len(final_n)]
      
      dt  <- data.table(s = state_buf, c = cost_buf)
      agg <- dt[, .(mean_cost = mean(c), count = .N), by = s]
      
      
      # --- 4(f) baseline adjust --------------------------------------
      existing    <- match(agg$s, V_states)
      is_exist    <- !is.na(existing)
      
      if (any(is_exist)) {
        ei  <- existing[is_exist]
        sub <- agg[is_exist,]
        w   <- sub$count / (N_states[ei] + sub$count)
        V_function[ei] <- (1 - w) * V_function[ei] + w * sub$mean_cost
        N_states[ei]   <- N_states[ei] + sub$count
      }
      
      if (any(!is_exist)) {
        new_rows   <- agg[!is_exist,]
        V_states   <- c(V_states,   new_rows$s)
        V_function <- c(V_function, new_rows$mean_cost)
        N_states   <- c(N_states,   new_rows$count)
      }
      
    }
    # --- 4(g) Execute & record -------------------------------------
    # future value per candidate
    Cost_future <- vapply(Set_state, function(sv){
      idx <- match(sv, V_states)
      ifelse(!is.na(idx), V_function[idx], 0)
    }, numeric(1))
    
    Cost_future[n_cand] <- 0
    
    # total expected cost per candidate
    Cost_imm_full <- Cost_imm + c(obstacle_df$disambiguation_cost[Set_obstacle[-n_cand]], 0)
    Cost_expected <- Cost_imm_full + Cost_future
    
    # scaling factor for info gain
    alpha_local <- if (n_cand > 2) {
      sd(Cost_expected[Cost_expected <= median(Cost_expected)]) / div_IG
    } else{
      0
    }
    
    # exact info gain per candidate
    gain_info  <- numeric(n_cand)
    I_sum_list <- numeric(n_cand)
    for (i in seq_len(n_cand - 1)) {
      ig <- compute_info_gain_exact(state           = Set_state[i],
                                    vertice_list    = vertice_list,
                                    obstacle_df     = obstacle_df,
                                    Record_obstacle = Record_obstacle,
                                    cov_prior       = cov_prior,
                                    noise_variance  = noise_variance,
                                    sensing_range   = sensing_range,
                                    I_sum           = I_sum,
                                    alpha           = alpha_local)
      gain_info[i]  <- ig$gain_info
      I_sum_list[i] <- ig$I_sum
    }
    
    # decision: minimise expected cost minus info gain
    cost_decision <- Cost_expected - gain_info
    
    # if goal is optimal — go directly
    if (which.min(cost_decision) == n_cand) {
      Record_path  <- c(Record_path, Path_record[[n_cand]])
      Total_length <- Total_length + Cost_imm_full[n_cand]
      s_current    <- g
      break
    }
    
    # among tied minima prefer lowest obstacle probability
    min_val <- min(cost_decision)
    tied    <- which(cost_decision == min_val)
    sel     <- tied[which.min(obstacle_df$prob[Set_obstacle[tied]])]
    
    s_next        <- Set_state[sel]
    obstacle_next <- Set_obstacle[sel]
    I_sum         <- I_sum_list[sel]
    
    # record executed step
    Record_path  <- c(Record_path, Path_record[[sel]])
    Total_length <- Total_length + Cost_imm_full[sel]
    
    if (s_next != g) {
      Total_cost              <- Total_cost +
        obstacle_df$disambiguation_cost[obstacle_next]
      Record_obstacle[obstacle_next] <- TRUE
      New_obstacle            <- obstacle_next
    } else {
      New_obstacle <- integer(0)
    }
    
    s_current <- s_next
    idx_step  <- idx_step + 1
    alpha_ig <- alpha_local
  }
  
  # OUTPUT
  output <- list(Record_path         = Record_path,
                 Path_cost           = Total_length, # already include disambiguation cost
                 Disambiguation_cost = Total_cost,
                 Record_obstacle     = Record_obstacle)
  return(output)
}
