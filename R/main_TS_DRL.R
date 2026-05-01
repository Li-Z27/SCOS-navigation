# =============================================================================
# main_TS_DRL.R
# Two-stage policy (MC)
# =============================================================================
Policy_TS_DRL <- function(obstacle_df, x, y, r, s, g, status, 
                          sensing_range, sense_noise,
                          noise_var_scale, l, sigma_f, replace_condi,
                          alpha_ig, div_IG){
  # INPUT
  # --- CAT I: Graph-related ----
  # obstacle_df: obstacle information
  # x_coor, y_coor
  # disambiguation_cost
  # prob
  # status
  # x, y: grid dimension
  # r: obstacle radius
  # s, g: start and goal vertex
  # status: ground truth status vector for all obstacles
  
  # --- CAT II: Agent-related ---
  # sensing_range: the range in which the sensor can get prob reading
  # sense_noise: lambda value (higher lambda, higher precision level)
  
  # --- CAT III: Belief-related ---
  # noise_var_scale: scale used based on variance from delta method
  # l: length scale
  # sigma_f: marginal standard deviation
  # replace_condi: support refinement threshold delta
  
  # --- CAT IV: InformationGain-related ---
  # alpha_ig: initial information gain scaling factor
  # div_IG: scaling factor adjusted based on observed values
  
  # --- STEP 0: Hyperparameters -------------------------------------
  n_rollouts    <- 30     # online rollouts per decision step (adjusted based on problem size)
  N_itr_max     <- 2000   # max offline MC iterations
  conv_tol      <- 0.005  # convergence tolerance
  conv_patience <- 10     # consecutive under-tolerance iters before stopping
  
  # --- STEP 1 Graph Construction -----------------------------------
  obstacle_df$status <- status
  obstacle_df        <- subset(obstacle_df, select = -prob_true)
  
  vertice_list  <- Lattice_Vertices(x, y)
  G_original    <- Graph_Discretized(x, y)
  intersect_mat <- Intersect_Info(G_original, x, y, obstacle_df, r)
  
  edge_df    <- get.data.frame(G_original, what = "edges")
  n_vertices <- vcount(G_original)
  
  # --- STEP 2 Belief Initialization --------------------------------
  noise_var      <- (64 / ((4+sense_noise) * (4-sense_noise) * 9))*noise_var_scale
  mean_prior     <- rep(0, nrow(obstacle_df))
  cov_prior      <- create_prior_cov(obstacle_df, sigma_f, l)
  noise_variance <- rep(noise_var, nrow(obstacle_df))
  # neighbors list for approximation info gain
  neighbors_list <- get_neighbors_list(obstacle_df, sensing_range + r)
  
  # --- STEP 3 Tracking Initialization ------------------------------
  s_current       <- s
  Record_obstacle <- rep(FALSE, nrow(obstacle_df)) # TRUE once disambiguated
  Reduce_obstacle <- integer(0) # obstacles confirmed as always-blocking/unnecessary
  New_obstacle    <- integer(0) # obstacle just disambiguated at recent step
  
  V_states    <- integer(0) # visited states with estimated distribution
  V_function  <- list() # support grid
  V_prob      <- list() # dirichlet parameters
  V_indicator <- list()  # indicator for support values - TRUE: experienced during simulation
  N_states    <- integer(0) # visit counts per state 
  
  I_sum        <- 0 # cumulative information
  Total_length <- 0 # path length (Euclidean)
  Total_cost   <- 0 # disambiguation cost accumulated
  Record_path  <- integer(0) # ordered vertex indices of executed path

  idx_step     <- 1
  
  # --- STEP 4 Main Loop --------------------------------------------
  while (s_current != g) {
    # stop once reach the goal vertex
    
    # --- 4(a) Sense & update belief --------------------------------
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
    # value convgernce track prepare
    n_cand <- length(Set_state)
    ChangeValue_state <- Set_state[-n_cand] # remove g (no need to track)
    
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
    
    # early stop condition
    if(n_cand == 1 && Set_state[1] == g){
      Record_path  <- c(Record_path, Path_record[[1]])
      Total_length <- Total_length + Cost_imm[1]
      s_current    <- g
      break
    }
    
    # --- 4(c) Initialize & Pre-compute info gain -------------------
    # initialise value states for new candidates
    new_states <- Set_state[!(Set_state %in% V_states)]
    if (length(new_states) > 0) {
      new_inits <- lapply(new_states, function(sv) {
        init_state_value(sv, g, obstacle_df, intersect_mat,
                         G_original, edge_df, replace_condi)
      })
      
      V_states    <- c(V_states,   new_states)
      V_function  <- c(V_function, lapply(new_inits, `[[`, "V_function"))
      V_prob      <- c(V_prob,     lapply(new_inits, `[[`, "V_prob"))
      V_indicator <- c(V_indicator,lapply(new_inits, `[[`, "V_indicator"))
      N_states    <- c(N_states,   rep(0, length(new_states)))
    }
    
    # info gain
    G_values <- compute_info_gain_approx(neighbors_list  = neighbors_list,
                                         Record_obstacle = Record_obstacle,
                                         cov_prior       = cov_prior,
                                         noise_variance  = noise_variance,
                                         alpha           = alpha_ig,
                                         I_sum           = I_sum)
    
    # --- 4(d) Offline loop -----------------------------------------
    if (idx_step == 1) {
      ChangeValue       <- rep(Inf, length(ChangeValue_state))
      ChangeSD          <- rep(Inf, length(ChangeValue_state))
      MeanValue         <- rep(0,   length(ChangeValue_state))
      SDValue           <- rep(0,   length(ChangeValue_state))
      
      extra_count       <- 0
      extra_triggered   <- FALSE
      n_offline         <- 0
      
      repeat{
        # sample environment
        sample_obs <- prediction_sample(obstacle_df       = obstacle_df, 
                                        pos_mean_log_odds = mean_prior, 
                                        pos_cov           = cov_prior)
        # weighted start state sampling
        nc  <- N_states[match(Set_state[-n_cand], V_states)]
        nc  <- ifelse(is.na(nc), 1, nc)
        w   <- 1 / (nc + 0.1)
        w <- w / sum(w)
        idx_start <- sample(n_cand - 1, 1, prob = w)
        
        # simulate with Bellman update and support refinement
        traj <- simulate_trajectory_DRL(s_start       = Set_state[idx_start],
                                        g             = g,
                                        obs_start     = Set_obstacle[idx_start],
                                        obstacle_df   = sample_obs,
                                        intersect_mat = intersect_mat,
                                        G_original    = G_original,
                                        edge_df       = edge_df,
                                        vertice_list  = vertice_list,
                                        V_states      = V_states,
                                        V_function    = V_function,
                                        V_prob        = V_prob,
                                        V_indicator   = V_indicator,
                                        N_states      = N_states,
                                        replace_condi = replace_condi,
                                        G_values      = G_values)
        V_states    <- traj$V_states
        V_function  <- traj$V_function
        V_prob      <- traj$V_prob
        V_indicator <- traj$V_indicator
        N_states    <- traj$N_states
        
        # update expected value and change value track
        for (ci in seq_along(ChangeValue_state)) {
          sv  <- ChangeValue_state[ci]
          idx <- match(sv, V_states)
          if (is.na(idx)){
            next
          } else {
            p_temp          <- V_prob[[idx]] / sum(V_prob[[idx]])
            cum_p           <- cumsum(p_temp)
            new_mean        <- sum(p_temp * V_function[[idx]])
            new_sd   <- sqrt(sum(p_temp * (V_function[[idx]] - new_mean)^2))
            ChangeValue[ci] <- max(
              abs(new_mean - MeanValue[ci]) / (abs(MeanValue[ci]) + 1e-8)
            )
            ChangeSD[ci] <- abs(new_sd   - SDValue[ci]) / (abs(SDValue[ci])   + 1e-8)
            
            MeanValue[ci] <- new_mean
            SDValue[ci]   <- new_sd
          }
        }
        n_offline <- n_offline + 1
        
        # convergence track
        if (max(ChangeValue) <= conv_tol && 
            max(ChangeSD) <= conv_tol && 
            n_offline >= nrow(obstacle_df)*15) {
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
    
    # --- 4(e) Online loop ------------------------------------------
    if (idx_step > 1){
      all_states <- list()
      all_costs  <- list()
      ptr        <- 1
      
      ChangeValue_online <- rep(Inf, length(ChangeValue_state))
      ChangeSD_online  <- rep(Inf, length(ChangeValue_state))
      MeanValue_online   <- vapply(ChangeValue_state, function(sv) {
        idx <- match(sv, V_states)
        p_temp <- V_prob[[idx]] / sum(V_prob[[idx]])
        sum(p_temp * V_function[[idx]])
      }, numeric(1))
      SDValue_online     <- vapply(ChangeValue_state, function(sv) {
        idx    <- match(sv, V_states)
        p_temp <- V_prob[[idx]] / sum(V_prob[[idx]])
        new_mean <- sum(p_temp * V_function[[idx]])
        sqrt(sum(p_temp * (V_function[[idx]] - new_mean)^2))
      }, numeric(1))
      
      
      
      n_online               <- 0
      max_online             <- 1000
      
      repeat {
        sample_obs <- prediction_sample(obstacle_df       = obstacle_df, 
                                        pos_mean_log_odds = mean_prior, 
                                        pos_cov           = cov_prior)
        
        # weighted sampling — prefer under-visited candidates
        nc  <- N_states[match(Set_state[-n_cand], V_states)]
        nc  <- ifelse(is.na(nc), 1, nc)
        w   <- 1 / (nc + 0.1) 
        w   <- w / sum(w)
        si  <- sample(n_cand - 1, 1, prob = w)
        
        traj <- simulate_trajectory_DRL_rollout(s_start       = Set_state[si],
                                                g             = g,
                                                obs_start     = Set_obstacle[si],
                                                obstacle_df   = sample_obs,
                                                intersect_mat = intersect_mat,
                                                G_original    = G_original,
                                                edge_df       = edge_df,
                                                vertice_list  = vertice_list,
                                                V_states      = V_states,
                                                V_function    = V_function,
                                                V_prob        = V_prob,
                                                V_indicator   = V_indicator,
                                                N_states      = N_states,
                                                replace_condi = replace_condi,
                                                G_values      = G_values)
        V_states    <- traj$V_states
        V_function  <- traj$V_function
        V_prob      <- traj$V_prob
        V_indicator <- traj$V_indicator
        N_states    <- traj$N_states
        # collect actual cumulative costs for support refinement
        traj_s    <- traj$trajectory_states
        traj_c    <- traj$trajectory_costs
        n_steps   <- length(traj_c)
        
        if (n_steps > 0) {
          costs_cum          <- rev(cumsum(rev(traj_c)))
          all_states[[ptr]]  <- traj_s[seq_len(n_steps)]
          all_costs[[ptr]]   <- costs_cum
          ptr                <- ptr + 1
        }
        
        # convergence check — same mechanism as offline
        for (ci in seq_along(ChangeValue_state)) {
          sv  <- ChangeValue_state[ci]
          idx <- match(sv, V_states)
          if (is.na(idx)) next
          p_temp                 <- V_prob[[idx]] / sum(V_prob[[idx]])
          cum_p                  <- cumsum(p_temp)
          new_mean               <- sum(p_temp * V_function[[idx]])
          new_sd   <- sqrt(sum(p_temp * (V_function[[idx]] - new_mean)^2))
          
          ChangeValue_online[ci] <- max(
            abs(new_mean - MeanValue_online[ci]) / (abs(MeanValue_online[ci]) + 1e-8)
          )
          ChangeSD_online[ci] <- abs(new_sd - SDValue_online[ci]) / (abs(SDValue_online[ci]) + 1e-8)
          MeanValue_online[ci]   <- new_mean
          SDValue_online[ci]     <- new_sd
          
        }
        
        n_online <- n_online + 1
        
        if (max(ChangeValue_online) <= conv_tol && 
            max(ChangeSD_online) <= conv_tol && 
            n_online >= max(n_cand*n_rollouts, nrow(obstacle_df)*10)) break
        if (n_online >= max_online)             break
      }
      # --- 4(f) Support refinement after all online rollouts ---------
      # flatten all trajectory observations into two vectors
      all_s <- unlist(all_states)
      all_c <- unlist(all_costs)
      
      # get unique states that need refinement
      unique_states <- unique(all_s)
      
      # process each unique state once — applying all its observations sequentially
      for (sv in unique_states) {
        idx <- match(sv, V_states)
        if (is.na(idx)) next
        
        # get all observed costs for this state across all rollouts
        obs_costs <- all_c[all_s == sv]
        
        # apply all observations sequentially
        for (val_obs in obs_costs) {
          upd                <- update_support(V_function[[idx]], V_prob[[idx]],
                                               V_indicator[[idx]], val_obs,
                                               replace_condi)
          V_function[[idx]]  <- upd$V_function
          V_prob[[idx]]      <- upd$V_prob
          V_indicator[[idx]] <- upd$V_indicator
        }
      }
    }
    
    # --- 4(g) Execute & record -------------------------------------
    # probability-weighted mean per candidate
    Cost_future <- vapply(Set_state, function(sv) {
      idx    <- match(sv, V_states)
      p_temp <- V_prob[[idx]] / sum(V_prob[[idx]])
      sum(p_temp * V_function[[idx]])
    }, numeric(1))
    Cost_future[n_cand] <- 0
    
    Cost_imm_full <- Cost_imm +
      c(obstacle_df$disambiguation_cost[Set_obstacle[-n_cand]], 0)
    Cost_expected <- Cost_imm_full + Cost_future
    
    # scaling factor for info gain
    alpha_local <- if (n_cand > 2) {
      sd(Cost_expected[Cost_expected <= median(Cost_expected)]) / div_IG
    } else {
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
    
    # tie-breaking by lowest obstacle probability
    min_val <- min(cost_decision)
    tied    <- which(cost_decision == min_val)
    sel     <- tied[which.min(obstacle_df$prob[Set_obstacle[tied]])]
    
    s_next        <- Set_state[sel]
    obstacle_next <- Set_obstacle[sel]
    I_sum         <- I_sum_list[sel]
    
    # record executed step
    Record_path  <- c(Record_path, Path_record[[sel]])
    Total_length <- Total_length + Cost_imm_full[sel]
    
    if (s_next != g){
      Total_cost                     <- Total_cost +
        obstacle_df$disambiguation_cost[obstacle_next]
      Record_obstacle[obstacle_next] <- TRUE
      New_obstacle                   <- obstacle_next
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
